#!/usr/bin/env python3

"""
    Utilities for plotting individual lightcurves.
    
    Command line arguments:
    --lc-path: The path containing .h5 lc files.
    --tic-gaia-fname: The output file to write tic to gaia mappings.
        --tic-gaia-file: if we already have it, specify the file path
    --toi-gaia-period-fname: The output file to write toi, gaia, period mappings.
    --lc-catalog: The input file names of all the GAIA catalogs that are created up to magfit.

    1. We have the LCs, want to extract their Gaia ids:
        lcs_to_gaia_ids(lc_path) -> gaia_ids
        It goes over the LCs, and extract their Gaia Ids from their names,
        and return the unique set of them.
    
    (*). We need the FOV and mag range of the LCs to be used in step 3.*** 
        find_fov_of_lcs(lc_catalog) -> fov_mag_range = {ra, dec, mag_{min, max}}
        It goes over the LC catalog and find the range (min, max) of ra, dec, and mag.

    3. To extract TICs from Gaia ids, we have three options:
        
        3.*. If we already have the tic_gaia mapping from a previous run,
        we can just read the related file and use it.

        query_vizier_n_times(gaia_ids, tic_gaia_fname) -> tic_gaia
        3.**. For each Gaia id, query vizier to find its TIC id
        We then write the tic_gaia mapping to a file for later use cases,
        and also to avoid querying vizier again if we already have it.
        (we may use this one for a small number of LCs, like < 1000)
        
        query_vizier_once(gaia_ids, tic_gaia_fname, fov_mag_range) -> tic_gaia
        3.***. Query vizier once based on the FOV of the LCs (ra, dec, mag range)
        that we got in (*) to find all TICs in that area.
        Then, we check the result by going over each row and check if the Gaia id
        in that row is in our list of Gaia ids from the LCs.
        If so, we save the TIC id and Gaia id in a dictionary.
        (we may use this one for a large number of LCs, like > 1000)

    4. Extract TOI IDs and their periods for the FOV of our LCs. 
        query_toi_in_fov(tic_gaia, cmdline_args.toi_gaia_period_fname, fov_mag_range) -> None
        Finally, it queries TOI table in exoplanet archive in the range in the FOV of the LCs
        and check if any of those TOI have a TIC id in our list. If so, it extracts the period.
        Also, it writes: toi_gaia_period.txt
    -------------------------------------------------
        
    Check this link for all tables that you may need:
    https://tapvizier.cds.unistra.fr/adql/
    more info in discord:
    https://discord.com/channels/581176949184135229/1291434287614398494/1356713434086903970

    Example usage:
    1. When you don't have tic_gaia.txt file:
            python gaia_to_toi.py --lc-path /path/to/lcs/ --lc-catalog /path/to/lc_catalog.fits
        2. When you already have tic_gaia.txt file:
            python gaia_to_toi.py --lc-path /path/to/lcs/ --lc-catalog /path/to/lc_catalog.fits --tic-gaia-file /path/to/tic_gaia.txt
"""


import os
import os.path
import glob
import logging
import tenacity
import pandas as pd
import pyvo
import astropy.units as u
from astropy.table import Table, vstack
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia
from configargparse import ArgumentParser, RawDescriptionHelpFormatter

# Configure logging for retry attempts
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

def parse_command_line():
    """Return the parsed command line arguments."""

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--lc-path',
        help='The path containing .h5 lc files.'
    )
    parser.add_argument(
        '--lc-catalogs',
        nargs='+',
        help='The input file names of all the GAIA catalogs.'
    )
    parser.add_argument(
        '--tic-gaia-fname',
        default='tic_gaia.txt',
        help='The output file to write tic to gaia mappings.'
    )
    parser.add_argument(
        '--tic-gaia-file',
        default=False,
        help='if we already have it, specify the file path'
    )
    parser.add_argument(
        '--toi-gaia-period-fname',
        default='toi_gaia_period.txt',
        help='The output file to write toi, gaia, period mappings.'
    )
    parser.add_argument(
        '--skip-gaia-variables',
        default=False,
        action='store_true',
        help='Skip querying Gaia variable classifications.'
    )

    # parser.add_argument(
    #     '--tables-to-query',
    #     default=['"IV/39/tic82"', 'toi'],
    #     help='The tables to query respectively.'
    # )


    return parser.parse_args()

# ------------------------------------------------------------------#
# - specialized for Vizier TAP service (ORIG-kept for reference) -- #
# ------------------------------------------------------------------#

@tenacity.retry(
    wait=tenacity.wait_exponential(multiplier=1, min=2, max=60),
    stop=tenacity.stop_after_attempt(10),
    reraise=True,
    before_sleep=tenacity.before_sleep_log(__import__('logging').getLogger(__name__), __import__('logging').WARNING)
)
def tap_vizier_query(url='http://tapvizier.u-strasbg.fr/TAPVizieR/tap/',
                     headers='*',
                     table_database='"J/ApJS/258/16/tess-ebs"', 
                     # TESS Eclipsing Binary stars. I. Sectors 1-26 (Prsa+, 2022) # 4.5 thousands rows
                     constraints=['1=1']): # a dummy constraint instead of None and later handle it
    """
        The TAP query to run. All queries must have valid ADQL, which takes the form:
        SELECT <column list> FROM <table> WHERE <constraints>

        Args:
            - url: The url to use for TAP query
            - headers (list) : The table column headers to use from given table database
            - table_database: The table database name to use (usually in quotations for vizier)
            - constraints (str list): The given constraints to use for TAP query

        Returns:
            - results: The TAP query results as a VOTable
        """

    service = pyvo.dal.TAPService(url)

    headers = ", ".join(item for item in headers)
    constraints = " WHERE " + " AND ".join(constraints)
    query = "SELECT " + headers + " FROM " + table_database + constraints

    results = service.search(query)

    return results


# ------------------------------------------------------------------#
# ------------- NEW ASTROQUERY-based functions ---------------------#
# ---------- Query a single Gaia ID in TIC catalog -----------------#
# ------------------------------------------------------------------#

@tenacity.retry(
    wait=tenacity.wait_exponential(multiplier=1, min=2, max=60),
    stop=tenacity.stop_after_attempt(10),
    reraise=True,
    before_sleep=tenacity.before_sleep_log(__import__('logging').getLogger(__name__), __import__('logging').WARNING)
)
def query_tic_by_gaia(gaia_id):
    """
    Query TIC 8.2 catalog for a specific Gaia ID.
    
    Args:
        - gaia_id: The Gaia ID to search for
        
    Returns:
        - result: Query result table or empty table if no match found
    """
    viz = Vizier(columns=['TIC', 'GAIA', 'Vmag', 'RAJ2000', 'DEJ2000'], timeout=300)
    viz.ROW_LIMIT = 1  # We expect at most one match per Gaia ID
    
    result = viz.query_constraints(catalog='IV/39/tic82', GAIA=gaia_id)
    return result[0] if len(result) > 0 else None

# ------------------------------------------------------------------#
# ------ Query TIC catalog in a spatial region (FOV) ---------------#
# ------------------------------------------------------------------#
@tenacity.retry(
    wait=tenacity.wait_exponential(multiplier=1, min=2, max=60),
    stop=tenacity.stop_after_attempt(10),
    reraise=True,
    before_sleep=tenacity.before_sleep_log(__import__('logging').getLogger(__name__), __import__('logging').WARNING)
)
def query_tic_in_region(fov_mag_range):
    """
    Query TIC 8.2 catalog in a spatial region defined by RA, Dec, and magnitude range.
    
    Args:
        - fov_mag_range: Dictionary with 'ra_min', 'ra_max', 'dec_min', 'dec_max', 'mag_max'
        
    Returns:
        - result: Query result table
    """
    viz = Vizier(
        columns=['TIC', 'GAIA', 'Vmag', 'RAJ2000', 'DEJ2000'],
        column_filters={"Vmag":f"< {fov_mag_range['mag_max']}"},
        # column_filters={"Vmag":"< 8.0"},
        # timeout=300
    )
    print(f'viz: {viz}')
    viz.ROW_LIMIT = -1  # No limit for region queries
    Vizier.SERVER = 'https://vizier.cds.unistra.fr'
    
    # Define the center and radius of the search region (approximate bounding box as circular region)
    ra_center = (fov_mag_range['ra_min'] + fov_mag_range['ra_max']) / 2
    dec_center = (fov_mag_range['dec_min'] + fov_mag_range['dec_max']) / 2
    
    # Approximate radius in degrees
    ra_width = abs(fov_mag_range['ra_max'] - fov_mag_range['ra_min'])
    dec_width = abs(fov_mag_range['dec_max'] - fov_mag_range['dec_min'])
    radius = max(ra_width, dec_width) / 2 + 0.5  # Add buffer
    
    coord = SkyCoord(ra=ra_center, dec=dec_center, unit=(u.deg, u.deg))
    
    # Vizier's query_region uses a circular region, which may include extra objects
    # We'll filter by magnitude and do bounds checking in the calling function
    print(f"Querying TIC catalog in region centered at (RA: {ra_center:.2f}, Dec: {dec_center:.2f}) with radius {radius:.2f} deg and Vmag < {fov_mag_range['mag_max']:.2f}")
    result = viz.query_region(
        coord,
        radius=radius * u.deg,
        catalog='IV/39/tic82'
    )
    with open('query_result.txt', 'w') as file:
        for row in result[0]:
            file.write(f"{row['GAIA']}\n")
    # Use len() instead of truthiness to avoid ambiguous Quantity comparison error
    return result[0] if len(result) > 0 else None


# ------------------------------------------------------------------#
# ------ Extract all GAIA IDs from light curves --------------------#
# ------------------------------------------------------------------#

def lcs_to_gaia_ids(lc_path):

    # Gather all lcs full paths:
    lcs = sorted(glob.glob(os.path.join(lc_path, '*.h5')))
    # lcs = ['/path/to/GDR3_1316714794020532096.h5', ...]

    # from lc full paths to the int(Gaia id):
    gaia_ids = {
        int(os.path.basename(lc).split('_')[1].split('.')[0])
        for lc in lcs
    }
    # gaia_ids = {1316714794020532096, ...}
    with open('gaia_ids.txt', 'w') as file:
        for gaia_id in gaia_ids:
            file.write(f"{gaia_id}\n")
    print(f'Found {len(gaia_ids)} unique gaia ids. Now querying TIC\n')

    return gaia_ids


# ------------------------------------------------------------------#
# --------- Extract all GAIA variables from LCs --------------------#
# ------------------------------------------------------------------#

def gaia_ids_to_variables(gaia_ids, batch_size=100):

    query_varialbles = """
    SELECT v.source_id,
           v.best_class_name,
           v.classification_score
    FROM gaiadr3.vari_classifier_result AS v
    JOIN user_table AS u
    ON v.source_id = u.source_id
    """

    gaia_ids_list = sorted(gaia_ids)
    results_batches = []

    for start in range(0, len(gaia_ids_list), batch_size):
        batch = gaia_ids_list[start:start + batch_size]
        df_gaia = pd.DataFrame({'source_id': batch}, dtype='int64')
        gaia_table = Table.from_pandas(df_gaia)
        upload_table_name = f"user_table_{start}_{len(batch)}"

        try:
            job = Gaia.launch_job_async(
                query=query_varialbles,
                upload_resource=gaia_table,
                upload_table_name=upload_table_name
            )
            batch_results = job.get_results()
        except Exception as exc:
            logging.warning(
                'Gaia variables query failed for batch %d-%d: %s',
                start,
                start + len(batch) - 1,
                exc
            )
            continue

        if batch_results is not None and len(batch_results) > 0:
            results_batches.append(batch_results)

    if not results_batches:
        print('No Gaia variable results returned.')
        return None

    results = vstack(results_batches)

    with open('gaia_variables.txt', 'w') as file:
        for row in results:
            file.write(
                f"{row['source_id']}, {row['best_class_name']}, {row['classification_score']}\n"
            )

    return None

# ------------------------------------------------------------------#
# ----------- Find the FOV of the LCs (catalog) --------------------#
# ------------------------------------------------------------------#

def find_fov_of_lcs(lc_catalogs):

    ra_min_all, ra_max_all = 1e10, -1e10
    dec_min_all, dec_max_all = 1e10, -1e10
    mag_min_all, mag_max_all = 1e10, -1e10

    for lc_catalog in lc_catalogs:
        with fits.open(lc_catalog) as hdul:
            lc_data = hdul[1].data

            ra_min, ra_max = lc_data['ra'].min(), lc_data['ra'].max()
            ra_min_all = min(ra_min_all, ra_min)
            ra_max_all = max(ra_max_all, ra_max)

            dec_min, dec_max = lc_data['dec'].min(), lc_data['dec'].max()
            dec_min_all = min(dec_min_all, dec_min)
            dec_max_all = max(dec_max_all, dec_max)
        
            mag_min, mag_max = lc_data['phot_g_mean_mag'].min(), lc_data['phot_g_mean_mag'].max()
            mag_min_all = min(mag_min_all, mag_min)
            mag_max_all = max(mag_max_all, mag_max)

    print(f'Ranges in lc catalogs:::\n'
            f'ra: ({ra_min_all:.2f}, {ra_max_all:.2f})\n'
            f'dec: ({dec_min_all:.2f}, {dec_max_all:.2f})\n'
            f'mag: ({mag_min_all:.2f}, {mag_max_all:.2f})\n'
    )

    return {'ra_min': ra_min_all,
            'ra_max': ra_max_all,
            'dec_min': dec_min_all,
            'dec_max': dec_max_all,
            'mag_min': mag_min_all,
            'mag_max': mag_max_all}


# ------------------------------------------------------------------#
# --- For each GAIA ID, query vizier to find its TIC -------------- #
# -------- (not efficient sometimes) ------------------------------ #
# --- ORIGINAL TAP VERSION (kept for reference) ------------------- #
# ------------------------------------------------------------------#

def query_vizier_n_times_tap(gaia_ids, tic_gaia_fname):

    tic_gaia = {}

    n = 0
    for gaia_id in gaia_ids:
        n += 1
        if n % 10 == 0: print(f'TIC Query result: passed {n}')

        query_gaia_to_tic = tap_vizier_query(
            url='http://tapvizier.u-strasbg.fr/TAPVizieR/tap/',
            headers='*', # or just TIC column?
            table_database='"IV/39/tic82"', 
            # TESS Input Catalog version 8.2 (TIC v8.2) (Paegert+, 2021) 1.7 Bilion rows
            constraints=['GAIA=' + str(gaia_id)]
        )
        if len(query_gaia_to_tic) > 0:
            tic_id = query_gaia_to_tic['TIC'].data[0]
            tic_gaia[tic_id] = gaia_id

    print(f'Found {len(tic_gaia)} unique TIC ids. Now querying TOI in FOV of LCs in the LC catalog\n')

    with open(tic_gaia_fname, 'w') as file:
        for tic_id, gaia_id in tic_gaia.items():
            file.write(f"{tic_id}, {gaia_id}\n")
    
    return tic_gaia


# ------------------------------------------------------------------#
# -- For each GAIA ID, query vizier to find its TIC (ASTROQUERY) -- #
# ------------------------------------------------------------------#

def query_vizier_n_times(gaia_ids, tic_gaia_fname):

    tic_gaia = {}

    n = 0
    for gaia_id in gaia_ids:
        n += 1
        if n % 10 == 0: print(f'TIC Query result: passed {n}')

        result = query_tic_by_gaia(gaia_id)
        
        if result is not None and len(result) > 0:
            tic_id = result['TIC'][0]
            tic_gaia[tic_id] = gaia_id

    print(f'Found {len(tic_gaia)} unique TIC ids. Now querying TOI in FOV of LCs in the LC catalog\n')

    with open(tic_gaia_fname, 'w') as file:
        for tic_id, gaia_id in tic_gaia.items():
            file.write(f"{tic_id}, {gaia_id}\n")
    
    return tic_gaia


# ------------------------------------------------------------------#
# --- Instead of query vizier n times, query it once based on the --#
# --- FOV of the LCs and then find the TICs from the results -------#
# --- ORIGINAL TAP VERSION (kept for reference) ------------------- #
# ------------------------------------------------------------------#

def query_vizier_once_tap(gaia_ids, tic_gaia_fname, fov_mag_range):

    tic_gaia = {}

    query_tic_in_fov = tap_vizier_query(
        url='http://tapvizier.u-strasbg.fr/TAPVizieR/tap/',
        headers=['TIC', 'GAIA', 'Vmag', 'RAJ2000', 'DEJ2000'],
        table_database='"IV/39/tic82"',
    # TESS Input Catalog version 8.2 (TIC v8.2) (Paegert+, 2021) 1.7 Bilion rows
        constraints=[
            f"RAJ2000 BETWEEN {fov_mag_range['ra_min']} AND {fov_mag_range['ra_max']}",
            f"DEJ2000 BETWEEN {fov_mag_range['dec_min']} AND {fov_mag_range['dec_max']}",
            f"Vmag < {fov_mag_range['mag_max']} + 1.0"
        ]
    )

    print(f'Found {len(query_tic_in_fov)} TICs in the specified region.'
        f' Now check this result for each LC to find matches.\n')

    for row in query_tic_in_fov:
        if row['GAIA'] in gaia_ids:
            tic_gaia[row['TIC']] = row['GAIA']


    print(f'Found {len(tic_gaia)} matched TICs.\n')


    with open(tic_gaia_fname, 'w') as file:
        for tic_id, gaia_id in tic_gaia.items():
            file.write(f"{tic_id}, {gaia_id}\n")

    return tic_gaia


# ------------------------------------------------------------------#
# --- Instead of query vizier n times, query it once based on the - #
# --- FOV of the LCs and then find the TICs from the results -------#
# --- ASTROQUERY VERSION ------------------------------------------ #
# ------------------------------------------------------------------#

def query_vizier_once(gaia_ids, tic_gaia_fname, fov_mag_range):

    tic_gaia = {}

    result = query_tic_in_region(fov_mag_range)
    
    if result is None:
        print('No TICs found in the specified region.')
        return tic_gaia

    print(f'Found {len(result)} TICs in the specified region.'
        f' Now check this result for each LC to find matches.\n')

    for row in result:
        # print(f'row: {row}')  # Debug print to check the structure of the row
        gaia_val = row['GAIA']
        # print(f'gaia_val: {gaia_val}')  # Debug print to check the GAIA value
        # Handle cases where GAIA might be None or masked
        if gaia_val is not None and gaia_val != '':
            try:
                gaia_id = int(gaia_val) if isinstance(gaia_val, str) else gaia_val
                # print(f'gaia_id: {gaia_id}')  # Debug print to check the parsed Gaia ID
                if gaia_id in gaia_ids:
                    tic_gaia[int(row['TIC'])] = gaia_id
            except (ValueError, TypeError):
                continue

    print(f'Found {len(tic_gaia)} matched TICs.\n')

    with open(tic_gaia_fname, 'w') as file:
        for tic_id, gaia_id in tic_gaia.items():
            file.write(f"{tic_id}, {gaia_id}\n")

    return tic_gaia


# ------------------------------------------------------------------#
# --- Instead of query n times, do it once based on the ------------#
# --- FOV of the LCs to find the TOI IDs in that FOV,   ------------#
# --- then check to see if any of them exist in our LCs ------------#
# --- ORIGINAL TAP VERSION (kept for reference) ------------------- #
# ------------------------------------------------------------------#

def query_toi_in_fov_tap(tic_gaia, toi_gaia_period_fname, fov_mag_range):

    query_toi_in_fov = tap_vizier_query(
        url='https://exoplanetarchive.ipac.caltech.edu/TAP',
        headers='*',
        table_database='toi',
        constraints=[
            f"RA BETWEEN {fov_mag_range['ra_min']} AND {fov_mag_range['ra_max']}",
            f"DEC BETWEEN {fov_mag_range['dec_min']} AND {fov_mag_range['dec_max']}",
            f"ST_TMAG < {fov_mag_range['mag_max']} + 1.0"
        ]
    )

    print(f'Found {len(query_toi_in_fov)} TOIs in the specified region.'
          f' Now check LCs for each TOI to find matches.\n')

    tic_gaia_period = {}

    for row in query_toi_in_fov:
        if row['tid'] in tic_gaia:
            tic_gaia_period[(row['tid'], tic_gaia[row['tid']])] = row['pl_orbper']

    print(
        f'Found {len(tic_gaia_period)} TOI with their periods. Now writing to file '
        f'as: TIC, Gaia id, Period'
    )

    with open(toi_gaia_period_fname, 'w') as file:
        for (tic_id, gaia_id), period in tic_gaia_period.items():
            file.write(f"{tic_id}, {gaia_id}, {period}\n")


# ------------------------------------------------------------------#
# --- Instead of query n times, do it once based on the ------------#
# --- FOV of the LCs to find the TOI IDs in that FOV,   ------------#
# --- then check to see if any of them exist in our LCs ------------#
# --- ASTROQUERY VERSION ------------------------------------------ #
# ------------------------------------------------------------------#

def query_toi_in_fov(tic_gaia, toi_gaia_period_fname, fov_mag_range):

    # Query TOI catalog from exoplanet archive
    viz = Vizier(columns=['TOI', 'TIC', 'pl_orbper'])
    viz.ROW_LIMIT = -1
    
    # Calculate region center and radius
    ra_center = (fov_mag_range['ra_min'] + fov_mag_range['ra_max']) / 2
    dec_center = (fov_mag_range['dec_min'] + fov_mag_range['dec_max']) / 2
    ra_width = abs(fov_mag_range['ra_max'] - fov_mag_range['ra_min'])
    dec_width = abs(fov_mag_range['dec_max'] - fov_mag_range['dec_min'])
    radius = max(ra_width, dec_width) / 2 + 0.5
    
    coord = SkyCoord(ra=ra_center, dec=dec_center, unit=(u.deg, u.deg))
    
    # Query TOI catalog (available through Vizier)
    try:
        result = viz.query_region(coord, radius=radius * u.deg, catalog='IV/38/toi')
        if len(result) > 0:
            query_toi_result = result[0]
        else:
            query_toi_result = None
    except Exception as e:
        print(f"Error querying TOI catalog: {e}")
        query_toi_result = None

    if query_toi_result is None:
        print('No TOIs found in the specified region.')
        return

    print(f'Found {len(query_toi_result)} TOIs in the specified region.'
          f' Now check LCs for each TOI to find matches.\n')

    tic_gaia_period = {}

    for row in query_toi_result:
        try:
            tic_id = int(row['TIC'])
            if tic_id in tic_gaia:
                period = row['pl_orbper']
                if period is not None:
                    tic_gaia_period[(tic_id, tic_gaia[tic_id])] = float(period)
        except (ValueError, TypeError, KeyError):
            continue

    print(
        f'Found {len(tic_gaia_period)} TOI with their periods. Now writing to file '
        f'as: TIC, Gaia id, Period'
    )

    with open(toi_gaia_period_fname, 'w') as file:
        for (tic_id, gaia_id), period in tic_gaia_period.items():
            file.write(f"{tic_id}, {gaia_id}, {period}\n")



def main():

    cmdline_args = parse_command_line()
    if not (cmdline_args.lc_path and cmdline_args.lc_catalogs):
        print("Error: --lc-path and --lc-catalog are required")
        exit(1)

    gaia_ids = lcs_to_gaia_ids(cmdline_args.lc_path)
    df_lcs = pd.DataFrame({'gaia_id': list(gaia_ids)})
    # print(f"df_lcs.head():\n{df_lcs.head()}")
    if cmdline_args.skip_gaia_variables:
        print('Skipping Gaia variable query.')
    else:
        _ = gaia_ids_to_variables(gaia_ids)
    # print(f'ebb: {len(ebb)}, rr: {len(rr)}')

    fov_mag_range = find_fov_of_lcs(cmdline_args.lc_catalogs)

    if cmdline_args.tic_gaia_file:
        tic_gaia = {}
        with open(cmdline_args.tic_gaia_file, 'r') as file:
            for line in file:
                tic_id, gaia_id = map(int, line.strip().split(', '))
                tic_gaia[tic_id] = gaia_id

    elif len(gaia_ids) < 1000:  ### should find the efficient threshold
        tic_gaia = query_vizier_n_times(gaia_ids, cmdline_args.tic_gaia_fname)
    else:
        tic_gaia = query_vizier_once(gaia_ids, cmdline_args.tic_gaia_fname, fov_mag_range)

    query_toi_in_fov(tic_gaia, cmdline_args.toi_gaia_period_fname, fov_mag_range)



if __name__ == '__main__':
    main()
