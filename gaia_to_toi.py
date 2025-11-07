#!/usr/bin/env python3

"""
    Utilities for plotting individual lightcurves.
    It goes over the LCs, extract the Gaia Ids, then run tap_vizier_query
    for querying TIC for those ids.
    Then, it queries TOI table in exoplanet archive and extract the objects
    in the FOV of the LCs and check if any of those TOI have a TIC id in our
    list. If so, it extracts the period.
    Also, it writes two files:
        1- tic_gaia.txt: tic id, gaia id (for later use cases)
        2- toi_gaia_period.txt: toi id, gaia id, period

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


import pyvo
import os
import os.path
import glob
import ast
import numpy
import tenacity
import pandas as pd
from astropy.io import fits

from configargparse import ArgumentParser, DefaultsFormatter

def parse_command_line():
    """Return the parsed command line arguments."""

    parser = ArgumentParser(
        description=__doc__
    )
    parser.add_argument(
        '--lc-path',
        default=None,
        help='The path containing .h5 lc files.'
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
        '--lc-catalog',
        default='lc_catalog.fits',
        help='The input file containing the light curve catalog.'
    )
    # parser.add_argument(
    #     '--tables-to-query',
    #     default=['"IV/39/tic82"', 'toi'],
    #     help='The tables to query respectively.'
    # )


    return parser.parse_args()

# general-purpose function for any TAP service
# @tenacity.retry(wait=tenacity.wait_fixed(2), stop=tenacity.stop_after_attempt(5), reraise=True)
# def tap_query(url,
#               headers,
#               table_database,
#               constraints=None):
#     """
#         The TAP query to run. All queries must have valid ADQL, which takes the form:
#         SELECT <column list> FROM <table> WHERE <constraints>

#         Args:
#             - url: The url to use for TAP query
#             - headers (list) : The table column headers to use from given table database
#             - table_database: The table database name to use
#             - constraints (str list): The given constraints to use for TAP query

#         Returns:
#             - results: The TAP query results as a VOTable
#     """

#     service = pyvo.dal.TAPService(url)

#     headers = ", ".join(str(item) for item in headers)
#     constraints = " AND ".join(str(item) for item in constraints)
#     query = "SELECT " + str(headers) + " FROM " + str(table_database) + " WHERE " + str(constraints)
    
#     results = service.search(query)

#     return results

# specialized for Vizier TAP service
@tenacity.retry(wait=tenacity.wait_fixed(2), stop=tenacity.stop_after_attempt(5), reraise=True)
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

    print(f'Found {len(gaia_ids)} unique gaia ids. Now querying TIC\n')

    return gaia_ids


# ------------------------------------------------------------------#
# ----------- Find the FOV of the LCs ------------------------------#
# ------------------------------------------------------------------#

def find_fov_of_lcs(lc_catalog):

    with fits.open(lc_catalog) as hdul:
        lc_data = hdul[1].data

        ra_min, ra_max = lc_data['ra'].min(), lc_data['ra'].max()
        dec_min, dec_max = lc_data['dec'].min(), lc_data['dec'].max()
        mag_min, mag_max = lc_data['phot_g_mean_mag'].min(), lc_data['phot_g_mean_mag'].max()

        print(f'Ranges in lc catalog:::\n'
              f'ra: ({ra_min:.2f}, {ra_max:.2f})\n'
              f'dec: ({dec_min:.2f}, {dec_max:.2f})\n'
              f'mag: ({mag_min:.2f}, {mag_max:.2f})\n'
        )

    return {'ra_min': ra_min,
            'ra_max': ra_max,
            'dec_min': dec_min,
            'dec_max': dec_max,
            'mag_min': mag_min,
            'mag_max': mag_max
            }


# ------------------------------------------------------------------#
# --- For each GAIA ID, query vizier to find its TIC (not efficient sometimes) ---#
# ------------------------------------------------------------------#

def query_vizier_n_times(gaia_ids, tic_gaia_fname):

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
# --- Instead of query vizier n times, query it once based on the --#
# --- FOV of the LCs and then find the TICs from the results -------#
# ------------------------------------------------------------------#

def query_vizier_once(gaia_ids, tic_gaia_fname, fov_mag_range):

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
# --- Instead of query n times, do it once based on the ------------#
# --- FOV of the LCs to find the TOI IDs in that FOV,   ------------#
# --- then check to see if any of them exist in our LCs ------------#
# ------------------------------------------------------------------#

def query_toi_in_fov(tic_gaia, toi_gaia_period_fname, fov_mag_range):

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
    # df_toi = query_toi_in_fov.to_table().to_pandas()
    # print(df_toi.columns)

    # Now, go over the query_result and for each column tid, check if we have it in
    # the values of gaia_tic
    tic_gaia_period = {}

    # for tic in query_toi_in_fov['tid']: # Instead of tic, we have tid in the query, and also toi
    for row in query_toi_in_fov:
        if row['tid'] in tic_gaia:
            # print(f'Found TOI with tid {row["tid"]} and GAIA id: {tic_gaia[row["tid"]]} in our LC catalog')
            # gaia_id = query_result[query_result['tid'] == tid]['gaia_id'].values[0]
            # tic_id = tid
            tic_gaia_period[(row['tid'], tic_gaia[row['tid']])] = row['pl_orbper']

    print(
        f'Found {len(tic_gaia_period)} TOI with their periods. Now writing to file '
        f'as: TIC, Gaia id, Period'
    )

    with open(toi_gaia_period_fname, 'w') as file:
        for (tic_id, gaia_id), period in tic_gaia_period.items():
            file.write(f"{tic_id}, {gaia_id}, {period}\n")



# -----------------------------------------------------------------#
# ---- For each TIC ID that we found before, ----------------------#
# ---- query it to see if it has a TOI (not efficient) ------------#
# -----------------------------------------------------------------#

# def ... :

    # gaia_tic_period = {}

    # n = 0
    # for gaia_id, tic_id in gaia_tic.items():
    #     n += 1
    #     if n % 10 == 0: print(f'TOI Query result: passed {n}')

    #     query_result = tap_query(
    #         url='https://exoplanetarchive.ipac.caltech.edu/TAP',
    #         headers='*',
    #         table_database='toi',
    #         constraints=['tid=' + str(tic_id)]
    #     )

    #     if len(query_result) > 0:
    #         gaia_tic_period[(gaia_id, tic_id)] = query_result['pl_orbper'].data[0]

    # print(
    #     f'Found {len(gaia_tic_period)} TOI with their periods. now writing to file '
    #     f'as: Gaia id, TOI id, Period'
    #       )

    # with open(cmdline_args.gaia_toi_file, 'w') as file:
    #     for (gaia_id, tic_id), period in gaia_tic_period.items():
    #         file.write(f"{gaia_id}, {tic_id}, {period}\n")

def main():

    cmdline_args = parse_command_line()
    if not (cmdline_args.lc_path and cmdline_args.lc_catalog):
        print("Error: --lc-path and --lc-catalog are required")
        exit(1)

    gaia_ids = lcs_to_gaia_ids(cmdline_args.lc_path)
    fov_mag_range = find_fov_of_lcs(cmdline_args.lc_catalog)

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


    ####### What is the point of this? which is not working anyway!
    
    # query_test = tap_vizier_query(
    # url='https://mast.stsci.edu/api/v0/tap',
    # headers='*',
    # table_database='tic_v8',
    # constraints=[
    #     f"RA BETWEEN 10 AND 12",
    #     f"DEC BETWEEN 15 AND 18",
    #     f"ST_TMAG < 10"
    #     ]
    # )
    # df_test = query_test.to_table().to_pandas()
    # print(df_test.columns)



if __name__ == '__main__':
    main()
