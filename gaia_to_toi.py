#!/usr/bin/env python3

"""
   Utilities for plotting individual lightcurves.
   It goes over the LCs, extract the Gaia Ids, then run tap_vizier_query
   for querying TIC for those ids. Finally, give out the TOI as their Gaia ids
   with their periods in a txt file like:
   [3345187475230297856, 3.817421]

   Check this link for all tables that you may need:
   https://tapvizier.cds.unistra.fr/adql/
   more info in discord:
   https://discord.com/channels/581176949184135229/1291434287614398494/1356713434086903970
   
   !!! This writes to the file all those 21000 tic that it finds which is not
   what we want. Compare with Angel's to see how to get rid of the first write to
   the file
   !!!

   !!! Maybe, there is(?) a problem that make it stops after
   passing 1000, frames while the Angel's script don't. This
   script tried to fix the issue of multiply querying of Angel's
   script, but seems to have a problem
   !!!
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
        '--gaia-tic-file',
        default='gaia_tic.txt',
        help='The output file to write gaia to tic mappings.'
    )
    parser.add_argument(
        '--gaia-toi-file',
        default='gaia_toi.txt',
        help='The output file to write gaia to toi mappings.'
    )


    return parser.parse_args()

# general-purpose function for any TAP service
@tenacity.retry(wait=tenacity.wait_fixed(2), stop=tenacity.stop_after_attempt(5), reraise=True)
def tap_query(url,
              headers,
              table_database,
              constraints=None):
    """
        The TAP query to run. All queries must have valid ADQL, which takes the form:
        SELECT <column list> FROM <table> WHERE <constraints>

        Args:
            - url: The url to use for TAP query
            - headers (list) : The table column headers to use from given table database
            - table_database: The table database name to use
            - constraints (str list): The given constraints to use for TAP query

        Returns:
            - results: The TAP query results as a VOTable
    """

    service = pyvo.dal.TAPService(url)

    constraints = " AND ".join(str(item) for item in constraints)
    headers = ", ".join(str(item) for item in headers)

    query = "SELECT " + str(headers) + " FROM " + str(table_database) + " WHERE " + str(constraints)

    #TODO add argparse argument for async or sync query
    #turn this on for synchronous query
    results = service.search(query)

    #turn this on for asynchronous query
    # results = service.run_async(query)

    return results

# specialized for Vizier TAP service
@tenacity.retry(wait=tenacity.wait_fixed(2), stop=tenacity.stop_after_attempt(5), reraise=True)
def tap_vizier_query(url='http://tapvizier.u-strasbg.fr/TAPVizieR/tap/',
                     headers='*',
                     table_database='"J/ApJS/258/16/tess-ebs"', 
                     # TESS Eclipsing Binary stars. I. Sectors 1-26 (Prsa+, 2022) # 4.5 thousands rows
                     constraints=None):
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

    if constraints is None:
        query = "SELECT " + str(headers) + " FROM " + str(table_database)

    elif len(constraints)==1:
        query = (
            "SELECT " + str(headers)
            + " FROM " + str(table_database)
            + " WHERE " + str(table_database)
            + "." + str(constraints[0])
            )

    #TODO this needs to get fixed for querying a list of ids
    # if id_list:
        # query = "SELECT headers FROM table_database WHERE table_database....#TODO add constraints#
        #           AND table_database.ID IN (%s)" % (', '.join(str(id) for id in id_list))

    else:
        query = ("SELECT " + str(headers) + " FROM " + str(table_database) + " WHERE " + str(table_database) + "."
                 + str(constraints[0])+ "".join([(" AND " + str(table_database) + "." + str(constraint))
                                                 for constraint in constraints[1:]]))

    # TODO add argparse argument for async or sync query
    # turn this on for synchronous query
    results = service.search(query)

    # turn this on for asynchronous query
    # results = service.run_async(query)

    return results


if __name__ == '__main__':

    cmdline_args = parse_command_line()
    if not cmdline_args.lc_path:
        print("Error: --lc-path is required")
        exit(1)

    # Gather all lcs full paths:
    lcs = sorted(glob.glob(os.path.join(cmdline_args.lc_path, '*.h5')))    
    # lcs = ['/path/to/GDR3_1316714794020532096.h5', ...]

    # from lc full paths to the int(Gaia id):
    gaia_ids = {
        int(os.path.basename(lc).split('_')[1].split('.')[0])
        for lc in lcs
    }
    # gaia_ids = {1316714794020532096, ...}

    print(f'Found {len(gaia_ids)} unique gaia ids. Now querying TIC\n')

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
        # it returns an empty result if nothing returned, not an error
        # so no need to try n except

        if len(query_gaia_to_tic) > 0:
            tic_id = query_gaia_to_tic['TIC'].data[0]
            tic_gaia[tic_id] = gaia_id


    # print(f'Found {len(gaia_tic)} unique tic ids. Now querying TIC for periods')
    print(f'Found {len(tic_gaia)} unique TIC ids. Now querying TOI in FOV of LCs in the LC catalog\n')


    with open(cmdline_args.gaia_tic_file, 'w') as file:
        for gaia_id, tic_id in tic_gaia.items():
            file.write(f"{gaia_id}, {tic_id}\n")

    gaia_tic_period = {}

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

#################################################################33
    # Try to first query all toi in a spacial constraint, then check if 
    # we have any LC in that list

    with fits.open('lc_catalog.fits') as hdul:
        lc_data = hdul[1].data

        ra_min, ra_max = lc_data['ra'].min(), lc_data['ra'].max()
        dec_min, dec_max = lc_data['dec'].min(), lc_data['dec'].max()
        mag_min, mag_max = lc_data['phot_g_mean_mag'].min(), lc_data['phot_g_mean_mag'].max()

        print(f'Ranges in lc catalog:::\n'
              f'ra: ({ra_min:.2f}, {ra_max:.2f})\n'
              f'dec: ({dec_min:.2f}, {dec_max:.2f})\n'
              f'mag: ({mag_min:.2f}, {mag_max:.2f})\n'
        )

    query_toi_in_fov = tap_query(
        url='https://exoplanetarchive.ipac.caltech.edu/TAP',
        headers='*',
        table_database='toi',
        constraints=[
            f"RA BETWEEN {ra_min} AND {ra_max}",
            f"DEC BETWEEN {dec_min} AND {dec_max}",
            f"ST_TMAG < {mag_max} + 1.0"
        ]
    )

    print(f'Found {len(query_toi_in_fov)} TOIs in the specified region.'
          f' Now check LCs for each TOI to find matches.\n')
    # df_toi = query_toi_in_fov.to_table().to_pandas()
    # print(df_toi.columns)

    # Now, go over the query_result and for each column tid, check if we have it in
    # the values of gaia_tic

    tic_gaia_period = {}

    for tic in query_toi_in_fov['tid']: # Instead of tic, we have tid in the query, and also toi
        if tic in tic_gaia:
            print(f'Found TOI with tid {tic} and GAIA id: {tic_gaia[tic]} in our LC catalog')
            # gaia_id = query_result[query_result['tid'] == tid]['gaia_id'].values[0]
            # tic_id = tid
            # tic_gaia_period[(tic, tic_gaia[tic])] = query['pl_orbper'].data[0]

    # print(
    #     f'Found {len(gaia_tic_period)} TOI with their periods. now writing to file '
    #     f'as: Gaia id, TOI id, Period'
    #       )

    with open(cmdline_args.gaia_toi_file, 'w') as file:
        for (gaia_id, tic_id), period in gaia_tic_period.items():
            file.write(f"{gaia_id}, {tic_id}, {period}\n")
