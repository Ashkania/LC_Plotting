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
        '--gaia-toi-file',
        default='output.txt',
        help='The output file to write gaia to tic mappings.'
    )

    return parser.parse_args()

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

    print(f'Found {len(gaia_ids)} unique gaia ids, now querying TIC')

    gaia_tic = {}

    n = 0
    for gaia_id in gaia_ids:
        n += 1
        if n % 100 == 0:
            print(f'TIC Query result: passed {n}')

        query_result = tap_vizier_query(
            url='http://tapvizier.u-strasbg.fr/TAPVizieR/tap/',
            headers='*', # or just TIC column?
            table_database='"IV/39/tic82"', 
            # TESS Input Catalog version 8.2 (TIC v8.2) (Paegert+, 2021) 1.7 Bilion rows
            constraints=['GAIA=' + str(gaia_id)]
        )
        # it returns an empty result if nothing returned, not an error
        # so no need to try n except

        if len(query_result) > 0:
            gaia_tic[gaia_id] = query_result['TIC'].data[0]


    print(f'Found {len(gaia_tic)} unique tic ids. now querying TIC for periods')


    gaia_tic_period = {}

    n = 0
    for gaia_id, tic_id in gaia_tic.items():
        n += 1
        if n % 100 == 0:
            print(f'TOI Query result: passed {n}')

        query_result = tap_query(
            url='https://exoplanetarchive.ipac.caltech.edu/TAP',
            headers='*',
            table_database='toi',
            constraints=['tid=' + str(tic_id)]
        )

        if len(query_result) > 0:
            gaia_tic_period[(gaia_id, tic_id)] = query_result['pl_orbper'].data[0]

    print(
        f'Found {len(gaia_tic_period)} TOI with their periods. now writing to file '
        f'as: Gaia id, TOI id, Period'
          )

    with open(cmdline_args.gaia_toi_file, 'w') as file:
        for (gaia_id, tic_id), period in gaia_tic_period.items():
            file.write(f"{gaia_id}, {tic_id}, {period}\n")
