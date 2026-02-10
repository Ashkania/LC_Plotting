# plotting

A: Running gaia_to_toi.py

Utilities for plotting individual lightcurves.
    
   Command line arguments:
    --lc-path: The path containing .h5 lc files.
    --tic-gaia-fname: The output file to write tic to gaia mappings.
        --tic-gaia-file: if we already have it, specify the file path
    --toi-gaia-period-fname: The output file to write toi, gaia, period mappings.
    --lc-catalog: The input file containing the light curve catalog.

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


B: Runnig plotting.py
    Plot the LC from its hdf5 file
