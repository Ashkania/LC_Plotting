# plotting

Utilities for plotting individual lightcurves.
   It goes over the LCs, extract the Gaia Ids, then run tap_vizier_query
   for querying TIC for those ids.  
   
   Then, it queries TOI table in exoplanet archive and extract the objects
   in the FOV of the LCs and check if any of those TOI have a TIC id in our
   list. If so, it extracts the period.
   
   Also, it writes two files:<br>
       1- tic_gaia.txt: tic id, gaia id (for later use cases)<br>
       2- toi_gaia_period.txt: toi id, gaia id, period<br>

   Check this link for all tables that you may need:
   https://tapvizier.cds.unistra.fr/adql/
