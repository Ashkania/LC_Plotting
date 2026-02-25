python3 -u ./gaia_to_toi.py \
  --lc-path ./LC \
  --lc-catalogs ./catalogs/* \
  --skip-gaia-variables \
  --tic-gaia-fname tic_gaia_astroquery.txt \
  --toi-gaia-period-fname toi_gaia_astroquery.txt \

#  > gaia_run.log 2>&1 &
