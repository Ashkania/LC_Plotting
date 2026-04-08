#!/usr/bin/env python3

import lightkurve as lk
import matplotlib.pyplot as plt


search = lk.search_lightcurve("Gaia DR3 3337025761759515776", mission="TESS")
# options: "TOI 700", "TIC 150428135", "Gaia DR3 4653803237223738240"


print(search)

lc = search.download()    # download the first available sector
# lc = search.download_all()    # multiple sectors

lc = lc.remove_nans().normalize()
lc = lc.flatten(window_length=401) # remove long-term trends

time_tess = lc.time.value      # already in BTJD
flux_tess = lc.flux.value


plt.figure(figsize=(10,5))

plt.scatter(time_tess, flux_tess, s=1, label="TESS", alpha=0.7)
# plt.scatter(time_my, flux_my, s=5, label="My Data", alpha=0.7)

plt.xlabel("Time (BTJD)")
plt.ylabel("Normalized Flux")
plt.legend()
plt.show()