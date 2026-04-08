#!/usr/bin/env python3

import lightkurve as lk
import matplotlib.pyplot as plt
import numpy as np


search = lk.search_lightcurve("Gaia DR3 3337025761759515776", mission="TESS")
# options: "TOI 700", "TIC 150428135", "Gaia DR3 4653803237223738240"


print(search)

lc = search.download()    # download the first available sector
# lc = search.download_all()    # multiple sectors

lc = lc.remove_nans().normalize()
lc = lc.flatten(window_length=401) # remove long-term trends

time_tess = lc.time.value      # already in BTJD
flux_tess = lc.flux.value

# Convert flux to magnitude (relative, since normalized flux)
mag_tess = -2.5 * np.log10(flux_tess)



# Phase-fold the light curve using the known period.
period = 1.7801492  # days
reference_epoch = 1469.56553  # choose a reference time; replace with a known epoch if available
phase = ((time_tess - reference_epoch) / period) % 1
phase = (phase + 0.25) % 1 # shift so that primary eclipse is at phase 0.25
sort_idx = np.argsort(phase)
phase = phase[sort_idx]
flux_fold = flux_tess[sort_idx]
mag_fold = mag_tess[sort_idx]

fig, ax = plt.subplots(2,1, figsize=(10,8))
ax[0].scatter(phase, flux_fold, s=1, label="TESS", alpha=0.7)
ax[0].set_xlabel("Phase")
ax[0].set_ylabel("Normalized Flux")
ax[0].set_title(f"Phase-folded TESS LC (P = {period} d)")
ax[0].legend()

ax[1].scatter(phase, mag_fold, s=1, label="TESS", alpha=0.7)
ax[1].set_xlabel("Phase")
ax[1].set_ylabel("Relative Magnitude")
ax[1].set_title("TESS LC: Gaia DR3 3337025761759515776, EB Villa, P= 1.7, Tmag=9.57")
ax[1].invert_yaxis()
ax[1].legend()

plt.tight_layout()
plt.show()