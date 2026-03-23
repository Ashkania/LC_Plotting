#!/usr/bin/env python3

# TODO: for all the photrefs, or some of them, like the
# best ones which should be 12, 20, 28
# (red channels for Maia60, Maia04, Asterope08) go over
# each apperture, extract the magnitudes, find the mad of them,
# plot the minimum value vs the magnitude of the star. It must
# be the same as the mag/mag plot from magfit step.

# For a given light curve, for each photref(= each camera and channel),
# we find the best aperture (the one with the lowest mad) and plot
# the mad vs ... of the star. This should match the mag/mag plot from the magfit step.

import numpy as np
import matplotlib.pyplot as plt
import argparse
import h5py

def parse_args():
    
    parser = argparse.ArgumentParser(
        description="Diagnostic plots for light curve analysis"
        )
    parser.add_argument(
        'lc_files',
        nargs='+',
        help="Input hdf5 LC files"
        )
    
    return parser.parse_args()

def load_data(file, photref=0):
    with h5py.File(file, 'r') as f:
        photrefs = f[f'AperturePhotometry/Aperture004/MagnitudeFitting/ConfigurationIndex'][:]
        
        data = {}
        mask = photrefs == photref
        for aperture in range(8):
            epd_all = f[f'AperturePhotometry/Aperture{aperture:03}/EPD/Magnitude'][:]
            epd = epd_all[mask]
            mad_epd = np.nanmedian(np.abs(epd - np.nanmedian(epd)))
            data[aperture] = mad_epd
    
    return data




def main():
    args = parse_args()
    
    for file in args.lc_files:
        print(f'Processing file: {file}')
        for photref in [12, 28, 20]:
        # for photref in range(40):  # check all photrefs up to 40
            data = load_data(file, photref=photref)
            MAD = np.min([data[aperture] for aperture in range(8)])
            ap = np.argmin([data[aperture] for aperture in range(8)])
            plt.plot(MAD, 'or', label=f'Photref {photref}')
            print(f'Photref: {photref}, Best Aperture: {ap},MAD: {MAD}')
        plt.show()

if __name__ == "__main__":
    main()
    