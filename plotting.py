#!/usr/bin/env python3

"""
    Plot the LC from its hdf5 file
    To do: how to select best aperture for different kinds of variability
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

def parse_arguments():
    
    parser = ArgumentParser(
        description="Plot Light Curve from HDF5 file"
        )
    parser.add_argument(
        "hdf5_file",
        type=str,
        help="Path to the HDF5 file containing the light curve data"
        )
    parser.add_argument(
        "--aperture",
        type=int,
        default=0,
        help="Aperture number to plot (default: 0)"
        )
    
    return parser.parse_args()

def read_lc_from_hdf5(hdf5_file, aperture=0):
    with h5py.File(hdf5_file, 'r') as f:
        bjd = f[f'SkyPosition/BJD'][:]
        mag = f[f'AperturePhotometry/Aperture{aperture:03}/MagnitudeFitting/Magnitude'][:]
    return bjd, mag

def plot_light_curve(bjd, mag, title="Light Curve", xlabel="BJD", ylabel="Magnitude"):
    plt.figure(figsize=(10, 5))
    plt.scatter(bjd, mag, s=1, color='blue')
    plt.gca().invert_yaxis()  # Magnitude axis is inverted
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.show()

if __name__ == "__main__":

    args = parse_arguments()
    hdf5_file = args.hdf5_file
    aperture = args.aperture

    bjd, mag = read_lc_from_hdf5(hdf5_file, aperture)
    plot_light_curve(bjd, mag, title=f"Light Curve for Aperture {aperture}")