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
    parser.add_argument(
        "--title",
        type=str,
        default="Light Curve",
        help="Title for the LC plot"
    )
    parser.add_argument(
        "--period",
        type=float,
        default=None,
        help="Period in days for phase folding (optional)"
    )
    parser.add_argument(
        "--epoch",
        type=float,
        default=None,
        help="Reference epoch (t0) for phase folding (defaults to first BJD)"
    )
    
    return parser.parse_args()

def read_lc_from_hdf5(hdf5_file, aperture=0):
    with h5py.File(hdf5_file, 'r') as f:
        bjd = f[f'SkyPosition/BJD'][:]
        mag = f[f'AperturePhotometry/Aperture{aperture:03}/MagnitudeFitting/Magnitude'][:]
        mag[mag < -1e5] = np.nan  # Replace Nan (large negative) values with NaN for better plotting
    return bjd, mag

def calculate_phase(bjd, period, epoch=None):
    """
    Calculate phase for light curve folding.
    
    Parameters:
    -----------
    bjd : array
        Barycentric Julian Date array
    period : float
        Period in days
    epoch : float, optional
        Reference epoch (t0). If None, uses the first BJD value.
    
    Returns:
    --------
    phase : array
        Phase values between 0 and 1
    """
    if epoch is None:
        epoch = bjd[0]
    
    phase = ((bjd - epoch) / period) % 1.0
    return phase

def plot_light_curve(bjd, mag, title="Light Curve", xlabel="BJD", ylabel="Magnitude", period=None):
    """
    Plot light curve. If period is provided, plots two cycles for continuity.
    """
    plt.figure(figsize=(10, 5))
    
    plt.scatter(bjd, mag, s=10, color='green')
    
    plt.gca().invert_yaxis()  # Magnitude axis is inverted
    plt.title(title)
    plt.xlabel(xlabel, fontdict={"fontweight":"bold", 'fontsize':14})
    plt.ylabel(ylabel, fontdict={"fontweight":"bold", 'fontsize':14})
    plt.grid(True)
    plt.savefig(f'{title}',dpi=400)
    plt.show()

if __name__ == "__main__":

    args = parse_arguments()
    hdf5_file = args.hdf5_file
    aperture = args.aperture
    title = args.title
    period = args.period
    epoch = args.epoch

    bjd, mag = read_lc_from_hdf5(hdf5_file, aperture)
    
    if period is not None:
        # Phase fold the light curve
        phase = calculate_phase(bjd, period, epoch)
        plot_light_curve(phase, mag, title=title, xlabel="Phase", period=period)
        print(f"Light curve folded with period = {period} days")
        if epoch:
            print(f"Reference epoch (t0) = {epoch}")
        else:
            print(f"Reference epoch (t0) = {bjd[0]} (first observation)")
    else:
        # Plot regular light curve
        plot_light_curve(bjd, mag, title=title)