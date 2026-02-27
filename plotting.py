#!/usr/bin/env python3

"""
    Plot the LC from its hdf5 file
    IMPORTANT: Edit the plot based on the construction of LCs:
        4 channels from all cameras comes together in LC.
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

def read_lc_from_hdf5(hdf5_file, aperture=4):

    times = {}
    magnitudes = {}

    with h5py.File(hdf5_file, 'r') as f:
        channels = f[f'FrameInformation/ConfigurationIndex'][:]   # [0 0 0 1 1 2 2 3 3]
        bjds = f[f'SkyPosition/BJD'][:]
        mag_epd = f[f'AperturePhotometry/Aperture{aperture:03}/EPD/Magnitude'][:]
        mag_epd[mag_epd < -1e5] = np.nan  # Replace Nan (large negative) values with NaN
        mag_magfit = f[f'AperturePhotometry/Aperture{aperture:03}/MagnitudeFitting/Magnitude'][:]
        mag_magfit[mag_magfit < -1e5] = np.nan  # Same for magfit

        for channel in range(4):
            times[channel] = bjds[channels == channel]
            for mag_type, mag_array in zip(['epd', 'magfit'], [mag_epd, mag_magfit]):
                magnitudes[f'{mag_type}_{channel}'] = mag_array[channels == channel]
            
    # times = {0: array([...]), 1: array([...]), ...}
    # magnitudes = {'epd_0': array([...]), 'magfit_0': array([...]), 'epd_1': array([...]), ...}
    return times, magnitudes

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

def plot_light_curve(bjd, mag_epd, mag_magfit, channel, title, xlabel="BJD", ylabel="Magnitude", period=None):
    """
    Plot light curve. If period is provided, plots two cycles for continuity.
    """
    plt.figure(figsize=(10, 5))
    
    plt.scatter(bjd, mag_epd, s=10, color='red', alpha=0.5, label='EPD')
    plt.scatter(bjd, mag_magfit, s=10, color='blue', alpha=0.5, label='MagFit')
    
    plt.gca().invert_yaxis()  # Magnitude axis is inverted
    plt.title(f'{title} - Channel: {channel}', fontdict={"fontweight":"bold", 'fontsize':16})
    plt.xlabel(xlabel, fontdict={"fontweight":"bold", 'fontsize':14})
    plt.ylabel(ylabel, fontdict={"fontweight":"bold", 'fontsize':14})
    plt.grid(True)
    plt.legend()
    plt.savefig(f'{title}_Channel_{channel}',dpi=400)
    plt.show()


def main():

    args = parse_arguments()
    hdf5_file = args.hdf5_file
    title = hdf5_file.split('/')[-1].split('.')[0]
    aperture = args.aperture
    period = args.period
    epoch = args.epoch

    times, magnitudes = read_lc_from_hdf5(hdf5_file, aperture)
    # times, magnitudes = {0: array([...]), 1: array([...]), 2: array([...]), 3: array([...])}

    # TODO: Plot all channels together with different colors and legends    
    for channel in sorted(times.keys()):
        bjd = times[channel]
        mag_epd = magnitudes[f'epd_{channel}']
        mag_magfit = magnitudes[f'magfit_{channel}']
    
        # If period is provided, phase fold the first channel's data
        if period is not None:
            phase = calculate_phase(bjd, period, epoch)
            plot_light_curve(phase, mag_epd, mag_magfit, channel=channel, title=title, xlabel="Phase", period=period)
            print(f"Light curve folded with period = {period} days")
            if epoch:
                print(f"Reference epoch (t0) = {epoch}")
            else:
                print(f"Reference epoch (t0) = {bjd[0]} (first observation)")
        else:
            # Plot regular light curve
            plot_light_curve(bjd, mag_epd, channel=channel, title=title)


if __name__ == "__main__":
    main()