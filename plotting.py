#!/usr/bin/env python3

"""
    Plot the LC from its hdf5 file
    To do: how to select best aperture for different kinds of variability
    example:
    python plotting.py LC/EB_Villa/GDR3_3337025761759515776.h5 
        --aperture 4 --period 1.7801492 --epoch 1469.56553 --fold-by 'photref'
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
        default=4,
        help="Aperture number to plot (default: 4)"
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
        help="Reference epoch (t0) for phase folding"
    )
    parser.add_argument(
        "--fold-by",
        type=str,
        choices=['channel', 'photref'],
        default='photref',
        help="Whether to fold by 'channel' or 'photref' (default: 'photref')"
    )
    parser.add_argument(
        "--plot-type",
        choices=['individual', 'combined'],
        default='combined',
        help="Plot all channels/photrefs in one figure"
    )

    
    return parser.parse_args()

def read_lc_from_hdf5(hdf5_file, aperture=4, fold_by='photref'):
    """
    Read light curve data from HDF5 file and organize it for plotting.
     Parameters:
    -----------
    hdf5_file : str
        Path to the HDF5 file containing the light curve data
    aperture : int, optional
        Aperture number to read (default: 4)
    fold_by : str, optional
        Whether to fold by 'channel' or 'photref' (default: 'photref')
    Returns:
    --------
    times : dict
        Dictionary of time arrays keyed by channel or photref
        ex: times = {0: array([...]), 1: array([...]), ...}
    magnitudes : dict
        Dictionary of magnitude arrays keyed by 'epd_channel', 'magfit_channel', etc.
        ex: magnitudes = {'epd_0': array([...]), 'magfit_0': array([...]), ...}
    """

    with h5py.File(hdf5_file, 'r') as f:
        channels = f[f'FrameInformation/ConfigurationIndex'][:]
        photrefs = f[f'AperturePhotometry/Aperture004/MagnitudeFitting/ConfigurationIndex'][:]
        
        bjds = f[f'SkyPosition/BJD'][:]
        
        mag_epd = f[f'AperturePhotometry/Aperture{aperture:03}/EPD/Magnitude'][:]
        mag_magfit = f[f'AperturePhotometry/Aperture{aperture:03}/MagnitudeFitting/Magnitude'][:]

        for mag_array in [mag_epd, mag_magfit]:
            mag_array[mag_array < -1e5] = np.nan  # Replace NaN (large negative) values with NaN
        
        times = {}
        magnitudes = {}

        if fold_by == 'channel':
            for chn in np.unique(channels):
                times[chn] = bjds[channels == chn]
                for mag_type, mag_array in zip(['epd', 'magfit'], [mag_epd, mag_magfit]):
                    magnitudes[f'{mag_type}_{chn}'] = (
                        mag_array[channels == chn] - np.nanmedian(mag_array[channels == chn])
                        )  # Detrend by subtracting median

        elif fold_by == 'photref':
            for phref in np.unique(photrefs):
                times[phref] = bjds[photrefs == phref]
                for mag_type, mag_array in zip(['epd', 'magfit'], [mag_epd, mag_magfit]):
                    magnitudes[f'{mag_type}_{phref}'] = (
                        mag_array[photrefs == phref] - np.nanmedian(mag_array[photrefs == phref])
                        )

    return times, magnitudes

def calculate_phase(bjd, period, epoch):
    """
    Calculate phase for light curve folding.
    
    Parameters:
    -----------
    bjd : array
        Barycentric Julian Date array
    period : float
        Period in days
    epoch : float, optional
        Reference epoch (t0).
    
    Returns:
    --------
    phase : array
        Phase values between 0 and 1
    """
    
    phase = ((bjd - epoch) / period) % 1.0
    return phase

def plot_light_curve(bjd, mag_epd, mag_magfit, fold_by, object_name, xlabel, period=None):
    """
    Plot light curve.
    """
    # plt.figure(figsize=(10, 5))
    plt.scatter(bjd, mag_epd, s=10, color='red', alpha=0.5, label='EPD')
    plt.scatter(bjd, mag_magfit, s=10, color='blue', alpha=0.5, label='MagFit')
    
    plt.gca().invert_yaxis()  # Magnitude axis is inverted
    plt.xlabel(xlabel, fontdict={"fontweight":"bold", 'fontsize':14})
    plt.ylabel("Magnitude", fontdict={"fontweight":"bold", 'fontsize':14})
    plt.grid(True)
    plt.legend()
    plt.savefig(f'{object_name}_FoldBy_{fold_by}',dpi=400)
    # if plot one by one: 
        # plt.ylim(1.1, -0.6)
        # plt.title(f'{object_name} - Sector/Chn: {fold_by}', fontdict={"fontweight":"bold", 'fontsize':16})
        # plt.show()

def main():

    args = parse_arguments()
    hdf5_file = args.hdf5_file
    object_name = hdf5_file.split('/')[-1].split('.')[0]
    aperture = args.aperture
    period = args.period
    epoch = args.epoch
    fold_by = args.fold_by
    plot_type = args.plot_type

    times, magnitudes = read_lc_from_hdf5(hdf5_file, aperture, fold_by)

    for channel in sorted(times.keys()):
        x_values = times[channel]
        mag_epd = magnitudes[f'epd_{channel}']
        mag_magfit = magnitudes[f'magfit_{channel}']
    
        if period:
            phase = calculate_phase(x_values, period, epoch)
            plot_light_curve(
                phase, mag_epd, mag_magfit,
                fold_by, object_name,
                xlabel="Phase", period=period
                )
        else:
            # Plot regular light curve
            bjd = x_values - 2450000
            plot_light_curve(
                bjd, mag_epd, mag_magfit,
                fold_by, object_name,
                xlabel="BJD - 2450000"
                )
    if plot_type == 'combined':
        plt.title(
            f'{object_name} - folded by: {fold_by}',
            fontdict={"fontweight":"bold", 'fontsize':16}
            )
        plt.ylim(1.5, -0.6)
        plt.show()


if __name__ == "__main__":
    main()