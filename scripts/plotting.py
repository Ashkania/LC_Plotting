#!/usr/bin/env python3

"""
    Plot the LC from its hdf5 file
    To do: how to select best aperture for different kinds of variability
    example:
    python plotting.py LC/EB_Villa/GDR3_3337025761759515776.h5 
        --aperture 4 --period 1.7801492 --epoch 1469.56553 --fold-by 'photref'
        --plot-type 'combined' --binning 'time' 10
"""

import h5py
import numpy as np
import pandas as pd
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
    parser.add_argument(
        "--binning",
        nargs=2,
        default=None,
        metavar=('METHOD', 'SIZE'),
        help="Binning method (time/phase) and size (minutes for time, points for phase)"
    )

    
    return parser.parse_args()

def read_lc_from_hdf5(hdf5_file, aperture=4, fold_by='photref'):
    """
    Read light curve data from HDF5 file and organize it by channel or photref.
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
    data : dict
        Dictionary of DataFrames keyed by channel or photref
        ex: data = {0: DataFrame, 1: DataFrame, ...}
        Each DataFrame has columns: 'bjd', 'mag_epd', 'mag_magfit'
    """

    with h5py.File(hdf5_file, 'r') as f:
        channels = f[f'FrameInformation/ConfigurationIndex'][:]
        photrefs = f[f'AperturePhotometry/Aperture004/MagnitudeFitting/ConfigurationIndex'][:]
        
        bjds = f[f'SkyPosition/BJD'][:]
        
        mag_magfit = f[f'AperturePhotometry/Aperture{aperture:03}/MagnitudeFitting/Magnitude'][:]
        mag_epd = f[f'AperturePhotometry/Aperture{aperture:03}/EPD/Magnitude'][:]

        for mag_array in [mag_epd, mag_magfit]:
            mag_array[mag_array < -1e5] = np.nan  # Replace NaN (large negative) values with NaN
        
        data = {}

        if fold_by == 'channel':
            for chn in np.unique(channels):
                mask = channels == chn
                med_epd = np.nanmedian(mag_epd[mask])
                med_magfit = np.nanmedian(mag_magfit[mask])
                data[chn] = pd.DataFrame({
                    'bjd': bjds[mask],
                    'mag_epd': mag_epd[mask] - med_epd,
                    'mag_magfit': mag_magfit[mask] - med_magfit
                })

        elif fold_by == 'photref':
            for phref in np.unique(photrefs):
                mask = photrefs == phref
                med_epd = np.nanmedian(mag_epd[mask])
                med_magfit = np.nanmedian(mag_magfit[mask])
                data[phref] = pd.DataFrame({
                    'bjd': bjds[mask],
                    'mag_epd': mag_epd[mask] - med_epd,
                    'mag_magfit': mag_magfit[mask] - med_magfit
                })

    return data

def apply_binning(df, x_values, method, bin_size):
    """
    Apply binning to light curve data.
    
    Parameters:
    -----------
    df : DataFrame
        DataFrame with 'mag_epd', 'mag_magfit' columns
    x_values : array
        Time or phase values (bjd or phase)
    method : str
        'time' or 'phase'
    bin_size : float
        Size of bins (minutes for time, number of phase divisions for phase)
    
    Returns:
    --------
    binned_x : array
        Binned x values (bin centers)
    binned_mag_epd : array
        Binned magnitudes (epd)
    binned_mag_magfit : array
        Binned magnitudes (magfit)
    """
    if method == 'time':
        # bin_size is in minutes, convert to days
        bin_size_days = bin_size / (24 * 60)
        bins = np.arange(x_values.min(), x_values.max() + bin_size_days, bin_size_days)
    elif method == 'phase':
        # Divide 0-1 phase into bin_size divisions
        bins = np.linspace(0, 1, int(bin_size) + 1)
    
    # Create a temporary DataFrame for binning
    temp_df = df.copy()
    temp_df['x'] = x_values
    temp_df['bin'] = pd.cut(x_values, bins=bins)
    
    # Group by bins and compute median
    binned = temp_df.groupby('bin', observed=True).agg({
        'x': 'median',
        'mag_epd': 'median',
        'mag_magfit': 'median'
    }).dropna()
    
    return binned['x'].values, binned['mag_epd'].values, binned['mag_magfit'].values

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

def plot_light_curve(bjd, mag_epd, mag_magfit,
                     fold_by, plot_type,
                     object_name, chn_or_phref,
                     xlabel, period=None
                    ):
    """
    Plot light curve.
    """
    # plt.figure(figsize=(10, 5))
    plt.scatter(bjd, mag_epd, s=10, color='red', alpha=0.5, label='EPD')
    # plt.scatter(bjd, mag_magfit, s=10, color='blue', alpha=0.5, label='MagFit')
    
    plt.gca().invert_yaxis()  # Magnitude axis is inverted
    plt.xlabel(xlabel, fontdict={"fontweight":"bold", 'fontsize':14})
    plt.ylabel("Magnitude", fontdict={"fontweight":"bold", 'fontsize':14})
    plt.ylim(1.5, -0.6)
    plt.grid(True)
    # plt.legend()
    # plt.savefig(f'{object_name}_Chn_or_photref_{chn_or_phref}',dpi=400)
    if plot_type == 'individual':
        plt.title(f'{object_name} - {fold_by}: {chn_or_phref}', fontdict={"fontweight":"bold", 'fontsize':16})
        plt.show()

def main():

    args = parse_arguments()
    hdf5_file = args.hdf5_file
    object_name = hdf5_file.split('/')[-1].split('.')[0]
    aperture = args.aperture
    period = args.period
    epoch = args.epoch
    fold_by = args.fold_by
    plot_type = args.plot_type
    binning_method = args.binning[0] if args.binning else None
    binning_size = float(args.binning[1]) if args.binning else None

    data = read_lc_from_hdf5(hdf5_file, aperture, fold_by)

    # If combined plot type with binning, merge all dataframes first
    if plot_type == 'combined' and binning_method:
        all_bjds = []
        all_mag_epd = []
        all_mag_magfit = []
        
        for chn_or_phref in sorted(data.keys()):
            df = data[chn_or_phref]
            all_bjds.append(df['bjd'].values)  # combination of separate arrays
            all_mag_epd.append(df['mag_epd'].values)
            all_mag_magfit.append(df['mag_magfit'].values)
        
        x_values = np.concatenate(all_bjds)  # make one big array of bjd values
        mag_epd = np.concatenate(all_mag_epd)
        mag_magfit = np.concatenate(all_mag_magfit)
        
        # Create combined DataFrame with all columns and sort by bjd
        combined_df = pd.DataFrame({
            'bjd': x_values,
            'mag_epd': mag_epd,
            'mag_magfit': mag_magfit
        })
        combined_df = combined_df.sort_values('bjd').reset_index(drop=True)
        
        # Apply phase folding if needed
        if period:
            x_plot = calculate_phase(combined_df['bjd'].values, period, epoch)
            xlabel = "Phase"
        else:
            x_plot = x_values - 2450000
            xlabel = "BJD - 2450000"
        
        # Apply binning
        x_plot, mag_epd_binned, mag_magfit_binned = apply_binning(
            combined_df, x_plot, binning_method, binning_size
        )
        
        plot_light_curve(
            x_plot, mag_epd_binned, mag_magfit_binned,
            fold_by, plot_type,
            object_name, 'combined',
            xlabel=xlabel, period=period
        )
    else:
        # Original logic for individual or combined without binning
        for chn_or_phref in sorted(data.keys()):
            df = data[chn_or_phref]
            x_values = df['bjd'].values
            mag_epd = df['mag_epd'].values
            mag_magfit = df['mag_magfit'].values
        
            if period:
                phase = calculate_phase(x_values, period, epoch)
                plot_light_curve(
                    phase, mag_epd, mag_magfit,
                    fold_by, plot_type,
                    object_name, chn_or_phref,
                    xlabel="Phase", period=period
                    )
            else:
                # Plot regular light curve
                bjd = x_values - 2450000
                plot_light_curve(
                    bjd, mag_epd, mag_magfit,
                    fold_by, plot_type,
                    object_name, chn_or_phref,
                    xlabel="BJD - 2450000"
                    )
    if plot_type == 'combined':
        plt.title(
            f'{object_name} - folded by: {fold_by}',
            fontdict={"fontweight":"bold", 'fontsize':16}
            )
        plt.ylim(1, -0.5)
        plt.savefig(f'{object_name}_folded_by_{fold_by}_combined', dpi=400)
        plt.show()


if __name__ == "__main__":
    main()