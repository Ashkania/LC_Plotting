#!/usr/bin/env python3

"""
    Plot the LC from its hdf5 file:
    
    1. Read the LC data from the hdf5 file. It includes:
       bjd, mag_epd, mag_magfit, channels index and photrefs index.
       Based on the folding method (by channel or photref), data dictionary
       is constructed to contain separate DataFrames of bjd, epd, magift for each channel/photref.
       Like:
         data = {
              0: DataFrame for channel/channel_camera 0,
              1: DataFrame for channel/channel_camera 1,
              ...
         }
    We then could plot one for each channel/photref separately (single mode)
    or combine all channels/photrefs into one plot. (folde mode)
    
    2. Optionally apply binning to the data (by time or phase).
    3. Optionally fold the light curve by a given period and epoch.
    4. Plot the light curve, either individually for each channel/photref or combined into one figure.
    example:
    python scripts/plotting.py LC/EB_Villa/GDR3_3337025761759515776.h5 
        --aperture 4 --period 1.7801492 --epoch 1469.56553 --sep-by 'photref'
        --mode 'single' --binning 'time' 10
"""

#TODO:
# 1. last part of the title in the photref case is kinda hardcoded
# 2. folded mode with binning is fine, but without binning, but if
# binning is not applied, the folded plot is not correct.

import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import lightkurve as lk
from argparse import ArgumentParser

def parse_arguments():
    
    parser = ArgumentParser(
        description="Plot Light Curve from HDF5 file"
        )
    parser.add_argument(
        "lc_file",
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
        "--tmag",
        type=float,
        default=None,
        help="TESS magnitude of the star (only for TESS plot annotation)"
    )
    parser.add_argument(
        "--sep-by",
        type=str,
        choices=['channel', 'photref', 'camera'],
        default='photref',
        help="Whether to separate by 'channel' or 'photref(=channel and camera)'" \
        "camera is not implemented yet. (default: 'photref')"
    )
    parser.add_argument(
        "--mode",
        choices=['single', 'folded'],
        default='single',
        help="Plot each channel/photrefs(=channel and camera) separately or" \
        "fold all channels/photrefs(=channel and cameras) on top of each other" \
        "(default: 'single')"
    )
    parser.add_argument(
        "--binning",
        nargs=2,
        default=None,
        metavar=('METHOD', 'SIZE'),
        help="Binning method (time/phase) and size (minutes for time, points for phase)" \
            "(optional, e.g. --binning time 10 or --binning phase 100)"
    )
    parser.add_argument(
        "--selected",
        nargs='*',
        type=int,
        default=None,
        help="List of photref/channel indices to plot (default: all)"
    )
    parser.add_argument(
        "--TESS-plot",
        default=False,
        action='store_true',
        help="Plot TESS light curve for the same object on top of the PANOPTES data"
    )
    parser.add_argument(
        "--save-plot",
        default=False,
        action='store_true',
        help="Save the plot as a PNG file"
    )

    return parser.parse_args()

def read_lc(hdf5_file, aperture=4, sep_by='photref'):
    """
    Read light curve data from HDF5 file and organize it by channel or photref.
     Parameters:
    -----------
    hdf5_file : str
        Path to the HDF5 file containing the light curve data
    aperture : int, optional
        Aperture number to read (default: 4)
    sep_by : str, optional
        Whether to separate by 'channel' or 'photref' (default: 'photref')
    Returns:
    --------
    data : dict
        Dictionary of DataFrames keyed by channel or photref
        ex: data = {0: DataFrame, 1: DataFrame, ...}
        Each DataFrame has columns: 'bjd', 'mag_epd', 'mag_magfit'
    photref_dict : dict or None
        Dictionary mapping photref index to processed name, only if sep_by='photref'
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
        photref_dict = {}

        if sep_by == 'channel':
            for chn in np.unique(channels):
                mask = channels == chn
                med_epd = np.nanmedian(mag_epd[mask])
                med_magfit = np.nanmedian(mag_magfit[mask])
                data[chn] = pd.DataFrame({
                    'bjd': bjds[mask],
                    'mag_epd': mag_epd[mask] - med_epd,
                    'mag_magfit': mag_magfit[mask] - med_magfit
                })

        elif sep_by == 'photref':
            photrefsnames = f[f'AperturePhotometry/Aperture000/MagnitudeFitting/Configuration/SinglePhotometricReference'][:]
            for i, name in enumerate(photrefsnames):
                name = name.decode('utf-8')
                part1 = name[name.find('DR')+3:name.find('DR')+7]
                # Find second occurrence of 'PAN' (not the one in 'PANOPTES')
                second_pan = name.find('PAN', name.find('PAN') + 1)
                part2 = name[second_pan - 3:second_pan - 1]
                part3 = name[-5:-3]
                photref_dict[i] = f"{part1}_{part2}_{part3}"
            
            for phref in np.unique(photrefs):
            # for phref in [12, 28, 20]:
                mask = photrefs == phref
                med_epd = np.nanmedian(mag_epd[mask])
                med_magfit = np.nanmedian(mag_magfit[mask])
                data[phref] = pd.DataFrame({
                    'bjd': bjds[mask],
                    'mag_epd': mag_epd[mask] - med_epd,
                    'mag_magfit': mag_magfit[mask] - med_magfit
                })
        # TODO: find the best way to separate by camera
        # elif sep_by == 'camera':
            # for cam in ...:
            #     for phref in ...:
            #     mask = (
            #         (photrefs == phref0) or
            #         (photrefs == phref1) or
            #         (photrefs == phref2) or
            #         (photrefs == phref3)
            #     )
            #     med_epd = np.nanmedian(mag_epd[mask])
            #     med_magfit = np.nanmedian(mag_magfit[mask])
            #     data[(chn, phref)] = pd.DataFrame({
            #         'bjd': bjds[mask],
            # 'mag_epd': mag_epd[mask] - med_epd,
            # 'mag_magfit': mag_magfit[mask] - med_magfit
                # })

    return data, photref_dict

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

def plot_TESS_lc(gaia_id, period, epoch, tmag=None):
    """
    Plot TESS light curve for the given Gaia ID, period, and epoch.
    Should we download all data or the first as now is enough?
    """
    gaia_id = 'Gaia DR3 ' + gaia_id
    search = lk.search_lightcurve(gaia_id, mission="TESS")
    # print(search)

    lc = search.download()
    lc = lc.remove_nans().normalize()
    lc = lc.flatten(window_length=401)

    time_tess = lc.time.value + 2457000
    # print(time_tess)
    # print(lc.time.jd)
    # print(lc.time.utc.iso)
    # print(lc.time.tdb.iso)
    # print(type(lc.time.format))
    phase = ((time_tess - epoch) / period) % 1
    flux_tess = lc.flux.value
    mag_tess = -2.5 * np.log10(flux_tess)

    plt.scatter(phase, mag_tess, s=1, color='green', alpha=0.7, label="TESS")
    plt.text(
        0.95, 0.04,
        f'TESS Mag: {tmag}', transform=plt.gca().transAxes,
        fontsize=10, ha='right', va='bottom'
        )

def plot_lc(bjd, mag_epd, mag_magfit,xlabel):
    """
    Plot light curve.
    """
    # plt.figure(figsize=(10, 5))
    plt.scatter(bjd, mag_magfit, s=10, color='blue', alpha=0.5, label='MagFit')
    plt.scatter(bjd, mag_epd, s=10, color='red', alpha=0.5, label='EPD')

    plt.gca().invert_yaxis()  # Magnitude axis is inverted
    plt.xlabel(xlabel, fontdict={"fontweight":"bold", 'fontsize':14})
    plt.ylabel("Magnitude", fontdict={"fontweight":"bold", 'fontsize':14})
    plt.ylim(0.7, -0.5)
    # plt.xlim(9900, 9907)
    plt.grid(True)
    # plt.legend()
    # plt.savefig(f'{object_name}_Chn_or_photref_{chn_or_phref}',dpi=400)
    # if mode == 'single':
    #     if photref_dict and chn_or_phref in photref_dict:
    #         title_part = f"{chn_or_phref}: {photref_dict[chn_or_phref]}"
    #     else:
    #         title_part = chn_or_phref
    #     plt.title(f'{object_name} - {sep_by}: {title_part}', fontdict={"fontweight":"bold", 'fontsize':12})
    #     plt.show()

def main():

    args = parse_arguments()
    lc_file = args.lc_file
    gaia_id = lc_file.split('/')[-1].split('.')[0].split('_')[1]
    aperture = args.aperture
    period = args.period
    epoch = args.epoch
    tmag = args.tmag
    sep_by = args.sep_by
    mode = args.mode
    TESS_plot = args.TESS_plot
    binning_method = args.binning[0] if args.binning else None
    binning_size = float(args.binning[1]) if args.binning else None

    data, photref_dict = read_lc(lc_file, aperture, sep_by)

    if args.selected:
        selected = args.selected
        data = {k: v for k, v in data.items() if k in selected}
        if photref_dict:
            photref_dict = {k: v for k, v in photref_dict.items() if k in selected}

    if mode == 'folded':
        all_bjds_list_of_lists = []
        all_mag_epd_list_of_lists = []
        all_mag_magfit_list_of_lists = []
        
        for chn_or_phref in sorted(data.keys()):
            df = data[chn_or_phref]
            all_bjds_list_of_lists.append(df['bjd'].values)  # combination of separate arrays
            all_mag_epd_list_of_lists.append(df['mag_epd'].values)
            all_mag_magfit_list_of_lists.append(df['mag_magfit'].values)
        
        all_bjds = np.concatenate(all_bjds_list_of_lists)  # make one big array of bjd values
        all_mag_epd = np.concatenate(all_mag_epd_list_of_lists)
        all_mag_magfit = np.concatenate(all_mag_magfit_list_of_lists)
        
        combined_df = pd.DataFrame({
            'bjd': all_bjds,
            'mag_epd': all_mag_epd,
            'mag_magfit': all_mag_magfit
        })

        if period:
            bjd_or_phase = calculate_phase(combined_df['bjd'].values, period, epoch)
            xlabel = "Phase"
            if TESS_plot:
                plot_TESS_lc(gaia_id, period, epoch, tmag)
        else:
            bjd_or_phase = combined_df['bjd'].values - 2450000
            xlabel = "BJD - 2450000"
        
        if binning_method:
            x_plot, mag_epd, mag_magfit = apply_binning(
                combined_df, bjd_or_phase, binning_method, binning_size
            )
        else:
            x_plot = bjd_or_phase
            mag_epd = combined_df['mag_epd'].values
            mag_magfit = combined_df['mag_magfit'].values
        
        plot_lc(x_plot, mag_epd, mag_magfit, xlabel)
        plt.title(
            f'Gaia {gaia_id} - separated by: {sep_by}',
            fontdict={"fontweight":"bold", 'fontsize':12}
            )
        plt.ylim(0.7, -0.5)
        if args.save_plot:
            plt.savefig(f'Gaia_{gaia_id}_sepby_{sep_by}_folded_aperture_{aperture}', dpi=400)
        plt.show()


    elif mode == 'single':
        for chn_or_phref in sorted(data.keys()):
            df = data[chn_or_phref]
            
            if period:
                bjd_or_phase = calculate_phase(df['bjd'].values, period, epoch)
                xlabel = "Phase"
                if TESS_plot:
                    plot_TESS_lc(gaia_id, period, epoch, tmag)
            else:
                bjd_or_phase = df['bjd'].values - 2450000
                xlabel = "BJD - 2450000"
            
            if binning_method:
                x_plot, mag_epd, mag_magfit = apply_binning(
                    df, bjd_or_phase, binning_method, binning_size)
            else:
                x_plot = bjd_or_phase
                mag_epd = df['mag_epd'].values
                mag_magfit = df['mag_magfit'].values
            
            plot_lc(x_plot, mag_epd, mag_magfit,xlabel)
            if photref_dict and chn_or_phref in photref_dict:
                title_part = f"{chn_or_phref}) {photref_dict[chn_or_phref]}"
            else:
                title_part = chn_or_phref
            plt.title(f'{gaia_id} - {sep_by}: {title_part}', fontdict={"fontweight":"bold", 'fontsize':12})
            if args.save_plot:
                plt.savefig(f'Gaia_{gaia_id}_sepby_{sep_by}_{title_part}_aperture_{aperture}', dpi=400)
            plt.show()

if __name__ == "__main__":
    main()