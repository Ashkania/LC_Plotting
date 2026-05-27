#!/usr/bin/env python3

"""
    Plot the LC from its hdf5 file:
    
    1. Read the LC data from the hdf5 file. It includes:
       BJD, mag_epd, mag_magfit, channels index and photrefs index.
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
from astroquery.gaia import Gaia
import sys

def parse_arguments():
    
    parser = ArgumentParser(
        description="Plot Light Curve from HDF5 file"
        )
    parser.add_argument(
        "lc_file",
        type=str,
        help="Path to the light curve HDF5 file"
        )
    parser.add_argument(
        "--mag-types",
        nargs='+',
        default=['magfit'],
        help="Which magnitude type to plot (magfit, epd, or tfa or multiple)"
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
        "--text",
        type=str,
        default=None,
        help="Printed text in the plot"
    )
    parser.add_argument(
        "--sep-by",
        type=str,
        choices=['channel', 'photref'],
        default='photref',
        help="Whether to separate data points by 'channel' or 'photref(=channel and camera)'" \
        "(default: 'photref')"
    )
    parser.add_argument(
        "--mode",
        choices=['single', 'folded'],
        default='single',
        help="Plot each channel/photrefs(usually based on channel, camera, exposure time, etc.)" \
        "separately or fold all channels/photrefs on top of each other" \
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
        "--GAIA-plot",
        default=False,
        action='store_true',
        help="Plot GAIA light curve for the same object on top of the PANOPTES data"
    )
    parser.add_argument(
        "--save-or-show-plot",
        default='show',
        choices=['save', 'show'],
        help="Whether to save the plot as a PNG file or display it (default: 'show')"
    )
    # parser.add_argument(
    #     "--folded-color",
    #     default=False,
    #     action='store_true',
    #     help="When used with --mode folded, plot each photref/channel in a separate color on the same folded figure"
    # )

    return parser.parse_args()

def read_lc(hdf5_file, aperture=4, sep_by='photref', mag_types=['magfit']):
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
    photref_dict : dict
        Dictionary mapping photref index to processed name
    """

    MAG_CONFIG = {
        'magfit': ('MagnitudeFitting/Magnitude', 'mag_magfit'),
        'epd':    ('EPD/Magnitude',              'mag_epd'),
        'tfa':    ('TFA/Magnitude',              'mag_tfa'),
    }

    with h5py.File(hdf5_file, 'r') as f:
        channels = f[f'FrameInformation/ConfigurationIndex'][:]
        photrefs = f[f'AperturePhotometry/Aperture004/MagnitudeFitting/ConfigurationIndex'][:]

        channel_photref_mapping = {
            ch: np.unique(photrefs[channels == ch])
            for ch in np.unique(channels)
        }
        
        bjds = f[f'SkyPosition/BJD'][:]
        
        mag_arrays = {}
        for mag in mag_types:
            hdf5_subpath, col_name = MAG_CONFIG[mag]
            try:
                arr = f[f'AperturePhotometry/Aperture{aperture:03}/{hdf5_subpath}'][:]
                arr[arr < -1e5] = np.nan
                mag_arrays[col_name] = arr
            except KeyError:
                print(f"KeyError: {hdf5_subpath} not found in {hdf5_file}")
                sys.exit(1)

        data = {}

        for phref in np.unique(photrefs):
            mask = photrefs == phref
            df_data = {'bjd': bjds[mask]}
            for col_name, arr in mag_arrays.items():
                df_data[col_name] = arr[mask] - np.nanmedian(arr[mask])
            data[phref] = pd.DataFrame(df_data)

        # --- Regroup if needed ---
        photref_dict = None

        if sep_by == 'channel':
            data = {
                chn: pd.concat([data[phref] for phref in phrefs], ignore_index=True)
                for chn, phrefs in channel_photref_mapping.items()
            }

        elif sep_by == 'photref':
            photrefsnames = f[f'AperturePhotometry/Aperture000/MagnitudeFitting/Configuration/SinglePhotometricReference'][:]
            photref_dict = {}
            for i, name in enumerate(photrefsnames):
                name = name.decode('utf-8')
                part1 = name[name.find('DR')+3:name.find('DR')+7]
                # Find second occurrence of 'PAN' (not the one in 'PANOPTES')
                second_pan = name.find('PAN', name.find('PAN') + 1)
                part2 = name[second_pan - 3:second_pan - 1]
                part3 = name[-5:-3]
                photref_dict[i] = f"{part1}_{part2}_{part3}"
                
    return data, photref_dict

def apply_binning(df, x_values, method, bin_size):   ############# edit with new mag_cols #############
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
    
    return binned['x'].values, binned['mag_epd'].values, binned['mag_magfit'].values ####################3 edit #######################################

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

def plot_TESS_lc(gaia_id, period, epoch, text=None):
    """
    Plot TESS light curve for the given Gaia ID, period, and epoch.
    Downloads the most recent SPOC light curve.
    """
    gaia_id = 'Gaia DR3 ' + gaia_id
    search = lk.search_lightcurve(gaia_id, mission="TESS", author="SPOC")
    if search is None or len(search) == 0:
        print(f"No TESS SPOC data found for {gaia_id}, trying any TESS data...")
        search = lk.search_lightcurve(gaia_id, mission="TESS")
        if search is None or len(search) == 0:
            print(f"No TESS data found for {gaia_id}. Skipping TESS plot.")
            return
    else:
        print(search)
    # Sort by t_min descending to get the most recent
    search = search[np.argsort(search.table['t_min'])[::-1]]
    print(f"Downloading most recent SPOC LC: {search[0]}")
    lc = search[0].download()
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

    plt.scatter(phase, mag_tess, s=1, color='green', alpha=0.7, label="TESS", zorder=2)
    plt.text(
        0.95, 0.04,
        f'TESS Mag: {text}', transform=plt.gca().transAxes,
        fontsize=10, ha='right', va='bottom'
        )

def plot_Gaia_lc(gaia_id, period, epoch, text=None):
    datalink = Gaia.load_data(
        ids=gaia_id,
        data_release="Gaia DR3",
        retrieval_type="EPOCH_PHOTOMETRY",
        format="votable"
    )
    table = datalink[f'EPOCH_PHOTOMETRY-Gaia DR3 {gaia_id}.xml'][0].to_table().to_pandas()

    t_gaia_bjd = np.array(table['g_transit_time']) + 2455197.5
    mag = np.array(table['g_transit_mag'])
    mask = np.isfinite(t_gaia_bjd) & np.isfinite(mag)
    t_gaia_bjd, mag = t_gaia_bjd[mask], mag[mask]

    mag = mag - np.mean(mag)
    phase = ((t_gaia_bjd - epoch) / period) % 1
    plt.scatter(phase, mag, s=70, marker='x', c='r', label="GAIA", zorder=2)
    plt.text(
        0.95, 0.04,
        f'{text}', transform=plt.gca().transAxes,
        fontsize=10, ha='right', va='bottom'
        )
def plot_lc(X, *mag_arrays, mag_cols, xlabel):
    """
    Plot light curve.

    Parameters:
    -----------
    X : array
        Time or phase axis values
    mag_arrays : arrays
        Magnitude arrays, one per entry in mag_cols
    mag_cols : list of str
        Column names corresponding to mag_arrays (e.g. ['mag_magfit', 'mag_epd'])
    xlabel : str
        X-axis label
    """
    MAG_DISPLAY = {
        'mag_magfit': ('MagFit', 'red'),
        'mag_epd':    ('EPD',    'blue'),
        'mag_tfa':    ('TFA',    'green'),
    }

    for i, (col, arr) in enumerate(zip(mag_cols, mag_arrays)):
        display_name, color = MAG_DISPLAY[col]
        plt.scatter(
            X, arr, s=10, color=color, alpha=0.5,
            label=display_name, zorder=len(mag_cols) - i,
        )

    plt.gca().invert_yaxis()  # Magnitude axis is inverted
    plt.xlabel(xlabel, fontdict={"fontweight":"bold", 'fontsize':14})
    plt.ylabel("Magnitude", fontdict={"fontweight":"bold", 'fontsize':14})
    plt.grid(True)
    plt.legend()

    all_mags = np.concatenate(mag_arrays)
    valid = ~np.isnan(all_mags)
    if valid.any():
        q1, q3 = np.percentile(all_mags[valid], [10, 90])
        iqr = q3 - q1
        plt.ylim(q3 + 5 * iqr, q1 - 3 * iqr)
    else:
        plt.ylim(0.7, -0.5)  # fallback

def main():

    args = parse_arguments()
    lc_file = args.lc_file
    mag_types = args.mag_types
    gaia_id = lc_file.split('/')[-1].split('.')[0].split('_')[1]
    aperture = args.aperture
    period = args.period
    epoch = args.epoch
    text = args.text
    sep_by = args.sep_by
    mode = args.mode
    TESS_plot = args.TESS_plot
    GAIA_plot = args.GAIA_plot
    binning_method = args.binning[0] if args.binning else None
    binning_size = float(args.binning[1]) if args.binning else None
    save_or_show = args.save_or_show_plot

    data, photref_dict = read_lc(lc_file, aperture, sep_by, mag_types)
    if args.selected:
        selected = args.selected
        data = {k: v for k, v in data.items() if k in selected}
        if photref_dict:
            photref_dict = {k: v for k, v in photref_dict.items() if k in selected}

    # Derive mag columns from data (anything that isn't 'bjd')
    sample_df = next(iter(data.values()))
    mag_cols = [c for c in sample_df.columns if c != 'bjd']

    if mode == 'folded':
        # if args.folded_color:
        #     cmap = plt.get_cmap('tab20')
        #     xlabel = "Phase" if period else "BJD - 2450000"

        #     for idx, chn_or_phref in enumerate(sorted(data.keys())):
        #         df = data[chn_or_phref]
        #         bjd_or_phase = (
        #             calculate_phase(df['bjd'].values, period, epoch)
        #             if period else df['bjd'].values - 2450000
        #         )

        #         if binning_method:
        #             x_plot, *mag_arrays = apply_binning(
        #                 df, bjd_or_phase, binning_method, binning_size
        #             )
        #         else:
        #             x_plot = bjd_or_phase
        #             mag_arrays = [df[col].values for col in mag_cols]

        #         title_part = (
        #             f"{chn_or_phref}) {photref_dict[chn_or_phref]}"
        #             if photref_dict and chn_or_phref in photref_dict
        #             else str(chn_or_phref)
        #         )

        #         plot_lc(
        #             x_plot,
        #             *mag_arrays, mag_cols,
        #             xlabel,
        #             color=cmap(idx % cmap.N),
        #             label=title_part,
        #             label_magfit_only=True
        #         )

        #     if period:
        #         if TESS_plot:
        #             plot_TESS_lc(gaia_id, period, epoch, text)
        #         elif GAIA_plot:
        #             plot_Gaia_lc(gaia_id, period, epoch, text)

        #     plt.title(
        #         f'Gaia {gaia_id} - separated by: {sep_by}',
        #         fontdict={"fontweight":"bold", 'fontsize':12}
        #     )
        #     plt.legend(fontsize=8, loc='best')
        #     if save_or_show == 'save':
        #         plt.savefig(f'Gaia_{gaia_id}_sepby_{sep_by}_folded_color_aperture_{aperture}', dpi=400)
        #     elif save_or_show == 'show':
        #         plt.show()
        # else:
        combined_df = pd.concat(
            [data[k] for k in sorted(data.keys())],
            ignore_index=True
        )

        bjd_or_phase = (
            calculate_phase(combined_df['bjd'].values, period, epoch)
            if period else combined_df['bjd'].values - 2450000
        )
        xlabel = "Phase" if period else "BJD - 2450000"
        
        if period:
            if TESS_plot:
                plot_TESS_lc(gaia_id, period, epoch, text)
            elif GAIA_plot:
                plot_Gaia_lc(gaia_id, period, epoch, text)
        
        if binning_method:
            x_plot, *mag_arrays = apply_binning(
                combined_df, bjd_or_phase, binning_method, binning_size, mag_cols
            )
        else:
            x_plot = bjd_or_phase
            mag_arrays = [combined_df[col].values for col in mag_cols]
        
        plot_lc(x_plot, *mag_arrays, mag_cols=mag_cols, xlabel=xlabel)
        plt.title(
            f'Gaia {gaia_id} - separated by: {sep_by}',
            fontdict={"fontweight":"bold", 'fontsize':12}
        )
        if save_or_show == 'save':
            plt.savefig(f'Gaia_{gaia_id}_sepby_{sep_by}_folded_aperture_{aperture}', dpi=400)
        elif save_or_show == 'show':
            plt.show()

    elif mode == 'single':
        for chn_or_phref in sorted(data.keys()):
            df = data[chn_or_phref]
            
            bjd_or_phase = (
                calculate_phase(df['bjd'].values, period, epoch)
                if period else df['bjd'].values - 2450000
            )
            xlabel = "Phase" if period else "BJD - 2450000"

            if period:
                xlabel = "Phase"
                if TESS_plot:
                    plot_TESS_lc(gaia_id, period, epoch, text)
                elif GAIA_plot:
                    plot_Gaia_lc(gaia_id, period, epoch, text)
            
            if binning_method:
                x_plot, *mag_arrays = apply_binning(
                    df, bjd_or_phase, binning_method, binning_size, mag_cols
                )
            else:
                x_plot = bjd_or_phase
                mag_arrays = [df[col].values for col in mag_cols]
            
            plot_lc(x_plot, *mag_arrays, mag_cols=mag_cols, xlabel=xlabel)
            title_part = (
                f"{chn_or_phref}) {photref_dict[chn_or_phref]}"
                if photref_dict and chn_or_phref in photref_dict
                else str(chn_or_phref)
            )
            plt.title(f'{gaia_id} - {sep_by}: {title_part}', fontdict={"fontweight":"bold", 'fontsize':12})
            if save_or_show == 'save':
                plt.savefig(f'Gaia_{gaia_id}_sepby_{sep_by}_{title_part}_aperture_{aperture}', dpi=400)
            elif save_or_show == 'show':
                plt.show()
            plt.clf()

if __name__ == "__main__":
    main()