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
    parser.add_argument(
        "--ylim",
        nargs=2,
        type=float,
        default=None,
        metavar=('YMIN', 'YMAX'),
        help="Y-axis limits for the plot (optional, e.g. --ylim -0.5 0.5)"
    )

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
    name_dict : dict or None
        Mapping from key in `data` to a human-readable name.
        - sep_by='photref': processed photref name
        - sep_by='channel': channel name read from FrameInformation/ColorChannel (e.g. 'R0','G1','G0','B0')
    channel_color_dict : dict or None
        Mapping from channel index to matplotlib color (only when sep_by='channel').
    """

    MAG_CONFIG = {
        'magfit': ('MagnitudeFitting/Magnitude', 'mag_magfit'),
        'epd':    ('EPD/Magnitude',              'mag_epd'),
        'tfa':    ('TFA/Magnitude',              'mag_tfa'),
    }
    CHANNEL_COLOR_MAP = {'R': 'red', 'G': 'green', 'B': 'blue'}

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
        name_dict = None
        channel_color_dict = None

        if sep_by == 'channel':
            data = {
                chn: pd.concat([data[phref] for phref in phrefs], ignore_index=True)
                for chn, phrefs in channel_photref_mapping.items()
            }
            try:
                channel_names = f['FrameInformation/Configuration/ColorChannel'][:]
                name_dict = {}
                channel_color_dict = {}
                for i, name in enumerate(channel_names):
                    if isinstance(name, bytes):
                        name = name.decode('utf-8')
                    name_dict[i] = name
                    channel_color_dict[i] = CHANNEL_COLOR_MAP.get(name[:1], 'gray')
            except KeyError:
                print("Warning: FrameInformation/ColorChannel not found; channel coloring disabled.")

        elif sep_by == 'photref':
            photrefsnames = f[f'AperturePhotometry/Aperture000/MagnitudeFitting/Configuration/SinglePhotometricReference'][:]
            name_dict = {}
            for i, name in enumerate(photrefsnames):
                name = name.decode('utf-8')
                part1 = name[name.find('DR')+3:name.find('DR')+7]
                # Find second occurrence of 'PAN' (not the one in 'PANOPTES')
                second_pan = name.find('PAN', name.find('PAN') + 1)
                part2 = name[second_pan - 3:second_pan - 1]
                part3 = name[-5:-3]
                name_dict[i] = f"{part1}_{part2}_{part3}"

    return data, name_dict, channel_color_dict

def apply_binning(df, x_values, method, bin_size, mag_cols):
    """
    Apply binning to light curve data.

    Parameters:
    -----------
    df : DataFrame
        DataFrame containing the magnitude columns listed in mag_cols
    x_values : array
        Time or phase values (bjd or phase)
    method : str
        'time' or 'phase'
    bin_size : float
        Size of bins (minutes for time, number of phase divisions for phase)
    mag_cols : list of str
        Magnitude column names to bin (e.g. ['mag_magfit', 'mag_epd', 'mag_tfa'])

    Returns:
    --------
    binned_x : array
        Binned x values (bin medians)
    *binned_mags : arrays
        One binned magnitude array per entry in mag_cols, in the same order
    """
    if method == 'time':
        # bin_size is in minutes, convert to days
        bin_size_days = bin_size / (24 * 60)
        bins = np.arange(x_values.min(), x_values.max() + bin_size_days, bin_size_days)
    elif method == 'phase':
        # Divide 0-1 phase into bin_size divisions
        bins = np.linspace(0, 1, int(bin_size) + 1)

    temp_df = df.copy()
    temp_df['x'] = x_values
    temp_df['bin'] = pd.cut(x_values, bins=bins)

    agg_spec = {'x': 'median', **{col: 'median' for col in mag_cols}}
    # agg_sepc = {'x': 'median', 'mag_magfit': 'median', ...}
    binned = temp_df.groupby('bin', observed=True).agg(agg_spec).dropna()

    return (binned['x'].values, *(binned[col].values for col in mag_cols))

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

def calculate_magnitude_offset(lc_phase, lc_magnitudes, tess_phase, tess_magnitudes):
    """
    Calculate the LC-to-TESS offset using their overlapping phase range.
    """

    lc_phase = np.asarray(lc_phase, dtype=float)
    lc_magnitudes = np.asarray(lc_magnitudes, dtype=float)
    tess_phase = np.asarray(tess_phase, dtype=float)
    tess_magnitudes = np.asarray(tess_magnitudes, dtype=float)

    lc_valid = np.isfinite(lc_phase) & np.isfinite(lc_magnitudes)
    tess_valid = np.isfinite(tess_phase) & np.isfinite(tess_magnitudes)

    if not lc_valid.any() or not tess_valid.any():
        return 0.0

    common_min = max(lc_phase[lc_valid].min(), tess_phase[tess_valid].min())
    common_max = min(lc_phase[lc_valid].max(), tess_phase[tess_valid].max())
    if common_min > common_max:
        print("Warning: LC and TESS phase ranges do not overlap; no offset applied.")
        return 0.0

    lc_in_range = lc_valid & (lc_phase >= common_min) & (lc_phase <= common_max)
    tess_in_range = tess_valid & (tess_phase >= common_min) & (tess_phase <= common_max)
    if not lc_in_range.any() or not tess_in_range.any():
        return 0.0

    return (
        np.nanmedian(lc_magnitudes[lc_in_range])
        - np.nanmedian(tess_magnitudes[tess_in_range])
    )

def plot_TESS_lc(gaia_id, period, epoch, lc_phase, lc_magnitudes, text=None):
    """
    Plot TESS light curve for the given Gaia ID, period, and epoch.
    Downloads and plots all available SPOC light curves.
    """
    gaia_id = 'Gaia DR3 ' + gaia_id
    search = lk.search_lightcurve(gaia_id, mission="TESS", author="SPOC")
    if search is None or len(search) == 0:
        print(f"No TESS SPOC data found for {gaia_id}, trying any TESS data...")
        search = lk.search_lightcurve(gaia_id, mission="TESS")
        if search is None or len(search) == 0:
            print(f"No TESS data found for {gaia_id}. Skipping TESS plot.")
            return 0.0
    else:
        print(search)
    print(f"Downloading {len(search)} TESS light curves...")
    light_curves = search.download_all()
    tess_phase_arrays = []
    tess_mag_arrays = []

    for index, lc in enumerate(light_curves):
        lc = lc.remove_nans().normalize()
        lc = lc.flatten(window_length=401)

        time_tess = lc.time.value + 2457000
        phase = ((time_tess - epoch) / period) % 1
        flux_tess = lc.flux.value
        mag_tess = -2.5 * np.log10(flux_tess)
        tess_phase_arrays.append(phase)
        tess_mag_arrays.append(mag_tess)

    tess_phase = np.concatenate(tess_phase_arrays) if tess_phase_arrays else np.array([])
    tess_magnitudes = np.concatenate(tess_mag_arrays) if tess_mag_arrays else np.array([])
    magnitude_offset = calculate_magnitude_offset(
        lc_phase, lc_magnitudes, tess_phase, tess_magnitudes
    )

    for index, (phase, mag_tess) in enumerate(zip(tess_phase_arrays, tess_mag_arrays)):

        plt.scatter(
            phase,
            mag_tess,
            s=1,
            color='blue',
            alpha=0.5,
            label="TESS" if index == 0 else None,
            zorder=2,
        )

    plt.text(
        0.95, 0.04,
        f'TESS Mag: {text}', transform=plt.gca().transAxes,
        fontsize=10, ha='right', va='bottom'
        )
    return magnitude_offset

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
    plt.scatter(phase, mag, s=40, color='green', label="GAIA", zorder=2)
    plt.text(
        0.95, 0.04,
        f'{text}', transform=plt.gca().transAxes,
        fontsize=10, ha='right', va='bottom'
        )



def _scatter_lc(X, *mag_arrays, mag_cols, color=None, label=None):
    """
    Scatter magnitude arrays onto the current axes.

    If `color` is given it overrides the per-mag-type default color
    (used when coloring by channel). If `label` is given it overrides
    the per-mag-type label and is applied only to the first mag array
    so the legend gets one entry per group.
    """

    MAG_DISPLAY = {
    'mag_magfit': ('MagFit', 'red'),
    'mag_epd':    ('EPD',    'blue'),
    'mag_tfa':    ('TFA',    'green'),
}
    
    for i, (col, arr) in enumerate(zip(mag_cols, mag_arrays)):
        default_name, default_color = MAG_DISPLAY[col]
        plt.scatter(
            X, arr, s=10,
            color=color if color is not None else default_color,
            alpha=0.8,
            label=(label if i == 0 else None) if label is not None else default_name,
            zorder=len(mag_cols) - i,
        )


def _finalize_lc_plot(xlabel, all_mags):
    plt.gca().invert_yaxis()  # Magnitude axis is inverted
    plt.xlabel(xlabel, fontdict={"fontweight":"bold", 'fontsize':14})
    plt.ylabel("Magnitude", fontdict={"fontweight":"bold", 'fontsize':14})
    plt.grid(True)
    plt.legend()

    valid = ~np.isnan(all_mags)
    if valid.any():
        q1, q3 = np.percentile(all_mags[valid], [10, 90])
        iqr = q3 - q1
        plt.ylim(q3 + 5 * iqr, q1 - 3 * iqr)
    else:
        plt.ylim(0.7, -0.5)  # fallback


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
    _scatter_lc(X, *mag_arrays, mag_cols=mag_cols)
    _finalize_lc_plot(xlabel, np.concatenate(mag_arrays))

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
    ylim = args.ylim

    data, name_dict, channel_color_dict = read_lc(lc_file, aperture, sep_by, mag_types)
    if args.selected:
        selected = args.selected
        data = {k: v for k, v in data.items() if k in selected}
        if name_dict:
            name_dict = {k: v for k, v in name_dict.items() if k in selected}
        if channel_color_dict:
            channel_color_dict = {k: v for k, v in channel_color_dict.items() if k in selected}

    # Derive mag columns from data (anything that isn't 'bjd')
    sample_df = next(iter(data.values()))
    mag_cols = [c for c in sample_df.columns if c != 'bjd']

    if mode == 'folded':
        xlabel = "Phase" if period else "BJD - 2450000"

        if period:
            if GAIA_plot:
                plot_Gaia_lc(gaia_id, period, epoch, text)

        if sep_by == 'channel' and channel_color_dict:
            # Same frame, but color each data point by its channel.
            all_mags_collected = []
            magnitude_offset = 0.0
            if period and TESS_plot:
                alignment_phases = [
                    calculate_phase(data[chn]['bjd'].values, period, epoch)
                    for chn in sorted(data.keys())
                ]
                alignment_magnitudes = [
                    data[chn][mag_cols[0]].values for chn in sorted(data.keys())
                ]
                magnitude_offset = plot_TESS_lc(
                    gaia_id, period, epoch,
                    np.concatenate(alignment_phases),
                    np.concatenate(alignment_magnitudes),
                    text,
                )
            for chn in sorted(data.keys()):
                df = data[chn]
                bjd_or_phase = (
                    calculate_phase(df['bjd'].values, period, epoch)
                    if period else df['bjd'].values - 2450000
                )
                if binning_method:
                    x_plot, *mag_arrays = apply_binning(
                        df, bjd_or_phase, binning_method, binning_size, mag_cols
                    )
                else:
                    x_plot = bjd_or_phase
                    mag_arrays = [df[col].values for col in mag_cols]
                mag_arrays = [array - magnitude_offset for array in mag_arrays]

                _scatter_lc(
                    x_plot, *mag_arrays, mag_cols=mag_cols,
                    color=channel_color_dict.get(chn, 'gray'),
                    label=name_dict.get(chn, str(chn)) if name_dict else str(chn),
                )
                all_mags_collected.extend(mag_arrays)

            _finalize_lc_plot(xlabel, np.concatenate(all_mags_collected))
            title = f'Gaia {gaia_id} - folded (by channel)'
        else:
            combined_df = pd.concat(
                [data[k] for k in sorted(data.keys())],
                ignore_index=True
            )
            bjd_or_phase = (
                calculate_phase(combined_df['bjd'].values, period, epoch)
                if period else combined_df['bjd'].values - 2450000
            )

            magnitude_offset = 0.0
            if period and TESS_plot:
                magnitude_offset = plot_TESS_lc(
                    gaia_id, period, epoch,
                    bjd_or_phase,
                    combined_df[mag_cols[0]].values,
                    text,
                )

            if binning_method:
                x_plot, *mag_arrays = apply_binning(
                    combined_df, bjd_or_phase, binning_method, binning_size, mag_cols
                )
            else:
                x_plot = bjd_or_phase
                mag_arrays = [combined_df[col].values for col in mag_cols]
            mag_arrays = [array - magnitude_offset for array in mag_arrays]

            plot_lc(x_plot, *mag_arrays, mag_cols=mag_cols, xlabel=xlabel)
            title = f'Gaia {gaia_id} - folded'

        plt.title(title, fontdict={"fontweight":"bold", 'fontsize':12})
        plt.ylim(ylim) if ylim else None
        # plt.xlim(0.83, 0.86) if period else None
        if save_or_show == 'save':
            plt.savefig(
                f'Gaia_{gaia_id}_folded_folded_aperture_{aperture}',
                dpi=400,
                format='pdf'
                )
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
                    magnitude_offset = plot_TESS_lc(
                        gaia_id, period, epoch,
                        bjd_or_phase,
                        df[mag_cols[0]].values,
                        text,
                    )
                else:
                    magnitude_offset = 0.0
                if GAIA_plot:
                    plot_Gaia_lc(gaia_id, period, epoch, text)
            else:
                magnitude_offset = 0.0
            
            if binning_method:
                x_plot, *mag_arrays = apply_binning(
                    df, bjd_or_phase, binning_method, binning_size, mag_cols
                )
            else:
                x_plot = bjd_or_phase
                mag_arrays = [df[col].values for col in mag_cols]
            mag_arrays = [array - magnitude_offset for array in mag_arrays]
            
            plot_lc(x_plot, *mag_arrays, mag_cols=mag_cols, xlabel=xlabel)
            title_part = (
                f"{chn_or_phref}) {name_dict[chn_or_phref]}"
                if name_dict and chn_or_phref in name_dict
                else str(chn_or_phref)
            )
            plt.title(f'{gaia_id} - {sep_by}: {title_part}', fontdict={"fontweight":"bold", 'fontsize':12})
            plt.ylim(ylim) if ylim else None
            if save_or_show == 'save':
                plt.savefig(f'Gaia_{gaia_id}_sepby_{sep_by}_{title_part}_aperture_{aperture}', dpi=400)
            elif save_or_show == 'show':
                plt.show()
            plt.clf()

if __name__ == "__main__":
    main()