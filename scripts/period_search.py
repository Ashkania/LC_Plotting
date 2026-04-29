#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
from plotting import read_lc
import pandas as pd
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

    return parser.parse_args()


def main():
    args = parse_arguments()
    lc = args.lc_file
    # data = {0: df0, 1: df1, ...} keys are either channels or photref
    data, photref_dict = read_lc(lc, aperture=4, sep_by='photref')

    all_bjds_list_of_lists = []
    # all_mag_epd_list_of_lists = []
    all_mag_magfit_list_of_lists = []

    for chn_or_phref in sorted(data.keys()):
        df = data[chn_or_phref]
        all_bjds_list_of_lists.append(df['bjd'].values)  # combination of separate arrays
        # all_mag_epd_list_of_lists.append(df['mag_epd'].values)
        all_mag_magfit_list_of_lists.append(df['mag_magfit'].values)

        all_bjds = np.concatenate(all_bjds_list_of_lists)  # make one big array of bjd values
        # all_mag_epd = np.concatenate(all_mag_epd_list_of_lists)
        all_mag_magfit = np.concatenate(all_mag_magfit_list_of_lists)

        combined_df = pd.DataFrame({
            'bjd': all_bjds,
            # 'mag_epd': all_mag_epd,
            'mag_magfit': all_mag_magfit
        })

    times = combined_df['bjd']
    mags = combined_df['mag_magfit']

    new_mask = (
        np.isfinite(times)
        &  np.isfinite(mags)
    )

    times = times[new_mask]
    mags = mags[new_mask]

    ls = LombScargle(
        times,
        mags
        )

    frequency, power = ls.autopower(
        # minimum_frequency=0.05,
        # maximum_frequency=100
    )

    power_sorted_ix = np.argsort(power)
    best_f = frequency[power_sorted_ix[-1]]
    best_f2 = frequency[power_sorted_ix[-2]]
    best_p = 1/best_f
    best_p2 = 1/best_f2


    # print("Best frequency:", best_f)
    # print("Best period:", best_p)

    # print("2nd Best frequency:", best_f2)
    # print("2nd Best period:", best_p2)

    print(best_p)
    print(best_p2)

if __name__ == '__main__':
    main()