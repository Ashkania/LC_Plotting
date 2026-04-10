#!/usr/bin/env python3

import lightkurve as lk
import matplotlib.pyplot as plt
import numpy as np
from argparse import ArgumentParser

def parse_arguments():
    parser = ArgumentParser(description="Plot TESS Light Curve")
    parser.add_argument(
        "--gaia-id",
        type=str,
        default="Gaia DR3 3337025761759515776",
        help="Gaia source ID (default: Gaia DR3 3337025761759515776)"
    )
    parser.add_argument(
        "--period",
        type=float,
        default=1.7801492,
        help="Period in days (default: 1.7801492)"
    )
    parser.add_argument(
        "--epoch",
        type=float,
        default=1469.56553,
        help="Reference epoch (t0) for phase folding (default: 1469.56553)"
    )
    parser.add_argument(
        "--mode",
        choices=['original', 'folded', 'both'],
        default='both',
        help="Plot mode: 'original' for time series, 'folded' for phase-folded, 'both' for both (default: both)"
    )
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    search = lk.search_lightcurve(args.gaia_id, mission="TESS")
    print(search)
    
    lc = search.download()
    lc = lc.remove_nans().normalize()
    lc = lc.flatten(window_length=401)
    
    time_tess = lc.time.value
    flux_tess = lc.flux.value
    mag_tess = -2.5 * np.log10(flux_tess)
    
    if args.mode in ['folded', 'both']:
        phase = ((time_tess - args.epoch) / args.period) % 1
        phase = (phase + 0.25) % 1
        sort_idx = np.argsort(phase)
        phase = phase[sort_idx]
        flux_fold = flux_tess[sort_idx]
        mag_fold = mag_tess[sort_idx]
    
    if args.mode == 'original':
        fig, ax = plt.subplots(1, 1, figsize=(12, 5))
        ax.scatter(time_tess, flux_tess, s=1, label="TESS", alpha=0.7)
        ax.set_xlabel("Time (BTJD)")
        ax.set_ylabel("Normalized Flux")
        ax.set_title(f"TESS LC: {args.gaia_id}")
        ax.legend()
        ax.grid(True)
        
    elif args.mode == 'folded':
        fig, ax = plt.subplots(2, 1, figsize=(10, 8))
        ax[0].scatter(phase, flux_fold, s=1, label="TESS", alpha=0.7)
        ax[0].set_xlabel("Phase")
        ax[0].set_ylabel("Normalized Flux")
        ax[0].set_title(f"Phase-folded TESS LC (P = {args.period} d)")
        ax[0].legend()
        ax[0].grid(True)
        
        ax[1].scatter(phase, mag_fold, s=1, label="TESS", alpha=0.7)
        ax[1].set_xlabel("Phase")
        ax[1].set_ylabel("Relative Magnitude")
        ax[1].set_title(f"TESS LC: {args.gaia_id}, P = {args.period} d")
        ax[1].invert_yaxis()
        ax[1].legend()
        ax[1].grid(True)
        
    elif args.mode == 'both':
        fig, ax = plt.subplots(2, 2, figsize=(14, 8))
        
        ax[0, 0].scatter(time_tess, flux_tess, s=1, label="TESS", alpha=0.7)
        ax[0, 0].set_xlabel("Time (BTJD)")
        ax[0, 0].set_ylabel("Normalized Flux")
        ax[0, 0].set_title(f"TESS LC: {args.gaia_id}")
        ax[0, 0].legend()
        ax[0, 0].grid(True)
        
        ax[0, 1].scatter(time_tess, mag_tess, s=1, label="TESS", alpha=0.7)
        ax[0, 1].set_xlabel("Time (BTJD)")
        ax[0, 1].set_ylabel("Relative Magnitude")
        ax[0, 1].set_title(f"TESS LC: {args.gaia_id} (Magnitude)")
        ax[0, 1].invert_yaxis()
        ax[0, 1].legend()
        ax[0, 1].grid(True)
        
        ax[1, 0].scatter(phase, flux_fold, s=1, label="TESS", alpha=0.7)
        ax[1, 0].set_xlabel("Phase")
        ax[1, 0].set_ylabel("Normalized Flux")
        ax[1, 0].set_title(f"Phase-folded TESS LC (P = {args.period} d)")
        ax[1, 0].legend()
        ax[1, 0].grid(True)
        
        ax[1, 1].scatter(phase, mag_fold, s=1, label="TESS", alpha=0.7)
        ax[1, 1].set_xlabel("Phase")
        ax[1, 1].set_ylabel("Relative Magnitude")
        ax[1, 1].set_title(f"Phase-folded TESS LC (Magnitude, P = {args.period} d)")
        ax[1, 1].invert_yaxis()
        ax[1, 1].legend()
        ax[1, 1].grid(True)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()