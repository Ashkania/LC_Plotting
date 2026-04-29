#!/usr/bin/env python3

# Given a Gaia id and the path to DR files, it goes over all the
# DR files, extract the x,y, fitted mag (it1) of that Gaia id,
# and report the number of points we find, which means in how many
# dr files this object appears.

# ex of cmdline usage:
# for tel in *; do for cam in $tel/*; do echo "$cam";
# /home/sxj190009/LC_Plotting/scripts/diagnostics/DR_vs_LC.py
# --gaia-id 3334901543952614528 --dr-path "$cam"; done ; done

import numpy as np
import matplotlib.pyplot as plt
import argparse
import h5py
import glob

def parse_args():
    
    parser = argparse.ArgumentParser(
        description="Compare the DR and LC for a given source across multiple files"
        )
    parser.add_argument(
        '--gaia-id',
        help="Gaia ID of the source to analyze"
        )
    parser.add_argument(
        '--dr-path',
        required=True,
        help="Path to DR files to analyze"
        )
    
    return parser.parse_args()

def load_data(dr, gaia_id):
    with h5py.File(dr, 'r') as f:
        try:
            gaia_ids = f['/ProjectedSources/Version000/source_id'][:]
            gaia_index = np.where(gaia_ids == int(gaia_id))[0]
        
            x = f['/ProjectedSources/Version000/x'][:][gaia_index]
            y = f['/ProjectedSources/Version000/y'][:][gaia_index]

            fitted_mag = f['/AperturePhotometry/Version000/Aperture000/FittedMagnitudes/Version000/Iteration001'][:][gaia_index]

        except:
            # print(f"Gaia ID {gaia_id} not found in DR file")
            return None, None, None
        
    return x, y, fitted_mag




def main():
    args = parse_args()
    X = []
    Y = []
    Fitted_Mag = []
    for dr in glob.glob(args.dr_path + "/*.h5"):
        # print(f"Processing DR file: {dr}")
        x, y, fitted_mag = load_data(dr, args.gaia_id)
        if x is not None and x.size > 0:
            X.append(x.item())
            Y.append(y.item())
            if fitted_mag is not None and fitted_mag.size > 0 and fitted_mag > -1000:
                Fitted_Mag.append(fitted_mag.item())
    if X and Y:
        print(f'X) min: {min(X)}, max: {max(X)}, len: {len(X)}')
        print(f'Y) min: {min(Y)}, max: {max(Y)}, len: {len(Y)}')
        print(f'FitMag) min: {min(Fitted_Mag)}, max: {max(Fitted_Mag)}, len: {len(Fitted_Mag)}')
    else:
        print("No valid data found for the given Gaia ID across the DR files.")

if __name__ == "__main__":
    main()
    