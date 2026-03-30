#!/usr/bin/env python3

# TODO: Given a Gaia id, go over all the DR files, extract the x,y
# of that source, and compare the number of points we find with 
# what we have in the LC

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
            # print(gaia_ids)
            # print(gaia_id)
            gaia_index = np.where(gaia_ids == int(gaia_id))[0]
            # print(gaia_index)
        
            x = f['/ProjectedSources/Version000/x'][:][gaia_index]
            y = f['/ProjectedSources/Version000/y'][:][gaia_index]
        except:
            # print(f"Gaia ID {gaia_id} not found in DR file")
            x, y = None, None
    
    return x.item(), y.item()




def main():
    args = parse_args()
    X = []
    Y = []
    for dr in glob.glob(args.dr_path + "/*.h5"):
        # print(f"Processing DR file: {dr}")
        x, y = load_data(dr, args.gaia_id)
        if x and y:
            X.append(x)
            Y.append(y)
    
    print(f'X) min: {min(X)}, max: {max(X)}, len: {len(X)}')
    print(f'Y) min: {min(Y)}, max: {max(Y)}, len: {len(Y)}')

if __name__ == "__main__":
    main()
    