#!/usr/bin/env python3

# Input: list of interested Gaia Ids, and a single DR file
# Output: for each Gaia Id, find the x,y in the DR file and
# generate the region file for ds9.

import numpy as np
import argparse
from h5py import File

def parse_args():
    
    parser = argparse.ArgumentParser(
        description="Generate DS9 region files for given Gaia IDs from DR files"
        )
    parser.add_argument(
        '--gaia-ids-file',
        help="File containing Gaia IDs of the sources to analyze"
        )
    parser.add_argument(
        '--dr-file',
        required=True,
        help="Path to DR file to analyze"
        )
    parser.add_argument(
        '--output-reg-file',
        default='reg.reg',
        help="Output DS9 region file name"
    )
    
    return parser.parse_args()


def main():
    args = parse_args()
    regfile = args.output_reg_file
    drfile = args.dr_file
    gaia_ids = np.loadtxt(args.gaia_ids_file, dtype=str)
    tot=0

    with open(regfile, 'w') as f_reg:
        with File(drfile, 'r') as f_dr:
            gaia_ids_arr = f_dr['/ProjectedSources/Version000/source_id'][:]
            # print(f"Total Gaia IDs in DR file: {gaia_ids_arr[:5]}")
            for gaia_id in gaia_ids:
                # print(f"Processing Gaia ID: {gaia_id}")
                try:
                    gaia_index = np.where(gaia_ids_arr == int(gaia_id))[0]
                    # print(f"Found Gaia ID {gaia_id} at index: {gaia_index}")
                    x = f_dr['/ProjectedSources/Version000/x'][:][gaia_index][0]
                    y = f_dr['/ProjectedSources/Version000/y'][:][gaia_index][0]
                    # print(x,y)
                    print(
                        'box({xc!r},{yc!r},{w!r},{h!r},0) # color=green'.format(
                            xc=x + 0.5,
                            yc=y + 0.5,
                            w=15,
                            h=15
                        ), file = f_reg
                    )
                    print(
                        'text({xc},{yc}) # text="{text}", color=green'.format(
                            xc=x + 0.5,
                            yc=y - 8.0,
                            text=gaia_id
                        ),
                        file=f_reg
                    )
                    tot+=1
                except:
                    print(f"Gaia ID {gaia_id} not found in DR file")
                    continue
    print(f'Total {tot} sources found and added to region file {regfile}')
if __name__ == "__main__":
    main()
    