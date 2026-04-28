#!/usr/bin/env bash

# Script to plot eclipsing binaries from villa_eclipsing_binaries.txt

# Check if the file exists

##############
# 1. prsa_eclipsing_binaries.txt  villa_eclipsing_binaries.txt  NEA_TOI.txt   NEA_exoplanets.txt
# 2. EB_10k  EB_Villa   TOI    Exo
file_name="NEA_exoplanets.txt"
lc_path="Exo"
##############

if [[ ! -f "../Results/objects_lists/${file_name}" ]]; then
    echo "Error: ${file_name}  not found"
    exit 1
fi

# Skip the header line and process each data line
tail -n +2 ../Results/objects_lists/${file_name} | while IFS=',' read -r gaia_id tic period tmag epoch note; do
    # Trim whitespace from each variable
    tmag=$(echo "$tmag" | xargs)
    gaia_id=$(echo "$gaia_id" | xargs)
    period=$(echo "$period" | xargs)
    epoch=$(echo "$epoch" | xargs)

    # Construct the LC filename
    lc_file="../Inputs/LC/${lc_path}/EPD_all_channels/GDR3_${gaia_id}.h5"

    # Check if the LC file exists
    if [[ ! -f "$lc_file" ]]; then
        echo "Warning: LC file not found: $lc_file - skipping"
        continue
    fi

    echo "Processing: $lc_file (Period: $period, Epoch: $epoch)"

    # Run the plotting script
    for ap in {0..7}; do
        python3 plotting.py "$lc_file" \
            --aperture "$ap" \
            --period "$period" \
            --epoch "$epoch" \
            --sep-by 'photref' \
            --mode 'folded' \
            --binning 'time' 4\
            --TESS-plot \
            --tmag "$tmag" \
            --save-or-show-plot "save"
    done
done

echo "Done processing all objects"
