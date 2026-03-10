#!/usr/bin/env bash

# Script to plot eclipsing binaries from villa_eclipsing_binaries.txt

# Check if the file exists
# if [[ ! -f "villa_eclipsing_binaries.txt" ]]; then
# if [[ ! -f "prsa_eclipsing_binaries.txt" ]]; then
# if [[ ! -f "toi_in_fov.txt" ]]; then
if [[ ! -f "exoplanets_in_fov.txt" ]]; then
    echo "Error: exoplanets_in_fov.txt not found"
    exit 1
fi

# Skip the header line and process each data line
tail -n +2 exoplanets_in_fov.txt | while IFS=',' read -r gaia_id tic period tmag epoch; do
    # Trim whitespace from each variable
    gaia_id=$(echo "$gaia_id" | xargs)
    period=$(echo "$period" | xargs)
    epoch=$(echo "$epoch" | xargs)
    
    # Construct the LC filename
    # lc_file="LC/EB_Villa/GDR3_${gaia_id}.h5"
    # lc_file="LC/EB_10k/GDR3_${gaia_id}.h5"
    # lc_file="LC/TOI/GDR3_${gaia_id}.h5"
    lc_file="LC/Exo/GDR3_${gaia_id}.h5"
    
    # Check if the LC file exists
    if [[ ! -f "$lc_file" ]]; then
        echo "Warning: LC file not found: $lc_file - skipping"
        continue
    fi
    
    echo "Processing: $lc_file (Period: $period, Epoch: $epoch)"
    
    # Run the plotting script
    python plotting.py "$lc_file" \
        --aperture 4 \
        --period "$period" \
        --epoch "$epoch" \
        --fold-by 'photref' \
        --plot-type 'combined' \
        --binning 'time' 10
    
done

echo "Done processing all objects"
