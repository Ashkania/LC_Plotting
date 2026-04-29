#!/bin/bash


for file_path in /home/ashkan/ASHKAN/Projects/LC_Plotting/Inputs/LC/GAIA_variables/GDR3_*_RS.h5; do
    mapfile -t periods < <(./period_search.py "$file_path" 2>/dev/null | grep -E '^[0-9.]+$')
    for p_ix in 0 1; do
        echo "${periods[$p_ix]}"
        python ~/ASHKAN/Projects/LC_Plotting/scripts/plotting.py "$file_path"\
        --aperture 3\
        --sep-by 'photref'\
        --mode 'folded'\
        --period "${periods[$p_ix]}"\
        --epoch 1705\
        --binning time 1\
        --text period_${p_ix}_${periods[$p_ix]} \
        --GAIA-plot\
        --save-or-show-plot save
        mv ./*png ./p${p_ix}
    done
done
