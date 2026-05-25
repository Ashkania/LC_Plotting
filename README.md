    Plot the LC from its hdf5 file:
    
    1. Read the LC data from the hdf5 file. It includes:
       bjd, mag_epd, mag_magfit, channels index and photrefs index.
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
