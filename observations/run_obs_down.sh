#!/bin/bash
#
# Downloading station observations from OGIMET:
# https://www.ogimet.com/synopsc.phtml.en 

# Data description:
# https://rdrr.io/cran/climate/man/meteo_ogimet.html

OBS_COUNTRY="Czech+Republic"
YEARS="2015 2016 2017 2018 2019"

for YEAR in $YEARS; do
    echo "Download obs for $YEAR"

    # Downloading data
    src/obs_download.R $OBS_COUNTRY $YEAR

    # Get sum of 24-hour precipitation
    src/obs_sum_prec.R $YEAR
done
