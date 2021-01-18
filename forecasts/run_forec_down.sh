#!/bin/bash

# Downloading the foreast data
OUTDIR="data/"
# Orography 
src/down_ecmwf_forecasts_geo.py $OUTDIR
# Pressure level 
src/down_ecmwf_forecasts_plv.py $OUTDIR
# Surface level
src/down_ecmwf_forecasts_srf.py $OUTDIR
# Data for deaccumulation of ff240h precipitation 
src/down_ecmwf_forecasts_tp_deaccum.py $OUTDIR
