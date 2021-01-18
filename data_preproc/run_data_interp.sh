#!/bin/bash
# Interpolation target and auxilliary variables at station locations

LEADTIMES="24 240"
TARGETS="t2m prec24"

# Interpolation target variables
for LEADTIME in $LEADTIMES; do
    for TARGET in $TARGETS; do
        # Target variables
        src/interpolation_target_${TARGET}.R $LEADTIME
    done
done

# Auxilliary features 
for LEADTIME in $LEADTIMES; do
    src/interpolation_aux_srf.R $LEADTIME
    src/interpolation_aux_plv.R $LEADTIME 850
    src/interpolation_aux_plv.R $LEADTIME 500
done

# Model orograpy metadata
src/interpolation_aux_geo.R

