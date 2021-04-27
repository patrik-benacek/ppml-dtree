#!/bin/bash
# Interpolation target and auxilliary variables at station locations

LEADTIMES="24 240"
TARGETS="t2m prec24"
PREC_DEACCUM=true

# Interpolation target variables
for LEADTIME in $LEADTIMES; do
    for TARGET in $TARGETS; do
        # Run deaccumulation of precipitation
        if [[ $PREC_DEACCUM = true && "$LEADTIME" = "240" && "$TARGET" = "prec24" ]]; then
            echo "Deaccumulate $TARGET for $LEADTIMES"
            ./run_data_deaccum.sh
        fi
        # Target variables
        src/interp_target_${TARGET}.R $LEADTIME
    done
done

# Auxilliary features 
for LEADTIME in $LEADTIMES; do
    src/interp_aux_srf.R $LEADTIME
    src/interp_aux_plv.R $LEADTIME 850
    src/interp_aux_plv.R $LEADTIME 500
done

# Model orograpy metadata
src/interpolation_aux_geo.R

