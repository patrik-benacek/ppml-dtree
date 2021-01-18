#!/bin/bash

# Merge dataset with station and model metadata
# Output:
#   * data/merged/data_$target_$leadtime.csv

LEADTIMES="24 240"
TARGETS="t2m prec24"

for LEADTIME in $LEADTIMES; do
    for TARGET in $TARGETS; do
        # Merge NetCDF forecast data: data_$target_$leadtime.csv
        src/merge_interp_data.R $LEADTIME $TARGET

        # Add metadata: data_wmeta_$target_$leadtime.csv
        src/merge_meta_data.R $LEADTIME $TARGET

        # Overwrite the dataset with the enhanced one
        mv data/merged/data_wmeta_${TARGET}_${LEADTIME}.csv data/merged/data_${TARGET}_${LEADTIME}.csv
    done
done
