#!/bin/bash

# Merge dataset with station and model metadata
# Output:
#   * data/merged/data_$target_$leadtime.csv

#LEADTIMES="24 240"
#TARGETS="t2m prec24"
LEADTIMES="24"
TARGETS="prec24"

for LEADTIME in $LEADTIMES; do
    for TARGET in $TARGETS; do
        echo "Merge interpolated data to data-frame ..."
        # Merge NetCDF forecast data: data_$target_$leadtime.csv
        src/merge_interp_data.R $LEADTIME $TARGET

        echo "Add metadata ..."
        # Add metadata: data_wmeta_$target_$leadtime.csv
        src/merge_meta_data.R $LEADTIME $TARGET

        echo "Create tarball..."
        # Overwrite the dataset with the enhanced one
        mv data/processed/data_wmeta_${TARGET}_ff${LEADTIME}h.csv data/processed/data_${TARGET}_ff${LEADTIME}h.csv

        # Zip processed files and clean
        zip data/processed/data_${TARGET}_ff${LEADTIME}h data/processed/data_${TARGET}_ff${LEADTIME}h.csv
        rm -f data/processed/data_${TARGET}_ff${LEADTIME}h.csv
    done
done
