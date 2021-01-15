#!/bin/bash

DATADIR=/home/patrik/Work/czechglobe/TIGGE/forecasts/data

# Run deaccumulation script
for YYYY in `seq 2015 2019`; do
    for filepath in `ls -C1 ${DATADIR}/${YYYY}/ecmwf_srf_*.nc`; do
        file=`basename ${filepath}` 
        # cut date from input file
        date=`echo ${file} | cut -c 11-18`
        # run deaccumulation
        echo "Deaccumulation file ${file} in ${date}"
        ./field_deaccumulation.R ${DATADIR}/${YYYY}/${file} $DATADIR/${YYYY}/ecmwf_tp_${date}.nc 240 24 $DATADIR/${YYYY}/ecmwf_srf_deacc_${date}.nc
    done
done

echo "Finish"
