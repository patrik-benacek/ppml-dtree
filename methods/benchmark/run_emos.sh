#!/bin/bash

# Run global- and local-EMOS methods for $TARGETS and $LEADTIMES. 
# Training period is from $START_TRAIN_YEAR to 2018.
# Testing period is 2019.

LEADTIMES="24 240"
TARGETS="t2m prec24"
START_TRAIN_YEARS="2015 2018"

for LEATIME in $LEADTIMES; do
    for TARGET in $TARGETS; do
        for START_TRAIN_YEAR in $START_TRAIN_YEARS; do
            src/emos_global_fixed_crch.R $TARGET $LEADTIME $START_TRAIN_YEAR &> log_emos_glb_${var}_f${fcst}_${year}.out 
            src/emos_local_fixed_crch.R $TARGET $LEADTIME $START_TRAIN_YEAR  &> log_emos_loc_${var}_f${fcst}_${year}.out 
        done
    done
done
