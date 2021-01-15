#!/bin/bash

fcst=24
#vars="t2m prec24"
vars="prec24"
for var in $vars; do
    #for year in 2015 2018; do
    for year in 2018; do
        src/emos_global_fixed_crch.R $var $year $fcst &> listings/emos_crch_global_${var}_f${fcst}_${year}.out 
        src/emos_local_fixed_crch.R $var $year $fcst  &> listings/emos_crch_local_${var}_f${fcst}_${year}.out 
    done
done
