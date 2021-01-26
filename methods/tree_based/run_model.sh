#!/bin/bash

# Grid-search model parameters
MODEL="qrf"          # qrf/xtr/ngb
TARGET="t2m"         # t2m/prec24
LEADTIME="24h"       # 24h/240h
START_TRAIN="2018"   # start training period
RUN_TUNE_MODEL=false # start model tuning

echo "---------------------"
echo "Run experiment:"
echo "---------------------"
# Prepare config file
sed -e "s/xxMODELxx/$MODEL/g" config/config_${TARGET}_${LEADTIME}_${START_TRAIN}.py > config/config.py
# Set experiment in config.py file
head -n4 config/config.py

# Activate eccodes environment
source /home/patrik/miniconda3/bin/activate nwp_pp

if [ $RUN_TUNE_MODEL = true ]; then
    # Hyperparameter tuning (GRID_PARAMS in config)
    echo "Tune model parameters"
    (time python src/run.py tune) &> models/summary_tune_${TARGET}_${LEADTIME}_${START_TRAIN}.out 
else
    # Train model (PARAMS in config)
    echo "Train model"
    (time python src/run.py train) &> models/summary_train_${TARGET}_${LEADTIME}_${START_TRAIN}.out
fi

# Test model (prediction)
echo "Test model"
(time python src/run.py test) &> models/summary_test_${MODEL}_${TARGET}_${LEADTIME}_${START_TRAIN}.out 

# Deactivate conda env
conda deactivate 
