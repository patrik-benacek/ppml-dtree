#!/bin/bash                                                                                                                                                                                  

# Grid-search model parameters
RUN_TUNE_MODEL=false
MODEL="ngb"         # qrf/xtr/ngb
TARGET="t2m"        # t2m/prec24
LEADTIME="24h"      # 24h/240h
START_TRAIN="2015"  # start training period

echo "---------------------"
echo "Run experiment:"
echo "---------------------"
# Prepare config file
sed -e "s/xxMODELxx/$MODEL/g" config/config_${TARGET}_${LEADTIME}h_${START_TRAIN}.py > config/config.py
# Set experiment in config.py file
head -n4 config/config.py
                                                                                                                                                                                             
# Activate eccodes environment                                                                                                                                                               
source /home/patrik/miniconda3/bin/activate nwp_pp

if [ $RUN_TUNE_MODEL = true ]; then
    # Hyperparameter tuning (GRID_PARAMS in config)
    echo "Tune model parameters"
    src/tune_model.py
else
    # Train model (PARAMS in config)
    echo "Train model"
    src/train_model.py
fi

# Test model (prediction)
echo "Test model"
src/test_model.py

# Deactivate conda env                                                                                                                                                                       
conda deactivate 
