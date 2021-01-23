# Model skgarden-qrf is too slow
# Use pyquantr instead!!

TARGET = "t2m"      # t2m/prec24
LEADTIME = "24h"    # 24h/240h
MODEL = "quantrf"   # qrf/xtr/ngb/quantrf
STATION="Cheb"      # all/station_name

DATA_DIR = "../../../data_preproc/data/processed/"
MODEL_DIR = "../models/"

TRAIN_PERIOD = [2015,2016,2017,2018]
TEST_PERIOD  = [2019]

# Model parameter grid-search setting
GRIDSEARCH_CV = 3
N_QUANTILES_PREDICT = 50
NUM_CORES = 5

# Set grid_params values as list()!!
GRID_PARAMS = {
    'model__n_estimators': [500],
    'model__max_features': [25],
    'model__min_samples_leaf': [2, 10],
    'model__bootstrap': [False]
}

PARAMS = {
    'model__n_estimators': 500,
    'model__max_features': 25,
    'model__min_samples_leaf': 2,
    'model__bootstrap': False
}
