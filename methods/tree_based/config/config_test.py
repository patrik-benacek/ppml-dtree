TARGET   = "prec24"    # t2m/prec24
LEADTIME = "24h"    # 24h/240h
MODEL    = "qrf"    # qrf/xtr/ngb
STATION  = "Cheb"   # all/station_name

DATA_DIR = "../../../data_preproc/data/processed/"
MODEL_DIR = "../models/"

TRAIN_PERIOD = [2015, 2016, 2017, 2018]
TEST_PERIOD  = [2019]

# Model parameter grid-search setting
GRIDSEARCH_CV = 3
N_QUANTILES_PREDICT = 50
NUM_CORES__GRID_SEARCH = 3
NUM_CORES__MODEL = 6

# Grid-search model parameter setting (tune_model.py):
GRID_PARAMS = {
    'ngb': {
        'model__Base__max_depth': [3],
        'model__n_estimators': [300, 500, 600],
        'model__minibatch_frac': [0.2],
        'model__learning_rate': [0.01]
    },
    'qrf': {
        'model__n_estimators': [500],
        'model__max_features': [25],
        'model__min_samples_leaf': [10],
        'model__bootstrap': [False]
    },
    'xtr': {
        'model__n_estimators': [500],
        'model__max_features': [25],
        'model__min_samples_leaf': [5],
        'model__bootstrap': [False]
    }
}

# Model parameters setting (train_model.py)
PARAMS = {
    'ngb': {
        'model__Base__max_depth': 3,
        'model__n_estimators': 200,
        'model__minibatch_frac': 0.2,
        'model__learning_rate': 0.03
    },
    'qrf': {
        'model__n_estimators': 500,
        'model__max_features': 25,
        'model__min_samples_leaf': 10,
        'model__bootstrap': False
    },
    'xtr': {
        'model__n_estimators': 500,
        'model__max_features': 25,
        'model__min_samples_leaf': 5,
        'model__bootstrap': False
    }
}
