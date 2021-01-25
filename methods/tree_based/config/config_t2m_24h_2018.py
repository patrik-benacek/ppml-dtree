TARGET   = "t2m"        # t2m/prec24
LEADTIME = "24h"        # 24h/240h
MODEL    = "xxMODELxx"  # qrf/xtr/ngb
STATION  = "all"        # all/station_name

DATA_DIR = "../../../data_preproc/data/processed/"
MODEL_DIR = "../models/"

TRAIN_PERIOD = [2018]
TEST_PERIOD  = [2019]

# Model parameter grid-search setting
GRIDSEARCH_CV = 3
N_QUANTILES_PREDICT = 50
NUM_CORES__GRID_SEARCH = 2
NUM_CORES__MODEL = 5

# Grid-search model parameter setting (tune_model.py):
GRID_PARAMS = {
    'ngb': {
        'model__Base__max_depth': [3, 5],
        'model__n_estimators': [200],
        'model__minibatch_frac': [0.2],
        'model__learning_rate': [0.03]
    },
    'qrf': {
        'model__n_estimators': [500],
        'model__max_features': [25],
        'model__min_samples_leaf': [2, 10],
        'model__bootstrap': [False]
    },
    'xtr': {
        'model__n_estimators': [500],
        'model__max_features': [25],
        'model__min_samples_leaf': [2, 10],
        'model__bootstrap': [False]
    }
}

# Model parameters setting (train_model.py)
PARAMS = {
    'ngb': {
        'model__Base__max_depth': 4,
        'model__n_estimators': 300,
        'model__minibatch_frac': 0.5,
        'model__learning_rate': 0.05
    },
    'qrf': {
        'model__n_estimators': 500,
        'model__max_features': 25,
        'model__min_samples_leaf': 2,
        'model__bootstrap': False
    },
    'xtr': {
        'model__n_estimators': 500,
        'model__max_features': 50,
        'model__min_samples_leaf': 2,
        'model__bootstrap': False
    }
}
