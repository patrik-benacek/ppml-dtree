import numpy as np
import pandas as pd 
import os

import joblib
import itertools
from multiprocessing.pool import ThreadPool as Pool
from functools import partial

from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from transformers import GetOrographyError
from skgarden import RandomForestQuantileRegressor, ExtraTreesQuantileRegressor
from skgarden import ExtraTreesQuantileRegressor
from ngboost import NGBRegressor
from quant_forest import QuantileRandomForestRegressor
import properscoring as ps

from make_dataset import read_dataset

from config import (TRAIN_PERIOD, TEST_PERIOD, 
                    MODEL, MODEL_DIR, PARAMS, GRID_PARAMS,
                    N_QUANTILES_PREDICT, GRIDSEARCH_CV, NUM_CORES)

def split_train_test(df):
    """Split dataset to train/test based on the predefined $train_period period."""

    # Get target and features
    x = df.drop(columns=['obs'])
    y = df['obs']

    # Split dataset
    x_train, y_train = x[str(TRAIN_PERIOD[0]):str(TRAIN_PERIOD[-1])], y[str(TRAIN_PERIOD[0]):str(TRAIN_PERIOD[-1])]
    x_test, y_test   = x[str(TEST_PERIOD[0]):str(TEST_PERIOD[-1])], y[str(TEST_PERIOD[0]):str(TEST_PERIOD[-1])]

    # Set training period
    print("Split dataset...")
    print("Train: {}-{}".format(x_train.index.min().year, x_train.index.max().year))
    print("Test:  {}-{}".format(x_test.index.min().year, x_test.index.max().year))
    print()
    return x_train, y_train, x_test, y_test

def build_model():
    """Building model""" 
    # Data preprocessor
    preprocessor = ColumnTransformer([('orog', GetOrographyError(), ["alt", "orog"])], remainder='passthrough')

    # Select Model
    if MODEL=="qrf":
        get_model = RandomForestQuantileRegressor(n_jobs=NUM_CORES, random_state=42)
    elif MODEL=="quantrf":
        get_model = QuantileRandomForestRegressor(nthreads=NUM_CORES, random_state=42)
    elif MODEL=="xtr":
        get_model = ExtraTreesQuantileRegressor(random_state=42)
    elif MODEL=="ngb":
        get_model = NGBRegressor()

    # Model pipeline
    model = Pipeline([
        ('preproc', preprocessor), 
        ('imputer', SimpleImputer(strategy='mean')),
        ('model', get_model)
    ])
    return model 

def crps_ensemble_score(test, pred):
    return ps.crps_ensemble(test, pred).mean()

def crps_gaussian_score(test, loc, scale):
    return ps.crps_gaussian(test, mu=loc, sig=scale).mean()

# Scikit-garden solution of hyperparameter tuning:
def predict_approx(model, X_test, quantiles=[0.05, 0.5, 0.95]):
    """
    Function to predict quantiles much faster than the default skgarden method
    This is the same method that the ranger and quantRegForest packages in R use
    Output is (n_samples, n_quantiles) or (n_samples, ) if a scalar is given as quantiles:
    https://stackoverflow.com/questions/51483951/quantile-random-forests-from-scikit-garden-very-slow-at-making-predictions
    """
    # Begin one-time calculation of random_values. This only depends on model, so could be saved.
    n_leaves = np.max(model.y_train_leaves_) + 1  # leaves run from 0 to max(leaf_number)
    random_values = np.zeros((model.n_estimators, n_leaves))
    for tree in range(model.n_estimators):
        for leaf in range(n_leaves):
            train_samples = np.argwhere(model.y_train_leaves_[tree, :] == leaf).reshape(-1)
            if len(train_samples) == 0:
                random_values[tree, leaf] = np.nan
            else:
                train_values = model.y_train_[train_samples]
                random_values[tree, leaf] = np.random.choice(train_values)
    # Optionally, save random_values as a model attribute for reuse later

    # For each sample, get the random leaf values from all the leaves they land in
    X_leaves = model.apply(X_test)
    leaf_values = np.zeros((X_test.shape[0], model.n_estimators))
    for i in range(model.n_estimators):
        leaf_values[:, i] = random_values[i, X_leaves[:, i]]

    # For each sample, calculate the quantiles of the leaf_values
    return np.quantile(leaf_values, np.array(quantiles), axis=1).transpose()

def run_estimator_(hpar, data, use_predict_approx=True):
    #from sklearn.base import clone
    # import warnings filter
    from warnings import simplefilter
    # ignore all future warnings
    simplefilter(action='ignore', category=FutureWarning)

    X_train, y_train, X_test, y_test = split_train_test(data)
    # Build model
    model = build_model()
    # Input grid-search values to dicitionary 
    dic_params = {param:value for param, value in zip(GRID_PARAMS.keys(), hpar)}
    # Set model parameters
    model.set_params(**dic_params)
    # estimator 
    model.fit(X_train, y_train)
    # Get prediction for particular quantiles
    nq = N_QUANTILES_PREDICT + 1
    if use_predict_approx:
        quantiles = np.arange(1/nq, nq/nq, 1/nq)
        X_test_ = model[:-1].fit_transform(X_test)
        y_pred = predict_approx(model[-1], X_test_, quantiles=quantiles)
    else:
        quantiles = 100 * np.arange(1/nq, nq/nq, 1/nq)
        y_pred = np.vstack([model.predict(X_test, quantile=q) for q in quantiles]).T
    # Evaluation
    score = crps_ensemble_score(y_test, y_pred)
    return score

def tune_model_():
    """Tune model parameters for skgarden methods: 
        * RandomForestQuantile
        * ExtraTreesQuantile
    These methods are time-demanding due to n_jobs issues (only n_jobs=1 is supported).
    Therefore, we apply multiprocessing for each grid parameter. Num of cores is defined
    by NUM_CORES in config file.
    """

    # Test for skgarden methods
    if not MODEL == "qrf" or not MODEL == "xtr":
        sys.exit("Function tune_model_ is supported only for skgarden methods. Change MODEL in the config file.")

    # Training data
    data, _ = read_dataset()
    #X_train, y_train, X_test, y_test = split_train_test(data)
    # Build model
    #model = build_model()

    # pairs for processors
    pool = Pool(processes=NUM_CORES)
    grid_pairs = itertools.product(*GRID_PARAMS.values())
    scores = pool.map_async(partial(run_estimator_, data=data), grid_pairs).get()
    #run workers
    pool.close()
    pool.join()
    # Find hyperparameters for minimal scores
    # TODO
    return scores

def run_model_grid_search(model, X_train, y_train, X_test, y_test):
    """Tune model parameters for non scikit-learn models.
    We need to separate model from transformers, otherwise error. """

    # Set model parameters
    nq = N_QUANTILES_PREDICT + 1
    quantiles = np.arange(1/nq, nq/nq, 1/nq)
    grid_pairs = itertools.product(*GRID_PARAMS.values())
    scores = list()
    for hpar in grid_pairs:
        # Input grid-search values to dicitionary 
        dic_params = {param:value for param, value in zip(GRID_PARAMS.keys(), hpar)}
        # Set model parameters
        model.set_params(**dic_params)
        # Prepare train/test datasets
        X_train_ = model[:-1].fit_transform(X_train)
        X_test_  = model[:-1].fit_transform(X_test)
        y_train_ = y_train.to_numpy()
        # Fit model
        model[-1].fit(X_train_, y_train_)
        # Predict model
        y_pred = model[-1].predict(X_test_, quantiles)
        # Evaluate model 
        scores.append(crps_ensemble_score(y_test, y_pred))
    return scores

def tune_model():
    """Grid search model parameters using the cross-validation method."""

    # Load dataset
    data, expname = read_dataset()
    # Train/test datasets
    X_train, y_train, X_test, y_test = split_train_test(data)
    # Build model
    model = build_model()

    # Grid parameters initialization
    grid_pairs = list(itertools.product(*GRID_PARAMS.values()))
    n_grid_pairs = len(grid_pairs)
    gs_scores = np.zeros([n_grid_pairs])

    # Cross-validation folds 
    nfolds  = GRIDSEARCH_CV 
    dataIDX = np.arange(X_train.shape[0])
    folds   = np.array_split(dataIDX, nfolds)

    for f in range(nfolds):
        # Set CV folds
        out_fold = folds[f]
        in_folds = np.concatenate([folds[fold] for fold in range(nfolds) if fold != f])
        # Run grid-search model parameters 
        X_train_ = X_train.iloc[in_folds]
        y_train_ = y_train.iloc[in_folds]
        X_val    = X_train.iloc[out_fold]
        y_val    = y_train.iloc[out_fold]
        gs_score = run_model_grid_search(model, X_train_, y_train_, X_val, y_val)
        gs_scores += gs_score
    # Mean CV-grid-search score 
    gs_scores = gs_scores/nfolds

    # Best model parameters
    best_grid_pair = grid_pairs[np.argmin(gs_score)]
    best_model_params = {param:value for param, value in zip(GRID_PARAMS.keys(), best_grid_pair)}
    print("Best GridSearchCV model parameters:\n{}".format(best_model_params))
    print("GridSearchCV score: {:.2f}".format(gs_score[np.argmin(gs_score)]))
    print()
    # Tune best model
    model.set_params(**best_model_params)
    # train/test 
    X_train_ = model[:-1].fit_transform(X_train)
    y_train_ = y_train.to_numpy()
    model[-1].fit(X_train_, y_train_)
    # Save model
    file_model = os.path.join(MODEL_DIR, f"model_{expname}.joblib")
    print("Save model to: {}".format(file_model))
    joblib.dump(model, file_model)

def train_model(print_params=True):
    # Read dataset
    data, expname = read_dataset()
    # Get training data
    X_train, y_train, _, _ = split_train_test(data)
    # Build model
    model = build_model()
    model.set_params(**PARAMS)

    if print_params:
        print(model.get_params())

    model.fit(X_train, y_train)
    # Save model
    joblib.dump(model, os.path.join(MODEL_DIR, f"model_{expname}.joblib"))

def test_model():
    # Read dataset
    data, expname = read_dataset()
    # Get testing data
    _, _, X_test, y_test = split_train_test(data)
    # Get trained model
    model = joblib.load(os.path.join(MODEL_DIR, f"model_{expname}.joblib"))

    # Make prediction
    nq = N_QUANTILES_PREDICT + 1
    quantiles = np.arange(1/nq, nq/nq, 1/nq)
    # train/test 
    X_test_  = model[:-1].fit_transform(X_test)
    y_pred = model[-1].predict(X_test_, quantiles) 

    # Make prediction for particular quantiles
    #y_pred = predict_approx(model, X_test, quantiles=np.arange(1/51, 51/51, 1/51))
    #how incorrporate to predict_approx: model.named_steps['model']

    print("CRPS on the test set: {:.2f}".format(
        crps_ensemble_score(y_test, y_pred)))