import numpy as np
import pandas as pd 
import os
import sys
import joblib
import itertools
from functools import partial
from multiprocessing.pool import ThreadPool as Pool

from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.tree import DecisionTreeRegressor
from ngboost import NGBRegressor
from ngboost.scores import CRPScore
from ngboost.distns import Normal
import properscoring as ps
import scipy

sys.path.append('config')
from core_trees import QuantileRandomForestRegressor, QuantileExtraTreesRegressor
from transformers import GetOrographyError
from utils import get_feature_names
from make_dataset import read_dataset
from config import (TRAIN_PERIOD, TEST_PERIOD, TARGET, STATION, LEADTIME,
                    MODEL, MODEL_DIR, PARAMS, GRID_PARAMS,
                    N_QUANTILES_PREDICT, GRIDSEARCH_CV, 
                    NUM_CORES__MODEL, NUM_CORES__GRID_SEARCH)

def split_train_test(df):
    """Split dataset to train/test based on the predefined $train_period period."""

    # Get target and features
    x = df.drop(columns=['obs','station_names'])
    y = df['obs']

    # Split dataset
    mask_train_set = x.index.year.isin(TRAIN_PERIOD)
    mask_test_set  = x.index.year.isin(TEST_PERIOD)
    x_train, y_train = x[mask_train_set], y[mask_train_set]
    x_test, y_test   = x[mask_test_set], y[mask_test_set]

    # Set training period
    print("Split dataset...")
    print("Train: {}".format(x_train.index.year.unique().to_list()))
    print("Test:  {}".format(x_test.index.year.unique().to_list()))
    print()
    return x_train, y_train, x_test, y_test

def build_model():
    """Building model""" 

    # Test input model
    if MODEL not in ['ngb', 'qrf', 'xtr']:
        sys.exit(f'Model {MODEL} is not supported. Check config file.')

    # Data preprocessor
    preprocessor = ColumnTransformer([('orog', GetOrographyError(), ["alt", "orog"])], remainder='passthrough')
    # Model estimator
    if MODEL=="qrf":
        get_model = QuantileRandomForestRegressor(
            nthreads=NUM_CORES__MODEL, 
            random_state=42
            )
    elif MODEL=="xtr":
        get_model = QuantileExtraTreesRegressor(
            nthreads=NUM_CORES__MODEL, 
            random_state=42
            )
    elif MODEL=="ngb":
        get_model = NGBRegressor(
	        Base=DecisionTreeRegressor(criterion='friedman_mse', random_state=42),
            Dist=Normal, 
            Score=CRPScore,
            random_state=42
            )

    # Data-model pipeline
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

def run_ngb_estimator(hpar, data):
    """Runing estimator for Natural Gradient Boosting."""

    # Input dataset
    X_train, y_train, X_val, y_val = data
    # Build model
    model = build_model()
    # Input grid-search values to dicitionary 
    dic_params = {param:value for param, value in zip(GRID_PARAMS[MODEL].keys(), hpar)}
    # Set model parameters
    model.set_params(**dic_params)
    # estimator 
    model.fit(X_train, y_train)
    # early stopping not work correctly. Problem with loss_eval values (@patrik)
    # model.fit(X_train, y_train, model__X_val=X_val.to_numpy(), model__Y_val=y_val, model__early_stopping_rounds=10)
    # Get (normal) distribution parameters prediction 
    X_val_ = model[:-1].transform(X_val)
    y_dist = model['model'].pred_dist(X_val_) 
    # Evaluation
    score = crps_gaussian_score(y_val, loc=y_dist.params['loc'], scale=y_dist.params['scale'])
    return score

def run_tree_estimator(hpar, data):
    """Runing estimator for Quantile Random Forest and Quantile Extra Trees."""

    # Input dataset
    X_train, y_train, X_val, y_val = data
    # Build model
    model = build_model()
    # Input grid-search values to dicitionary 
    dic_params = {param:value for param, value in zip(GRID_PARAMS[MODEL].keys(), hpar)}
    # Set model parameters
    model.set_params(**dic_params)
    # Prepare train/validation datasets
    X_train_ = model[:-1].fit_transform(X_train)
    X_val_  = model[:-1].transform(X_val)
    y_train_ = y_train.to_numpy()
    # Fit model
    model['model'].fit(X_train_, y_train_)
    # Predict model
    nq = N_QUANTILES_PREDICT + 1
    y_pred = model['model'].predict(X_val_, np.arange(1/nq, nq/nq, 1/nq))
    # Evaluate model 
    score = crps_ensemble_score(y_val, y_pred)
    return score

def tune_model():
    """Grid search model parameters using the cross-validation method."""

    # Load dataset
    data = read_dataset()
    # Train/test datasets
    X_train, y_train, X_test, y_test = split_train_test(data)
    # Build model
    model = build_model()
    # Grid parameters initialization
    grid_pairs = list(itertools.product(*GRID_PARAMS[MODEL].values()))
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
        # Run model grid-search with multiprocessing
        pool = Pool(processes=NUM_CORES__GRID_SEARCH)
        if MODEL == 'ngb':
            gs_score = pool.map_async(partial(run_ngb_estimator, data=(X_train_, y_train_, X_val, y_val)), grid_pairs).get()
        elif MODEL in ['qrf', 'xtr']:
            gs_score = pool.map_async(partial(run_tree_estimator, data=(X_train_, y_train_, X_val, y_val)), grid_pairs).get()
        pool.close()
        pool.join()
        gs_scores += gs_score
    # Mean CV-grid-search score 
    gs_scores = gs_scores/nfolds

    # Best model parameters
    best_grid_pair = grid_pairs[np.argmin(gs_scores)]
    best_model_params = {param:value for param, value in zip(GRID_PARAMS[MODEL].keys(), best_grid_pair)}
    print("Best GridSearchCV model parameters:\n{}".format(best_model_params))
    print("GridSearchCV score: {:.2f}".format(gs_scores[np.argmin(gs_scores)]))
    print()
    # Tune best model
    model.set_params(**best_model_params)
    # train/test 
    X_train_ = model[:-1].fit_transform(X_train)
    y_train_ = y_train.to_numpy()
    model['model'].fit(X_train_, y_train_)
    # Save model
    file_model = os.path.join(MODEL_DIR, f"model_{MODEL}_{TARGET}_ff{LEADTIME}_{STATION}_{TRAIN_PERIOD[0]}.joblib")
    print("Save model to: {}".format(file_model))
    joblib.dump(model, file_model)

def train_model():
    # Read dataset
    data = read_dataset()
    # Get training data
    X_train, y_train, _, _ = split_train_test(data)
    # Build model
    model = build_model()
    model.set_params(**PARAMS[MODEL])
    # Fit model
    print("Train model: {}".format(MODEL))
    X_train_ = model[:-1].fit_transform(X_train)
    y_train_ = y_train.to_numpy()
    model["model"].fit(X_train_, y_train_)
    # Save model
    file_model = os.path.join(MODEL_DIR, f"model_{MODEL}_{TARGET}_ff{LEADTIME}_{STATION}_{TRAIN_PERIOD[0]}.joblib")
    print("Save model to: {}".format(file_model))
    joblib.dump(model, file_model)

def test_model(use_approx=False):
    # Read dataset
    data = read_dataset()
    # Get testing data
    _, _, X_test, y_test = split_train_test(data)
    # Load trained model
    file_model = os.path.join(
        MODEL_DIR, 
        f"model_{MODEL}_{TARGET}_ff{LEADTIME}_{STATION}_{TRAIN_PERIOD[0]}.joblib")
    print("Load trained model: {}".format(file_model))
    model = joblib.load(file_model)
    # Get model features
    column_names = get_model_features(model['preproc'])
    # Get quantile prediction
    nq = N_QUANTILES_PREDICT + 1
    quantiles = np.arange(1/nq, nq/nq, 1/nq)
    # Make prediction
    X_test_  = pd.DataFrame(model[:-1].transform(X_test), columns=column_names)

    if MODEL == 'ngb':
        # Get normal-distribution parameters prediction 
        y_dist = model['model'].pred_dist(X_test_)
        results = pd.DataFrame({
            'station': data.loc[X_test.index.unique(), 'station_names'].values,
            #'obs': y_test.values,
            'mean': y_dist.params['loc'], 
            'std': y_dist.params['scale']
            }, index=X_test.index)
        # Evaluate prediction
        crps_score_model = crps_gaussian_score(
            y_test, loc=y_dist.params['loc'], scale=y_dist.params['scale'])
        # Backward square-root transformation for precipitation
        #if TARGET=='prec24':
        #    y_pred = scipy.stats.norm.ppf(quantiles, loc=y_dist.params['loc'], scale=y_dist.params['scale'])
        #    y_pred = y_pred ** 2
        #    results = pd.concat([
        #        pd.DataFrame({'station': data.loc[X_test.index.unique(), 'station_names'].values}, index=X_test.index),
        #        pd.DataFrame(y_pred.round(4), columns=quantiles.round(3), index=X_test.index)
        #        ], axis=1)

    elif MODEL in ['qrf', 'xtr']:
        y_pred = model['model'].predict(X_test_, quantiles, use_approx=use_approx) 
        # Backward square-root transformation for precipitation
        results = pd.concat([
            pd.DataFrame({'station': data.loc[X_test.index.unique(), 'station_names'].values}, index=X_test.index),
            pd.DataFrame(y_pred.round(4), columns=quantiles.round(4), index=X_test.index)
            ], axis=1)
        # Evaluate prediction
        crps_score_model = crps_ensemble_score(y_test, y_pred)
        # Backward square-root transformation for precipitation
        #if TARGET=='prec24':
        #    results = pd.concat([results['station'], results.drop(columns='station')**2], axis=1)

    # Show prediction scores:
    # baseline reference
    y_loc_base   = X_test_[TARGET+'_mean']
    y_std_base = np.sqrt(X_test_[TARGET+'_var'])
    crps_score_baseline = crps_gaussian_score(y_test, loc=y_loc_base, scale=y_std_base)
    # For precipitation (not gaussian)
    #idx_nonzero_scale = (y_std_base>0).values
    #crps_score_baseline = crps_gaussian_score(
    #    y_test[idx_nonzero_scale], loc=y_loc_base[idx_nonzero_scale], scale=y_std_base[idx_nonzero_scale])

    print("\nCRPS raw-forecast: {:.2f}".format(crps_score_baseline))
    print("CRPS {}-model: {:.2f}".format(MODEL, crps_score_model))
    print("Skill-Score: {:.2f}".format((crps_score_baseline - crps_score_model) / crps_score_baseline))

    # Save prediction
    if STATION=='all':
        file_results = os.path.join(MODEL_DIR, f"pred_{MODEL}_{TARGET}_ff{LEADTIME}_{TRAIN_PERIOD[0]}.csv")
    else:
        file_results = os.path.join(MODEL_DIR, f"pred_{MODEL}_{TARGET}_ff{LEADTIME}_{STATION}_{TRAIN_PERIOD[0]}.csv")
    print("\nSave prediction to: {}".format(file_results))
    results.reset_index().to_csv(file_results, index=False)
