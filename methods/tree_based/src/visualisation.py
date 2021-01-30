import pandas as pd
import numpy as np
import joblib
import sys
import os

from sklearn.inspection import permutation_importance
from model import read_dataset, split_train_test
import matplotlib.pyplot as plt

sys.path.append('config')
from config import (MODEL, MODEL_DIR, TRAIN_PERIOD, NUM_CORES__MODEL)

def run_permutation_importance(n_features=10):
    """Permutation importance measures the mean absolute error of
    0.5 quantile prediction (median) on the test set. Visualisation of 
    the first most important n_features.
    Input:
        * n_features: num. of most important features to plot 
    """
    # Read dataset
    data, expname = read_dataset()
    # Get testing data
    _, _, X_test, y_test = split_train_test(data)
    # Get trained model
    file_model = os.path.join(MODEL_DIR, f"model_{MODEL}_{expname}_{TRAIN_PERIOD[0]}.joblib")
    print("Load trained model: {}".format(file_model))
    model = joblib.load(file_model) 
    X_test_ = model[:-1].fit_transform(X_test)

    # Get permutation importance for the first most important n_features 
    result = permutation_importance(
        model['model'], X_test_, y_test, 
        scoring='neg_mean_absolute_error', 
        n_repeats=10, random_state=42,
        n_jobs=NUM_CORES__MODEL)
    sorted_idx = result.importances_mean.argsort()[-n_features:]

    # Plot permutation importance
    print("Save plot to: {}".format(MODEL_DIR))
    fig, ax = plt.subplots()
    ax.boxplot(result.importances[sorted_idx].T,
            vert=False, labels=X_test.columns[sorted_idx])
    ax.set_title("Permutation Importances (test set)")
    fig.tight_layout()
    plt.savefig(
        os.path.join(MODEL_DIR, f'feature_perimp_{MODEL}_{expname}_{TRAIN_PERIOD[0]}.png'), 
        dpi = 600,
        bbox_inches = 'tight'
    )
    plt.close()