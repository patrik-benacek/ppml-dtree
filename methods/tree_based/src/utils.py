"""
Utility functions for NWP_PP.

Author: Patrik Benacek 
"""
import os
from scipy.stats import norm
import numpy as np
np.random.seed(42)
from datetime import datetime
import pandas as pd
import pickle
import matplotlib.pyplot as plt

def create_results_df(dates, station_ids, means, stds):
    """
    """
    df = pd.DataFrame({
        'date': dates,
        'station_id': station_ids,
        'mean': means,
        'std': stds
        })
    return df


def save_pickle(data_dir, fn, train_dates=['2015-01-01', '2016-01-01'],
                add_current_error=False, current_error_len=1,
                test_dates=['2016-01-01', '2017-01-01']):
    """Load and pickle dataset"""
    sets = get_train_test_sets(
        data_dir, train_dates, test_dates, aux_dict=aux_dict,
        add_current_error=add_current_error, current_error_len=current_error_len
    )
    with open(data_dir + fn, 'wb') as f:
        pickle.dump(sets, f)


def plot_fc(data_set, idx, distr='pdf', preds=None):
    fc = data_set.features[idx, :2] * data_set.scale_factors
    obs = data_set.targets[idx]

    x = np.linspace(fc[0] - 5 * fc[1], fc[0] + 5 * fc[1], 100)
    if distr == 'pdf':
        y = norm.pdf(x, fc[0], fc[1])
    elif distr == 'cdf':
        y = norm.cdf(x, fc[0], fc[1])
    else:
        raise Exception
    plt.plot(x, y, label='raw ensemble')
    plt.axvline(obs, color='red', label='obs')
    if preds is not None:
        p = preds[idx]
        if distr == 'pdf':
            y = norm.pdf(x, p[0], p[1])
        elif distr == 'cdf':
            y = norm.cdf(x, p[0], p[1])
        plt.plot(x, y, label='prediction')
    plt.xlabel('Temperature [C]')
    plt.legend()
    plt.show()

def plot_shap_values(X_train, model, label, features, nfirst=15, lshap=True, savepath=None):
    """
    Shap value detection.
    Input:
	X_train :: training data
	model   :: best model estimator
	label   :: model name [by user]
	nfirst  :: first n most important features is ploted
        lshap   :: True: ensemble tree, False: linear regression
        save_figure :: save figure as png [True/False]
    """
    import numpy as np
    import pandas as pd
    import shap
    import os
    import matplotlib.pyplot as plt

    varX = features
    dfall = pd.DataFrame(index=varX)
    if lshap:
        shaps = shap.TreeExplainer(model).shap_values(X_train)
        dfall = pd.DataFrame(shaps, columns=varX).apply(lambda x: np.mean(np.abs(x)), axis=0)
        dfall = dfall.sort_values(ascending=False).head(nfirst)
    else:
        dfall = pd.Series(model.coef_, index=varX).apply(lambda x: np.abs(x)).sort_values(ascending=False).head(nfirst)

    ax = dfall.sort_values(ascending=True).plot(kind='barh', log=False, stacked=False, rot=0)
    if lshap:
        ax.set_ylabel("Střední hodnota(|SHAP|)")
    else:
        ax.set_ylabel("|\u03B2-koeficient|")

    ax.set_title('Významnost prediktorů pro %s předpověď ' % label)
    ax.set_xlabel("Prediktor")
    #ax.tick_params(axis="x", labelsize=15)
    #ax.tick_params(axis="y", labelsize=15)
    plt.tight_layout()

    if savepath is not None:
        plt.savefig(savepath)
    else:
        plt.show(block=False)

    plt.close()

# This is ultimately slow                                                                                                                                                                    
# see https://stackoverflow.com/questions/51483951/quantile-random-forests-from-scikit-garden-very-slow-at-making-predictions                                                                
# y_pred = qrf.predict(X_test, quantile = 95)                                                                                                                                                
def predict_approx(model, X_test, quantiles=[0.05, 0.5, 0.95]):                                                                                                                              
    """                                                                                                                                                                                      
    Function to predict quantiles much faster than the default skgarden method                                                                                                               
    This is the same method that the ranger and quantRegForest packages in R use                                                                                                             
    Output is (n_samples, n_quantiles) or (n_samples, ) if a scalar is given as quantiles                                                                                                    
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
