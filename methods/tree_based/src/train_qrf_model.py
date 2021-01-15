#!/usr/bin/env python3

"""
Tuning Decision Tree ensemble methods (QRF, QXT, NGB) with respect to CRPS scores. 
This is performed by oputna hyperparameter optimization framework:
See https://optuna.org/
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import multiprocessing
import itertools 
import time
import yaml
import sys                                                                                                                                                                                   
                                                                                                                                                                                             
if len(sys.argv)<3:                                                                                                                                                                          
    print("Put target and expname as input arguments!") 
    sys.exit(0)                                                                                                                                                                              
else:                                                                                                                                                                                        
    target = sys.argv[1]
    expname = sys.argv[2]   
    if target not in ['t2m', 'prec24']:
        print(f"{target} is not supported yet.")                                                                                                                                             
        sys.exit(0) 

# Setting
model = 'qrf'
leadtime = 24
start_train = 2018

train_date = ("%s-01-01"%start_train,"2019-01-01")
eval_date  = ("2019-01-01","2019-12-31")
exp_suffix = "{}_{}_{}_{}_{}_{}".format(model, target, expname, leadtime, train_date[0][:4], train_date[1][:4])

in_data_dir = "/home/benacek.p/TIGGE/pp_tree_ensemble_methods/data"
out_data_dir = "/home/benacek.p/TIGGE/pp_tree_ensemble_methods/results/ff%sh"%leadtime

# Open model configuration file                                                                                                                                                              
config = yaml.safe_load(open("config_%s.yml"%target)) 

# Initialize model parameter                                                                                                                                                                 
model_params = config["%sh"%leadtime][model]['model']                                                                                                                                        
num_cores = config["%sh"%leadtime][model]['ncores']   

if not os.path.exists(out_data_dir):                                                                                                                                                         
    os.makedirs(out_data_dir)                                                                                                                                                                
                                                                                                                                                                                             
# Print summary                                                                                                                                                                              
print(30*"=")                                                                                                                                                                                
print("Run {} model for the {}h forecast".format(model, leadtime))                                                                                                                           
print(30*"=")                                                                                                                                                                                
print("Paths:")                                                                                                                                                                              
print("--> Input : %s" % in_data_dir)                                                                                                                                                        
print("--> Output: %s" % out_data_dir)                                                                                                                                                       
print()                                                                                                                                                                                      
print("Time period:")                                                                                                                                                                        
print("--> Train: {} to {}".format(train_date[0][:4], train_date[1][:4]))                                                                                                                    
print("--> Test : {} to {}".format(eval_date[0][:4], eval_date[1][:4]))                                                                                                                      
print()                                                                                                                                                                                      
print("Model settins:")                                                                                                                                                                      
print("--> Params of {}: {}".format(model, model_params))                                                                                                                                    
print("--> Num of cores: {}".format(num_cores))
print()

#------------------------
# GET DATA 
#------------------------
def read_dataset():
    # Read data from CSV file
    data = pd.read_csv(os.path.join(in_data_dir,"data_wmeta_{}_ff{}h_2015_2019.csv".format(target, leadtime)), parse_dates=True, index_col="date").iloc[:, 1:]
    # Model orography error
    data["orogerr"] = data.alt - data.orog
    data.drop(columns=["orog", "alt"], inplace=True)

    # Remove rows with NA (in observations)
    data.dropna(inplace=True)
    return(data)

#------------------------
# PREPARE DATA FOR ML
#------------------------
def get_train_test(data):
    # Metadata attributes
    metadata_vars = ["station_names", "lon", "lat"]

    X = data.drop(columns=['obs'] + metadata_vars)
    y = data['obs']

    print("Dimension of X {} and Y {}".format(X.shape, y.shape))

    X_train, y_train = X[train_date[0]:train_date[1]], y[train_date[0]:train_date[1]]
    X_test, y_test   = X[eval_date[0]:eval_date[1]], y[eval_date[0]:eval_date[1]]

    print("Dimension of TRAIN {} and TEST {}".format(X_train.shape, X_test.shape))
    return(X_train, X_test, y_train, y_test)

#------------------------------
# FINE-TUNE MODEL BY SKLEARN 
#------------------------------
def run_estimator(hpar):
    from skgarden import RandomForestQuantileRegressor
    import properscoring as ps
    # import warnings filter
    from warnings import simplefilter
    # ignore all future warnings
    simplefilter(action='ignore', category=FutureWarning)
    
    # estimator 
    estimator_obj = RandomForestQuantileRegressor(
        n_estimators=hpar[0], 
        max_features=hpar[1], 
        min_samples_leaf=hpar[2], 
        bootstrap=hpar[3],
        random_state=42, 
        n_jobs=1
    )

    estimator_obj.fit(X_train, y_train)

    # Get prediction for particular quantiles
    y_pred = predict_approx(estimator_obj, X_test, quantiles=np.arange(1/51, 51/51, 1/51))
    if target=='prec24':
        crps_score = ps.crps_ensemble(y_test**2, y_pred**2).mean()
    else:
        crps_score = ps.crps_ensemble(y_test, y_pred).mean()

    return((hpar[0], hpar[1], hpar[2], hpar[3], crps_score))

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

# Prepare resulting data frame
def expand_grid(dictionary):
    result = pd.DataFrame([row for row in itertools.product(*dictionary.values())], columns=dictionary.keys())
    return(result)

if __name__=='__main__':

    # read data
    data = read_dataset()
    # get train and test sets 
    X_train, X_test, y_train, y_test = get_train_test(data)
    # pairs for processors
    pool = multiprocessing.Pool(processes=num_cores)
    pairs = itertools.product(*model_params.values())

    scores = []
    t0 = time.time()
    results = pool.map_async(run_estimator, pairs).get()
    for r in results:
        print('%d,%d,%d,%s,%f\n' % r)
        scores.append(r[4])

    # run workers
    pool.close()
    pool.join()
    t1 = time.time()
    print("Total time: {:.2f}s".format(t1-t0))
    print()

    df_grid = expand_grid(model_params)
    df_grid['score'] = np.asarray(scores, dtype=np.float32)
    df_grid.to_csv(os.path.join(out_data_dir, "results_%s.csv" % exp_suffix), index=False)

    print("Finish sucessfully.")
