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
import time                                                                                                                                                                                  
import utils
import yaml

# Setting
model        = 'ngboost'
leadtime     = '24'       # 24h/240h
start_train  = '2018'      # 2015/2018
target       = 'prec24'    # t2m/prec24

# Experimental setting
# * False: load model if exist
# * True : fit new model on dataset
lrerun_model = True

train_date = ("%s-01-01"%start_train, "2019-01-01") 
eval_date  = ("2019-01-01","2019-12-31")                                                                                                                                                     
suffix = '{}_ff{}h_{}_{}'.format(target, leadtime, train_date[0][:4], train_date[1][:4])

in_data_dir  = "/home/patrik/Work/czechglobe/TIGGE/data_preproc/generate_dataset/data"
out_data_dir = os.path.join("/home/patrik/Work/czechglobe/TIGGE/methods/tree_ensemble/results", "ff%sh"%leadtime)

# Open model configuration file
config = yaml.safe_load(open("config/config_%s.yml"%target))

# Initialize model parameter
model_params = config["%sh"%leadtime][model]['model']
learner_params = config["%sh"%leadtime][model]['learner']

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
print("--> Params of {}: {}".format('learner', learner_params))
print()

#------------------------                                                                                                                                                                    
# GET DATA                                                                                                                                                                                   
#------------------------                                                                                                                                                                    
# Read data from CSV file                                                                                                                                                                    
data = pd.read_csv(os.path.join(in_data_dir,"data_wmeta_{}_ff{}h_2015_2019.csv".format(target, leadtime)), parse_dates=True, index_col="date")

# Metadata attributes 
metadata_vars = ["station", "lon", "lat"]
                                                                                                                                                                                             
# Model orography error                                                                                                                                                                      
data["orogerr"] = data.alt - data.orog                                                                                                                                                       
data.drop(columns=["orog", "alt", "station_names"], inplace=True)
                                                                                                                                                                                             
# Remove rows with NA (in observations)                                                                                                                                                      
data.dropna(inplace=True)

#------------------------
# PREPARE DATA FOR ML
#------------------------
X = data.drop(columns=['obs'] + metadata_vars)
y = data['obs']

print("Dataset dimension:")
print("--> X {} and Y {}".format(X.shape, y.shape))

X_train, y_train = X[train_date[0]:train_date[1]], y[train_date[0]:train_date[1]]
X_test, y_test   = X[eval_date[0]:eval_date[1]], y[eval_date[0]:eval_date[1]]

print("--> Train {} and Test {}".format(X_train.shape, X_test.shape))
print()

# Convert dataframe to numpy array
X_train, X_test = X_train.to_numpy(), X_test.to_numpy()
y_train, y_test = y_train.to_numpy(), y_test.to_numpy()

#------------------------------
# FIT MODEL 
#------------------------------
from ngboost import NGBRegressor
from sklearn.tree import DecisionTreeRegressor
from ngboost.scores import CRPScore
from ngboost.distns import Normal
import properscoring as ps
from sklearn.externals import joblib

# Read model if exist and flag is not set to rerun
file_model = os.path.join(out_data_dir, "model_{}_{}.pkl".format(model, suffix))
if os.path.exists(file_model) and not lrerun_model:
    print("Load existing model ...")
    estimator = joblib.load(file_model)
    skip_save_model = True
else:
    # Base learner - tree: hyperparameters based on the best from fpp01
    learner = DecisionTreeRegressor(**learner_params, random_state=42)
    estimator = NGBRegressor(
        **model_params,
        Dist=Normal,                                                                                                                                                                         
        Score=CRPScore, 
        Base=learner
    ) 
    print("Fit model ...")
    # Fit model
    estimator.fit(X_train, y_train)
    skip_save_model = False

print("Run prediction ...")
# Get prediction
y_preds = estimator.predict(X_test)    # point (mean) prediction
y_dists = estimator.pred_dist(X_test)  # distribution parameters (normal) prediction

# Transformation for precipitation
if target=='prec24':
    param_loc = y_dists.params['loc']
    #param_loc[param_loc<0] = 0 # negative to zero
    param_loc = param_loc**2
    param_scale = y_dists.params['scale']
    #y_preds[y_preds<0] = 0 # negative to zero
    y_preds = y_preds**2
    y_test = y_test**2
else:
    param_loc = y_dists.params['loc']
    param_scale = y_dists.params['scale']

print("Save model results ...")
# Save prediction
from utils import create_results_df
res_df = data[eval_date[0]:eval_date[1]].reset_index()
#res_df['station_names'] = res_df.station_names.apply(lambda x: x.strip().replace(" / ", "/"))
res = create_results_df(res_df.date, res_df.station, param_loc, param_scale)
res.to_csv(os.path.join(out_data_dir, "pred_{}_{}.csv".format(model, suffix)), index=False)

# Save model if loaded
if not skip_save_model:
    joblib.dump(estimator, os.path.join(out_data_dir, "model_{}_{}.pkl".format(model, suffix)))

#------------------------------
# MODEL EVALUATION
#------------------------------
def plot_feature_imp():
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    df_loc = pd.DataFrame({'feature':X.columns, 'importance':feature_importance_loc}).sort_values('importance',ascending=False)
    df_scale = pd.DataFrame({'feature':X.columns, 'importance':feature_importance_scale}).sort_values('importance',ascending=False)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20,8))
    fig.suptitle("Feature importance plot for distribution parameters", fontsize=17)
    sns.barplot(x='importance',y='feature',ax=ax1,data=df_loc[:15], color="skyblue").set_title('loc param')
    sns.barplot(x='importance',y='feature',ax=ax2,data=df_scale[:15], color="skyblue").set_title('scale param')
    plt.savefig(os.path.join(out_data_dir, 'feature_imp_{}_{}.png'.format(model, suffix)))
    plt.close()

print("Model evaluation ...")
# Calculate CRPScore
import properscoring as ps
crps_score = ps.crps_gaussian(y_test, mu=param_loc, sig=param_scale).mean()
print("--> CRPScore: {:.4f}".format(crps_score))
print()

## Feature importance for loc and scale trees
feature_importance_loc = estimator.feature_importances_[0]
feature_importance_scale = estimator.feature_importances_[1]

print("Plot feature importance ...")
plot_feature_imp()
print("Plot SHAP importance ...")
from utils import plot_shap_values 
plot_shap_values(X_train, estimator, label=model.upper(), features=X.columns, lshap=True, savepath=os.path.join(out_data_dir, 'shap_imp_{}_{}.png'.format(model, suffix)))

print("Finish")
