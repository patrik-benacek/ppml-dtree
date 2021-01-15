#!/usr/bin/env Rscript                                                                                                                                                                       

##-----------------------------------------------------------------------------
## Implementation of quantile regression forest approach to post-processing
## proposed by Taillardat et al (2016, Monthly Weather Review)
## partly based on R code provided by Maxime Taillardat
##-----------------------------------------------------------------------------
## --> Hyperparameter tunning

rm(list=ls())

# Input arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Input arguments: [start_train_period [year]] [forecast leadtime [hours]]", call.=FALSE)
} else if (length(args)==1) {
  # default lead time 
  args[2] = "24"
}

start_train_year = args[1]
leadtime_hours = args[2]

# Load packages
library(quantregForest)
library(scoringRules)
library(lubridate)
library(tidyverse)

# Settings
exp_name <- "qrf_global_tune"
leadtime <- paste0("ff", leadtime_hours,"h")
train_date <- c(paste0(start_train_year, "-01-01"),"2019-01-01")
eval_date <- c("2019-01-01","2019-12-31")

# Summary
print("===========================================================")
print(paste0("Tune model: ", exp_name))
print("===========================================================")
print("-------------------")
print("Setting")
print("-------------------")
print(paste0("Leadtime: ", args[2], "h"))
print(paste0("Train period: ", train_date[1], " to ", train_date[2]))
print(paste0("Test period : ", eval_date[1], " to ", eval_date[2]))

in_data_dir <- "/home/patrik/Work/czechglobe/TIGGE/data_preproc/generate_dataset/data"
out_data_dir <- "/home/patrik/Work/czechglobe/TIGGE/postproc_benchmark/results/qrf_tuning"

# Name of experiment
exp_suffix = paste(leadtime, year(ymd(train_date[1])), year(ymd(train_date[2])), sep='_')

# Read station measurements
file_obs = paste0("data_wmeta_", leadtime, "_2015_2019.csv")
data <- read.csv(file.path(in_data_dir, file_obs))
data$date = as.Date(data$date)

#-----------------------
# Feature engineering 
#-----------------------
print("Prepare dataset ...")

# Metadata
metadata_vars <- c("date", "station", "station_names", "lon", "lat", 'obs')
  
# Model orography error
data$orogerr <- data$alt - data$orog
data$orog <- NULL
data$alt <- NULL

# Time features
data$month <- month(data$date)

# Feature selection
#data <- data[, -which(!(names(data) %in% c("obs", "date", "station", "t2m_mean", "t2m_var")))]

#----------------------
# Train/Test datasets
#----------------------
# Train dataset
train_start <- as.Date(paste(train_date[1], "00:00 UTC")) 
train_end <- as.Date(paste(train_date[2], "00:00 UTC")) - days(2)

data_train_all <- subset(data, date >= train_start & date <= train_end)
data_train_all <- data_train_all[complete.cases(data_train_all), ]

# Eval dataset
eval_start <- as.Date(paste(eval_date[1], "00:00 UTC"))
eval_end <- as.Date(paste(eval_date[2], "00:00 UTC"))
eval_dates <- seq(eval_start, eval_end, by = "1 day")

data_eval_all <- subset(data, date >= eval_start & date <= eval_end)

#-------------------------
# Build model
#-------------------------
qt_levels <- seq(1/51, 50/51, by = 1/51)
qts_save <- matrix(NA, 
                   nrow = nrow(data_eval_all), 
                   ncol = length(qt_levels))


# see below for description

fit_qrf <- function(ntree, nodesize, mtry, replace){
  
  data_train <- data_train_all
  data_train <- data_train[complete.cases(data_train),]
  
  data_eval <- data_eval_all
  # need to delete cases with NA obs in eval data to avoid error in qrf package
  # (even though not needed in predict function...)
  data_eval$obs <- NULL
  
  # estimate model on training data
  qrF_model <- quantregForest(x = data_train[, !(names(data_train) %in% metadata_vars)], 
                              y = data_train$obs,
                              nodesize = nodesize,
                              ntree = ntree,
                              mtry = mtry,
                              replace = replace,
                              nthreads=4)
  
  # compute quantiles on evaluation data
  qrF_prediction <-   predict(qrF_model,
                              newdata = data_eval[, !(names(data_eval) %in% metadata_vars)],
                              what = qt_levels,
                              all = TRUE)
  
  ind_noNAinRow <- which(!apply(qrF_prediction, 1, anyNA))
  ind_use <- intersect(which(!is.na(data_eval_all$obs)), ind_noNAinRow)
  
  qrf_crps <- crps_sample(y = data_eval_all$obs[ind_use],
                          dat = qrF_prediction[ind_use,])
  
  savename <- paste0(exp_name, "_", exp_suffix, "_ntree", ntree, "_nodesize", nodesize, "_mtry", mtry, "_repl", replace, ".Rdata")
  fname <- file.path(out_data_dir, savename)
  
  save(qrf_crps, file = fname)
}

# ntree: Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.
#     default: 500
# nodesize: Minimum size of terminal nodes
#     default: 5
# mtry: Number of variables randomly sampled as candidates at each split
#     default: max(floor(ncol(x)/3), 1) = 11
# replace: Should sampling of cases be done with or without replacement?
#     default: 1 (or 0) [only use 1 due to results from local model!]

ntree.try <- c(125,250,500,1000)
nodesize.try <- c(5,10,15,20,50)
mtry.try <- c(10,15,20,25)
replace.try <- c(1)

#ntree.try <- c(1000)
#nodesize.try <- c(5)
#mtry.try <- c(30, 40, 50)
#replace.try <- c(1)

pars <- expand.grid(ntree.try, nodesize.try, mtry.try, replace.try)

nrow(pars) # --> 100 tasks

for (id in seq_along(1:nrow(pars))){
  print(paste("Tune nb.", id))
  print(paste(c("ntree:", "nodesize:", "mtry:", "replace"), paste(pars[id,], sep="_")))
  fit_qrf(ntree = pars[id,1], nodesize = pars[id,2], mtry = pars[id,3], replace = pars[id,4])
}
