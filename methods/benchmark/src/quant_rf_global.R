#!/usr/bin/env Rscript                                                                                                                                                                       

##-----------------------------------------------------------------------------
## Implementation of quantile regression forest approach to post-processing
## proposed by Taillardat et al (2016, Monthly Weather Review)
## partly based on R code provided by Maxime Taillardat
##-----------------------------------------------------------------------------
## * chosen tuning parameters: ntree=1000, nodesize=5, mtry=25

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
exp_name <- "qrf_global"
leadtime <- paste0("ff", leadtime_hours,"h")
train_date <- c(paste0(start_train_year, "-01-01"),"2019-01-01")
eval_date <- c("2019-01-01","2019-12-31")

# Summary
print("===========================================================")
print(paste0("Run benchmark experiment: ", exp_name))
print("===========================================================")
print("-------------------")
print("Setting")
print("-------------------")
print(paste0("Leadtime: ", args[2], "h"))
print(paste0("Train period: ", train_date[1], " to ", train_date[2]))
print(paste0("Test period : ", eval_date[1], " to ", eval_date[2]))

in_data_dir <- "/home/patrik/Work/czechglobe/TIGGE/data_preproc/generate_dataset/data"
out_data_dir <- "/home/patrik/Work/czechglobe/TIGGE/postproc_benchmark/results"

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

data_train <- data_train_all
data_eval <- data_eval_all

#-------------------------
# Build model
#-------------------------
print("")
print("-------------------")
print("Build model")
print("-------------------")
t1 <- Sys.time()
qt_levels <- seq(1/51, 50/51, by = 1/51)
qts_save <- matrix(NA, 
                   nrow = nrow(data_eval_all), 
                   ncol = length(qt_levels))

# model (default: ntree=500, )
# based on hyperparameter tuning (ntree=1000, nodesize=5, mtry=25)
qrF_model <- quantregForest(x = data_train[, !(names(data_train) %in% metadata_vars)], 
                            y = data_train$obs,
                            ntree = 1000, 
                            nodesize = 5, 
                            mtry = 25,
                            nthreads=4)

# compute quantiles on evaluation data
qrF_prediction <-   predict(qrF_model,
                            newdata = data_eval[, !(names(data_eval) %in% metadata_vars)],
                            what = qt_levels,
                            all = TRUE)

print(qrF_model)

t2 <- Sys.time()

print(paste0("Computational time: ", round(t2-t1, 1), "s"))

print("")
print("-------------------")
print("Evaluation")
print("-------------------")
print("CRPS scores:")
score_raw = crps_norm(data_eval$obs, mean = data_eval$t2m_mean, sd = data_eval$t2m_var)

ind_noNAinRow <- which(!apply(qrF_prediction, 1, anyNA))
ind_use <- intersect(which(!is.na(data_eval$obs)), ind_noNAinRow)

score_qrf <- crps_sample(y = data_eval$obs[ind_use],
                         dat = qrF_prediction[ind_use,])

print(paste0("--> Raw Forecast: "  , round(mean(score_raw, na.rm=TRUE), 3)))
print(paste0("--> ", exp_name, ": ", round(mean(score_qrf, na.rm=TRUE), 3)))
print("")

print("-------------------")
print("CRPS score summary")
print("-------------------")
print("--> Raw Forecast:") 
print(summary(score_raw))
print("")
print(paste0("--> ", exp_name, ":"))
print(summary(score_qrf))
print("")

# Variable importance
feat_imp <- importance(qrF_model, scale=TRUE)
df_imp <- data.frame(feature = rownames(feat_imp), importance = feat_imp, row.names = NULL)
colnames(df_imp)[2] <- "importance"

df_imp_sort <- df_imp %>% 
  arrange(importance) 

# Save feature importance table
write.csv(df_imp_sort,
          file = file.path(out_data_dir, paste0("feat_imp_", exp_name, "_", exp_suffix, ".csv")),
          row.names = FALSE)

# Plot feature importance
pl <- df_imp_sort %>% 
  tail(15) %>% 
  ggplot(aes(x=importance, y=factor(feature, levels=feature))) +
  geom_col() + 
  theme_bw()

# Save feat-imp
ggsave(filename = file.path(out_data_dir, paste0("feat_imp_", exp_name, "_", exp_suffix, ".png")), pl)

#----------------------
# Save results
#----------------------
write.csv(qrF_prediction,
          file = file.path(out_data_dir, paste0(exp_name, "_", exp_suffix, ".csv")),
          row.names = FALSE)

print("Finished Successfully")
