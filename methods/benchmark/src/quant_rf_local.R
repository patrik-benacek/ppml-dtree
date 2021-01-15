#!/usr/bin/env Rscript                                                                                                                                                                       

##-----------------------------------------------------------------------------
## Implementation of quantile regression forest approach to post-processing
## proposed by Taillardat et al (2016, Monthly Weather Review)
## partly based on R code provided by Maxime Taillardat
## here: local = one model for every station

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
exp_name <- "qrf_local"
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
#data$orogerr <- data$alt - data$orog
data$orog <- NULL
data$alt <- NULL

# Time features
#data$month <- month(data$date)

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
print("")
print("-------------------")
print("Build model")
print("-------------------")
t1 <- Sys.time()
qt_levels <- seq(1/51, 50/51, by = 1/51)
qts_save <- matrix(NA, 
                   nrow = nrow(data_eval_all), 
                   ncol = length(qt_levels))

stations <- unique(data$station)
for(this_station in stations){
  # progress indicator
  progind <- which(stations == this_station)
  if(progind %% 10 == 0){
    cat(progind, "of", length(stations), "started at", paste(Sys.time()), "\n")
  }
  
  data_train <- subset(data_train_all, station == this_station)
  data_eval  <- subset(data_eval_all, station == this_station)
                      
  # need to delete cases with NA obs in eval data to avoid error in qrf package
  # (even though not needed in predict function...)
  data_eval$obs <- NULL
  
  # omit stations with too few observation to successfully estimate model
  if(nrow(data_train) <= 10){next}
  
  # estimate model on training data
  qrF_model <- quantregForest(x = data_train[, !(names(data_train) %in% metadata_vars)], 
                              y = data_train$obs,
                              nthreads=4)
  
  # compute quantiles on evaluation data
  qrF_prediction <-   predict(qrF_model,
                              newdata = data_eval[, !(names(data_eval) %in% metadata_vars)], 
                              what = qt_levels,
                              all = TRUE)
  
  # position of fc cases for eval data
  # to save qrF prediction to row in qts_save
  qts_save[which(data_eval_all$station == this_station),] <- qrF_prediction
}

t2 <- Sys.time()

print(paste0("Computational time: ", round(t2-t1, 1), "s"))

print("")
print("-------------------")
print("Evaluation")
print("-------------------")
# Non-defined values detection
ind_noNAinRow <- which(!apply(qts_save, 1, anyNA))
ind_use <- intersect(which(!is.na(data_eval_all$obs)), ind_noNAinRow)
print(paste("All NA cases in evaluation dataset:", dim(qts_save)[1] - length(ind_use)))
print(paste("NA cases due to model estimation:", dim(qts_save)[1] - length(ind_use) - sum(is.na(data_eval_all$obs))))
print("")

print("CRPS scores:")
score_raw = crps_norm(data_eval_all$obs, mean = data_eval_all$t2m_mean, sd = data_eval_all$t2m_var)
score_qrf = crps_sample(y = data_eval_all$obs[ind_use],
                        dat = qts_save[ind_use,])

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

# nodesize = 10, ntree = 250, --> mean CRPS 0.9964
# nodesize = 5, ntree = 250, --> mean CRPS 0.9906
# nodesize = 5, maxnodes = 20, ntree = 250, --> mean CRPS 1.058
# all to default: --> mean CRPS 0.9972
# nodesize = 5, ntree = 500 --> mean CRPS 0.9882
# nodesize = 15, ntree = 500 --> mean CRPS 1.0083
# nodesize = 5, ntree = 1000 --> 0.9875
