#!/usr/bin/env Rscript                                                                                                                                                                       
#--------------------------------------------------------------------
# Prepare evaluation scores for an experiment.
# Usage: ./calc_scores.R prediction_file[csv] family [normal/sample]
# Author: Patrik Benacek
#----------------------------------------------
rm(list=ls())

suppressMessages(library(lubridate))
suppressMessages(library(scoringRules))

# Input arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("Input arguments: file[prediction] family[normal/sample]", call.=FALSE)
}

file   = args[1]
family = args[2]

# Outpath
datadir_obs <- "/home/patrik/Work/czechglobe/TIGGE/evaluation/results/"
outpath <- datadir_obs

# Read experiment predictions
filename = basename(file)
meta = strsplit(filename, "_")[[1]]

if (length(meta)==6){
  exp_name    = meta[2]
  target      = meta[3]
  fc_time     = meta[4]
  start_train = meta[5]
}else{
  stop("Input filename has not valid name. Use e.g.: pred_EMOS-global_t2m_ff24h_2015_2019.csv") 
}

data_ens <- read.csv(file)
data_ens$date = as.Date(data_ens$date)

# Read observations
data_obs <- readRDS(file.path(datadir_obs, paste0("eval_obs_", target,"_", fc_time,"_2019.RData")))

# Check consistency
if (!nrow(data_ens)==nrow(data_obs)){warning("Num of observation and prediction is different. Use inner merge!")}
mdata <- merge(data_obs, data_ens, by.x = c('date', 'station_id'), by.y = c('date','station_id'))
rm(data_obs)
rm(data_ens)

# Check num of test data
print(paste0("Num. of data for evalution: ", nrow(mdata)))

df_res <- mdata[,c('date', 'station_id')]
df_res$exp <- exp_name

# Calculate scores given observations and the predicted parameters of a local-scale normal distribution
if (family=='normal'){
  df_res$crps <- crps_norm(y = mdata$obs, location = mdata$mean, scale = mdata$std)
  df_res$pit  <- pnorm(mdata$obs, mdata$mean, mdata$std) 
  df_res$ae   <- abs(mdata$obs - mdata$mean)
  
# Calculate scores given observations and draws from the predictive distribution
}else if (family=='sample'){
  max_sample_bin = 51
  fc_matrix <- as.matrix(mdata[,4:ncol(mdata)])
  df_res$crps <- crps_sample(y = mdata$obs, dat = fc_matrix)
  df_res$pit  <- sapply(seq_along(1:nrow(fc_matrix)),
                    function(i) rank(c(mdata$obs[i], fc_matrix[i,]))[1])/max_sample_bin
  df_res$ae   <- abs(mdata$obs - apply(fc_matrix, 1, median))
}else{
  stop(paste("Family", family, "is not defined. Stop."))
}

saveRDS(df_res, file = file.path(outpath, paste0("eval_scores_", exp_name, "_", target, "_", fc_time,"_", start_train,"_2019.Rdata")))
