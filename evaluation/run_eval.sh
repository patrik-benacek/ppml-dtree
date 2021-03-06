#!/bin/bash
leadtime=24        # 24/240
start_train=2018
variable="t2m"
reference_models="Raw-Fcst"

leval=FALSE
lplot=TRUE

#----------------------------------
# Prepare ScoreFiles for evaluation
#----------------------------------
if [ $leval == "TRUE" ]; then
  echo "Prepare Raw forecasts prediction"
  # Prepare raw forecast prediction 
  src/prep_eval_obs_raw.R $variable $leadtime

  # Reference Raw-Forec
  echo "Prepare ScoreFiles for Raw forecasts"
  src/calc_scores.R results/prediction/pred_Raw-Fcst_${variable}_ff${leadtime}h_2019.csv sample
  mv results/eval_scores_Raw-Fcst_${variable}_ff${leadtime}h_2019.Rdata results/eval_scores_Raw-Fcst_${variable}_ff${leadtime}h_${start_train}.Rdata

  echo "Prepare ScoreFiles for Parametric models (normal cdf)"
  src/calc_scores.R ../methods/benchmark/results/pred_EMOS-global_${variable}_ff${leadtime}h_${start_train}.csv normal
  src/calc_scores.R ../methods/benchmark/results/pred_EMOS-local_${variable}_ff${leadtime}h_${start_train}.csv normal
  src/calc_scores.R ../methods/tree_based/models/pred_ngb_${variable}_ff${leadtime}h_${start_train}.csv normal

  echo "Prepare ScoreFiles for Non-parametric models (sample cdf)"
  src/calc_scores.R ../methods/tree_based/models/pred_qrf_${variable}_ff${leadtime}h_${start_train}.csv sample
  src/calc_scores.R ../methods/tree_based/models/pred_xtr_${variable}_ff${leadtime}h_${start_train}.csv sample
fi

#----------------------------------
# Run Visualisation
#----------------------------------
if [ $lplot == "TRUE" ]; then
  echo "Run visualisation ..."
  src/plot_scores.R results/eval_scores_*_${variable}_ff${leadtime}h_${start_train}.Rdata
  rm -f Rplots.pdf
fi

echo "Finish."

