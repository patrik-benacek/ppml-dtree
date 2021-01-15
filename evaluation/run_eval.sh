#!/bin/bash
leadtime=24
start_train=2018
variable="prec24"
reference_models="Raw-Fcst"
normal_dist_models="EMOS-global EMOS-local NGBoost"
sample_dist_models="QRF XTR"

leval=TRUE
lplot=TRUE

if [ $leval == "TRUE" ]; then
#----------------------------------
# Prepare ScoreFiles for evaluation
#----------------------------------
# Reference Raw-Forec
echo "Prepare ScoreFiles for Raw forecasts"
for model in $reference_models; do
    src/calc_scores.R results/prediction/pred_${model}_${variable}_ff${leadtime}h_2019_2019.csv normal
    ln -s eval_scores_${model}_${variable}_ff${leadtime}h_2019_2019.Rdata results/eval_scores_${model}_${variable}_ff${leadtime}h_${start_train}_2019.Rdata
done

echo "Prepare ScoreFiles for Parametric models"
# Parametric models (normal cdf)
for model in $normal_dist_models; do
    src/calc_scores.R results/prediction/pred_${model}_${variable}_ff${leadtime}h_${start_train}_2019.csv normal
done

echo "Prepare ScoreFiles for Non-parametric models"
# Non-parameteric models (sample cdf)
for model in $sample_dist_models; do
    src/calc_scores.R results/prediction/pred_${model}_${variable}_ff${leadtime}h_${start_train}_2019.csv sample
done
fi

if [ $lplot == "TRUE" ]; then
#----------------------------------
# Run Visualisation
#----------------------------------
echo "Run visualisation ..."
src/plot_scores.R results/eval_scores_*_${variable}_ff${leadtime}h_${start_train}_2019.Rdata
fi

echo "Finish."

