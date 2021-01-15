#!/bin/bash                                                                                                                                                                                  
                                                                                                                                                                                             
# Activate eccodes environment                                                                                                                                                               
source /home/patrik/miniconda3/bin/activate nwp_pp

#what='vte_amper_tigge'   # vte/final/vte_amper/vte_amper_tigge
what='vte_amper'   # vte/final/vte_amper/vte_amper_tigge
exp='d2'
                                                                                                                                                                                             
# Run scripts                                                                                                                                                                                
src/run_ngb_pred.py $what &> log_ngb_${what}_$exp.out
src/run_qrf_pred.py $what &> log_qrf_${what}_$exp.out  
src/run_xtr_pred.py $what &> log_xrt_${what}_$exp.out  
                                                                                                                                                                                             
# Deactivate conda env                                                                                                                                                                       
conda deactivate 
