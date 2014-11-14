#!/bin/bash
#$ -N compute_rf_pair_classification
#$ -e /home/zchen/matlab_code/trunk/chen_dev/error/$JOB_NAME_$JOB_ID_$TASK_ID.e
#$ -o /home/zchen/matlab_code/trunk/chen_dev/error/$JOB_NAME_$JOB_ID_$TASK_ID.o 
matlab -nodisplay -nojvm -r "compute_rf_pair_classification($SGE_TASK_ID, 'abnormal')" > /home/zchen/matlab_code/trunk/chen_dev/logs/compute_rf_pair_classification$JOB_ID.$SGE_TASK_ID.log
