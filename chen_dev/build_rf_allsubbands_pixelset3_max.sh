#!/bin/bash  
#$ -N build_rf_allsubbands3_max 
#$ -e /home/zchen/matlab_code/trunk/chen_dev/error/$JOB_NAME_$JOB_ID_$TASK_ID.e
#$ -o /home/zchen/matlab_code/trunk/chen_dev/error/$JOB_NAME_$JOB_ID_$TASK_ID.o 
matlab -nodisplay -nojvm -r "build_rf_allsubbands_pixelset3_max($JOB_ID, $SGE_TASK_ID, 2e5, 20, 10)" > logs/build_rf_max3$JOB_ID.$SGE_TASK_ID.log