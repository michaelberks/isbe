#!/bin/bash   
#$ -N build_rf_contralateral27 
#$ -e /home/zchen/matlab_code/trunk/chen_dev/error/$JOB_NAME_$JOB_ID_$TASK_ID.e
#$ -o /home/zchen/matlab_code/trunk/chen_dev/error/$JOB_NAME_$JOB_ID_$TASK_ID.o 
matlab -nodisplay -nojvm -r "build_rf_contralateral_image_linux(27)" > logs/build_rf_contralateral_image_linux$JOB_ID.$SGE_TASK_ID.log