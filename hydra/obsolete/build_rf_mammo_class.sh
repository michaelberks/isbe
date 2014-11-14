#!/bin/bash
#
matlab -nodisplay -nojvm -r "build_rf_mammo_class" > logs/mammo_rf_$JOB_ID.log$SGE_TASK_ID