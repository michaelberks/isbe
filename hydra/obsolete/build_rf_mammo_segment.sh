#!/bin/bash
#
matlab -nodisplay -nojvm -r "build_rf_mammo_segment" > logs/build_rf_mammo_segment$JOB_ID.log$SGE_TASK_ID