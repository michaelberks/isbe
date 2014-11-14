#!/bin/bash
#
matlab -nodisplay -nojvm -r "build_rf_pixel" > logs/build_rf_pixel$JOB_ID.log$SGE_TASK_ID