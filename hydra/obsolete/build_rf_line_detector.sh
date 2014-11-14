#!/bin/bash
#
matlab -nodisplay -nojvm -r "build_rf_line_detector" > logs/build_rf_$JOB_ID.log$SGE_TASK_ID