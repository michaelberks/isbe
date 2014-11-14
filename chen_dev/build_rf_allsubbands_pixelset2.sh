#!/bin/bash
#
matlab -nodisplay -nojvm -r "build_rf_allsubbands_pixelset2($JOB_ID, $SGE_TASK_ID, 2e5, 20, 8)" > logs/build_rf$JOB_ID.$SGE_TASK_ID.log