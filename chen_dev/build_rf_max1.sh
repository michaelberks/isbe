#!/bin/bash
#
matlab -nodisplay -nojvm -r "build_rf_max1($JOB_ID, $SGE_TASK_ID, 2e5, 20, 3)" > logs/build_rf$JOB_ID.$SGE_TASK_ID.log