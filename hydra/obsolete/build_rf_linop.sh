#!/bin/bash
#
matlab -nodisplay -nojvm -r "build_rf_linop($SGE_TASK_ID, 2e5, 13, 3, 5, 8)" > logs/build_rf_linop$SGE_TASK_ID.log