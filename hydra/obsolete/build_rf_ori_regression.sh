#!/bin/bash
#
matlab -nodisplay -nojvm -r "build_rf_ori_regression($JOB_ID, $SGE_TASK_ID, 2e5, 10, 20)" > logs/rf_ori_reg$SGE_TASK_ID.log