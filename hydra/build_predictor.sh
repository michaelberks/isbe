#!/bin/bash
#
matlab -nodisplay -nojvm -r "build_predictor" > logs/build_predictor_$JOB_ID.log$SGE_TASK_ID