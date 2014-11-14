#!/bin/bash
#
matlab -nodisplay -nojvm -r "build_vessel_predictor" > logs/build_vessel_predictor_$JOB_ID.log$SGE_TASK_ID