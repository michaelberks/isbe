#!/bin/bash
#
matlab -nodisplay -nojvm -r "build_rf_radial_class" > logs/radial_rf_$JOB_ID.log$SGE_TASK_ID