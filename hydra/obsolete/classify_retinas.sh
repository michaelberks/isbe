#!/bin/bash
#
matlab -nodisplay -nojvm -r "classify_retinas" > logs/classify_retinas_$JOB_ID-$SGE_TASK_ID.log