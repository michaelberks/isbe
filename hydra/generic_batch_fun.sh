#!/bin/bash
#
matlab -nodisplay -nojvm -r "generic_batch_fun($SGE_TASK_ID)" > logs/generic_batch_fun$JOB_ID.log$SGE_TASK_ID