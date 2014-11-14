#!/bin/bash
#
matlab -nodisplay -nojvm -r "pt_testbench" > logs/testbench_$JOB_ID-$SGE_TASK_ID.log