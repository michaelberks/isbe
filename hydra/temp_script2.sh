#!/bin/bash
#
matlab -nodisplay -nojvm -r "temp_script2($SGE_TASK_ID)" > logs/temp_script2_$SGE_TASK_ID.log