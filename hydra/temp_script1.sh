#!/bin/bash
#
matlab -nodisplay -nojvm -r "temp_script1($SGE_TASK_ID)" > logs/temp_script1_$SGE_TASK_ID.log