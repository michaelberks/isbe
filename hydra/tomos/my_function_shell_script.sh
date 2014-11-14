#!/bin/bash
#
matlab -nodisplay -nojvm -r "my_function_wrapper($SGE_TASK_ID)" > logs/my_function_$SGE_TASK_ID.log