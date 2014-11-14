#!/bin/bash
#
matlab -nodisplay -nojvm -r "pyr_synthesis($SGE_TASK_ID, 'normal512')" > logs/pyr_synthesis$SGE_TASK_ID.log