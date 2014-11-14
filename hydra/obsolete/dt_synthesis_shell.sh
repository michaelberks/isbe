#!/bin/bash
#
matlab -nodisplay -nojvm -r "dt_synthesis_oris($SGE_TASK_ID, 'normal512')" > logs/dt_synthesis$SGE_TASK_ID.log