#!/bin/bash
#
matlab -nodisplay -nojvm -r "dt_synthesis($SGE_TASK_ID, 'normal512', 'normal512_k10_w1_5_w2_0')" > logs/dt_synthesis$SGE_TASK_ID.log