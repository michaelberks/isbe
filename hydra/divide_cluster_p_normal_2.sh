#!/bin/bash
#
matlab -nodisplay -nojvm -r "divide_cluster_pyramid($SGE_TASK_ID, 'normal_2')" > logs/divide_cluster_p_normal$SGE_TASK_ID.log