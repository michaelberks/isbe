#!/bin/bash
#
matlab -nodisplay -nojvm -r "band_cluster_pyramid2($SGE_TASK_ID, 'normal_2')" > logs/band_cluster_p_normal$SGE_TASK_ID.log