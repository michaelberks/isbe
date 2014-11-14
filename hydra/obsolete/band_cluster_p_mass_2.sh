#!/bin/bash
#
matlab -nodisplay -nojvm -r "band_cluster_pyramid_lower($SGE_TASK_ID, 'normal_2')" > logs/band_cluster_p_normal_lower$SGE_TASK_ID.log