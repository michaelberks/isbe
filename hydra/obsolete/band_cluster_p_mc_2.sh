#!/bin/bash
#
matlab -nodisplay -nojvm -r "band_cluster_pyramid($SGE_TASK_ID, 'mc_2')" > logs/band_cluster_p_mc$SGE_TASK_ID.log