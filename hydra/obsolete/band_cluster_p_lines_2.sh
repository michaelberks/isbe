#!/bin/bash
#
matlab -nodisplay -nojvm -r "band_cluster_pyramid($SGE_TASK_ID, 'lines')" > logs/band_cluster_p_lines$SGE_TASK_ID.log