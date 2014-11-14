#!/bin/bash
#
matlab -nodisplay -nojvm -r "band_cluster_pyramid($SGE_TASK_ID, 'normal512', 256, 11, 5)" > logs/band_cluster_pyr$SGE_TASK_ID.log