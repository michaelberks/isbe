#!/bin/bash
#
matlab -nodisplay -nojvm -r "knn_contralateral_region($SGE_TASK_ID, 25)" > logs/knn_dist_map$SGE_TASK_ID.log