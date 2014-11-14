#!/bin/bash
#
matlab -nodisplay -nojvm -r "band_cluster_dual_tree($SGE_TASK_ID, 'normal512', 256, 5, 0)" > logs/dual_tree$SGE_TASK_ID.log