#!/bin/bash
#
matlab -nodisplay -nojvm -r "band_cluster_dual_tree_all($SGE_TASK_ID, 'normal512', 256, 3, 1, 10)" > logs/cluster_dual_tree$SGE_TASK_ID.log