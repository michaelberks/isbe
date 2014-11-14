#!/bin/bash
#
matlab -nodisplay -nojvm -r "depth_cluster_dual_tree('lecb', 512, 10, $SGE_TASK_ID)" > logs/depth_cluster_dual_tree$SGE_TASK_ID.log