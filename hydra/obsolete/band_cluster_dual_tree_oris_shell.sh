#!/bin/bash
#
matlab -nodisplay -nojvm -r "band_cluster_dual_tree_oris2($SGE_TASK_ID, 'mass_2', 256, 5, 1)" > logs/dual_tree$SGE_TASK_ID.log