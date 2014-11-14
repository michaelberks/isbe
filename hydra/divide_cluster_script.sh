#!/bin/bash
#
matlab -nodisplay -nojvm -r "divide_cluster_wrapper($SGE_TASK_ID, 'mc')" > logs/divide_cluster$SGE_TASK_ID.log