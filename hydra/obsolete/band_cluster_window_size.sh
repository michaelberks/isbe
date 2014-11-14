#!/bin/bash
#
matlab -nodisplay -nojvm -r "band_cluster_pyramid_window_size($SGE_TASK_ID, 'normal512')" > logs/window_size$SGE_TASK_ID.log