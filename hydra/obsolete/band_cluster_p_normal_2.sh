#!/bin/bash
#
matlab -nodisplay -nojvm -r "band_cluster_pyramid_window_size($SGE_TASK_ID, 'mass_2')" > logs/window_size$SGE_TASK_ID.log