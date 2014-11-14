#!/bin/bash
#
matlab -nodisplay -nojvm -r "band_cluster_g_pyramid($SGE_TASK_ID, 'mass_2')" > logs/band_cluster_gp_mass$SGE_TASK_ID.log