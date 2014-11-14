#!/bin/bash
#
matlab -nodisplay -nojvm -r "band_cluster_pyramid_cls($SGE_TASK_ID, 'mass_2')" > logs/band_cluster_cls_mass$SGE_TASK_ID.log