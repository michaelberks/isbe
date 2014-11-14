#!/bin/bash
#
matlab -nodisplay -nojvm -r "compute_rad_maps_batch($DATA_TYPE)" > logs/compute_rad_maps_batch$JOB_ID.log$SGE_TASK_ID