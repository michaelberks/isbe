#!/bin/bash
#

echo " "
echo "Job ran on and started at: "
hostname
date
echo " "

matlab -nodisplay -nojvm -r "compute_k_maps_batch($DATA_TYPE, $METHOD_TYPE)" > logs/compute_k_maps_batch$JOB_ID.log$SGE_TASK_ID

echo " "
echo "Job finished at: "
date