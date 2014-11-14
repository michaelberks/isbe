#!/bin/bash
matlab -nodisplay -nojvm -r "compute_rf_pair_classification_mb($SGE_TASK_ID, 'abnormals', 2)" > logs/compute_rf_pair_classification$JOB_ID.$SGE_TASK_ID.log
