#!/bin/bash
#
matlab -nodisplay -nojvm -r "add_hydra_rfs($FOREST_JOB, $SGE_TASK_ID)" > logs/combine_rf_$FOREST_JOB.log