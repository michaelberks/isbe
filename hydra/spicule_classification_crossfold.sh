#!/bin/bash
#
matlab -nodisplay -nojvm -r "spicule_classification_crossfold($SGE_TASK_ID, 10, 3, 1, 'all')" > logs/spic_cross_1_$SGE_TASK_ID.log
matlab -nodisplay -nojvm -r "spicule_classification_crossfold($SGE_TASK_ID, 10, 4, 1, 'all')" > logs/spic_cross_1_$SGE_TASK_ID.log
matlab -nodisplay -nojvm -r "spicule_classification_crossfold($SGE_TASK_ID, 10, 5, 1, 'all')" > logs/spic_cross_1_$SGE_TASK_ID.log