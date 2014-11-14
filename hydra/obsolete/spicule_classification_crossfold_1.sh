#!/bin/bash
#
matlab -nodisplay -nojvm -r "spicule_classification_crossfold($SGE_TASK_ID, 10, 6, 1, 'all', 'dt')" > logs/spic_cross_1_$SGE_TASK_ID.log