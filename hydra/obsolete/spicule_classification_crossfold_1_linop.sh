#!/bin/bash
#
matlab -nodisplay -nojvm -r "spicule_classification_crossfold($SGE_TASK_ID, 10, 6, 1, 'all', 'linop')" > logs/spic_cross_1_linop_$SGE_TASK_ID.log