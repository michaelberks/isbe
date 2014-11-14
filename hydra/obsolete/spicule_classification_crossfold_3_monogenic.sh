#!/bin/bash
#
matlab -nodisplay -nojvm -r "spicule_classification_crossfold($SGE_TASK_ID, 10, 6, 3, 'all', 'monogenic')" > logs/spic_cross_3_monogenic_$SGE_TASK_ID.log