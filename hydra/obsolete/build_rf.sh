#!/bin/bash
#
matlab -nodisplay -nojvm -r "build_rf($SGE_TASK_ID, 2e5, 8, 10)" > logs/build_rf$SGE_TASK_ID.log