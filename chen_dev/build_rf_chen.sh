#!/bin/bash
#
matlab -nodisplay -nojvm -r "build_rf($SGE_TASK_ID, 2e4, 20, 10)" > logs/build_rf$SGE_TASK_ID.log