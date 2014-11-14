#!/bin/bash
#
matlab -nodisplay -nojvm -r "build_rf_breiman('bg+bar_128/', 4e5, 100, 10)" > logs/build_rf$SGE_TASK_ID.log