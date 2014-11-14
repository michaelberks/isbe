#!/bin/bash
#
matlab -nodisplay -nojvm -r "combine_hydra_rfs($MODEL_PATH)" > logs/combine_rf_$MODEL_PATH.log