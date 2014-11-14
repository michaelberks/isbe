#!/bin/bash
#
matlab -nodisplay -nojvm -r "predict_image_set('model_id', $MODEL_PATH)" > logs/predict_image_set$JOB_ID.log$SGE_TASK_ID