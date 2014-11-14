#!/bin/bash
#
matlab -nodisplay -nojvm -r "classify_image_set($FOREST_JOB, $TEST_IMAGE_DIR)" > logs/class_images_$JOB_ID.log$SGE_TASK_ID