%% -------------------------------------------------------------------------
%Increasing noise experiment: see increasing_noise_experiment.m

% Orientatio prediction...
%Build forests
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[0.25]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[0.50]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="gaussian" NOISE_PARAMS="[0.25]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="gaussian" NOISE_PARAMS="[0.50]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="gaussian" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="gaussian" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="dt" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="dt" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
%Combine forests
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5671'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5672'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5763'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5764'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5675'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5676'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5677'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5678'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5679'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'10210'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10210 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'10211'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10211 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
%Test images are generated and predicted on the fly in the experiment
%script

%% -------------------------------------------------------------------------
%Lines on edge backgrounds: see edge_vs_line_experiment.m

% Line detection...
%Build forests
DECOMP_TYPE="g2d" WIN_SIZE=1 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="edge" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="g2d" WIN_SIZE=3 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="edge" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="dt" WIN_SIZE=1 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="edge" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="dt" WIN_SIZE=3 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="edge" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh

DECOMP_TYPE="g2d" WIN_SIZE=1 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="g2d" WIN_SIZE=3 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="dt" WIN_SIZE=1 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="dt" WIN_SIZE=3 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh

%Combine forests
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10273'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10273 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10274'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10274 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10275'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10275 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10276'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10276 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10277'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10277 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10278'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10278 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10279'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10279 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10280'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10280 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh

%Classify images
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10273'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10290 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10274'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10291 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10275'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10292 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10276'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10293 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10277'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10294 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10278'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10295 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10279'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10296 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10280'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10297 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh

%% ------------------------------------------------------------------------
%Synthetic lines: comparison of representations

%Line detection
DECOMP_TYPE="dt" WIN_SIZE=3 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="g2d" WIN_SIZE=3 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="g12d" WIN_SIZE=3 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="mono" WIN_SIZE=3 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="linop" WIN_SIZE=3 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="haar" WIN_SIZE=3 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="gabor" WIN_SIZE=3 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh

MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10400'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10400 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10401'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10401 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10402'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10402 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10403'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10403 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10833'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10833 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10405'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10405 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'11391'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 11391 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh

% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10400'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10407 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10401'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10408 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10402'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10409 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10403'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10410 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10833'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10834 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10405'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10412 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% MODEL_ROOT="models/synthetic_lines/detection/rf_classification" MODEL_PATH="'10406'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10413 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh

%Orientation prediction
DECOMP_TYPE="dt" WIN_SIZE=3 OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="g2d" WIN_SIZE=3 OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="g12d" WIN_SIZE=3 OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="mono" WIN_SIZE=3 OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="linop" WIN_SIZE=3 OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="haar" WIN_SIZE=3 OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="gabor" WIN_SIZE=3 OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh

MODEL_ROOT="models/synthetic_lines/orientation/rf_regression" MODEL_PATH="'10850'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10850 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/orientation/rf_regression" MODEL_PATH="'10851'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10851 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/orientation/rf_regression" MODEL_PATH="'10852'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10852 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/orientation/rf_regression" MODEL_PATH="'10853'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10853 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/orientation/rf_regression" MODEL_PATH="'10854'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10854 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/orientation/rf_regression" MODEL_PATH="'10855'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10855 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/orientation/rf_regression" MODEL_PATH="'11322'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 11322 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh

% MODEL_ROOT="models/synthetic_lines/orientation/rf_regression" MODEL_PATH="'10850'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10856 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% MODEL_ROOT="models/synthetic_lines/orientation/rf_regression" MODEL_PATH="'10851'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10857 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% MODEL_ROOT="models/synthetic_lines/orientation/rf_regression" MODEL_PATH="'10852'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10858 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% MODEL_ROOT="models/synthetic_lines/orientation/rf_regression" MODEL_PATH="'10853'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10859 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% MODEL_ROOT="models/synthetic_lines/orientation/rf_regression" MODEL_PATH="'10854'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10860 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
% MODEL_ROOT="models/synthetic_lines/orientation/rf_regression" MODEL_PATH="'10855'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10861 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh

%Width...
DECOMP_TYPE="dt" WIN_SIZE=3 OUTPUT_TYPE="width" PREDICTION_TYPE="rf_regression" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="g2d" WIN_SIZE=3 OUTPUT_TYPE="width" PREDICTION_TYPE="rf_regression" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="g12d" WIN_SIZE=3 OUTPUT_TYPE="width" PREDICTION_TYPE="rf_regression" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="mono" WIN_SIZE=3 OUTPUT_TYPE="width" PREDICTION_TYPE="rf_regression" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="linop" WIN_SIZE=3 OUTPUT_TYPE="width" PREDICTION_TYPE="rf_regression" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="haar" WIN_SIZE=3 OUTPUT_TYPE="width" PREDICTION_TYPE="rf_regression" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="gabor" WIN_SIZE=3 OUTPUT_TYPE="width" PREDICTION_TYPE="rf_regression" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh


MODEL_ROOT="models/synthetic_lines/width/rf_regression" MODEL_PATH="'10438'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10438 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/width/rf_regression" MODEL_PATH="'11378'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 11378 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/width/rf_regression" MODEL_PATH="'11379'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 11379 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/width/rf_regression" MODEL_PATH="'11380'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 11380 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/width/rf_regression" MODEL_PATH="'11381'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 11381 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/width/rf_regression" MODEL_PATH="'11382'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 11382 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/width/rf_regression" MODEL_PATH="'11383'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 11383 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh

% MODEL_ROOT="models/synthetic_lines/width/rf_regression" MODEL_PATH="'10438'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/comparing_representations" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10439 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
