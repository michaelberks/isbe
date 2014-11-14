# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions
if [ "$USER" == "ptresadern" ]; then
	alias qsub='qsub -m e -M philip.tresadern@manchester.ac.uk'
fi       

# %Generic settings
export USE_NAG=0
export RAND_SEED=""
export PX_PER_MM="100/9"
export THRESH=""
export QUIET=1
export SAVE_TYPE="normal"
export IMAGE_FORMAT=".png"

# %Batch arguments
export NUM_JOBS=20
export CUSTOM_ID=""

# %Directory locations
export BG_PROB_DIR=""
export BLANK_BG_DIR="blank_bgs"
export CLASS_DIR=""
export DATA_ROOT="/san/images/asym/asym/data/"
export DATA_TYPE="abnormals"
export FG_MASK_DIR="vessel_masks"
export FOV_MASK_DIR="fov_masks"
export IMAGE_DIR="images"
export IMAGE_ROOT="retinograms/DRIVE/training"
export LINE_DIR="line_maps"
export MAM_NAMES=""
export MASK_DIR=""
export MODEL_PATH=""
export MODEL_ROOT=""
export ORI_DIR="orientations"
export PECTORAL_DIR=""
export PREDICTION_DIR="predictions"
export PRESENCE_DIR="vessel_mask"
export PROBABILITY_DIR=""
export RADIAL_DIR="k_stellate_maps"
export RELEVANCE_DIR=""
export TEMPLATE_DIR="template_maps"
export WIDTH_DIR="width_maps"

# % Predictor parameters
export PREDICTOR_NAME="predictor"
export PREDICTION_TYPE="rf_classification";

# % Boosted predictor arguments
export BOOST_N_LEVELS=100
export BOOST_WEAK_LEARNER="piecewise_constant"
export BOOST_OUTPUT_TYPE=raw
export BOOST_N_BINS=24
export BOOST_SHRINKAGE=0.05

# %forest construction arguments
export NUM_TREES=200
export SPLIT_CRITERION_C="gdi"
export VAR_CRITERION_C="mabs"
export SPLIT_CRITERION_R="ssq"
export VAR_CRITERION_R="ssq"
export SPLIT_MIN=100
export END_CUT_MIN=25
export DO_UBOUND=0
export DO_CIRCULAR=""
export W_PRIOR=0
export IMPURE_THRESH="1e-4"
export D=""
export OVERWRITE=0

# % Image feature/decomposition parameters
export NUM_LEVELS=5
export RGB_CHANNEL="rgb"
export NORMALISE=0
export WIN_SIZE=3
export PCA_FILENAME=""
export DO_MAX=0
export ROTATE=0
export DECOMP_TYPE="dt"
        
# %DT-CWT representation arguments
export FEATURE_SHAPE="rect"
export FEATURE_TYPE="conj"

# % Gaussian derivative arguments
export SIGMA_RANGE="[1 2 4 8 16]"

# %Addition linop representation arguments
export NUM_ANGLES=8

# %Addition monogenic representation arguments
export MIN_WAVELENGTH=4
export ONF=0.65

# % Raw pixel arguments
export SUBTRACT_MEAN="true"
        
# %General arguments for sampling data
export SAMPLING_METHOD="generate_training_data"
export NUM_SAMPLES=200000
export PTS_PER_IMAGE=500
export FG_RATIO=0.5
export BG_RATIO=1
export MAKE_RESAMPLING_MAPS=0.5
export MAX_N_IMAGES=""
export SHRINK_FOV="false"
export IMAGE_TYPE="real"

# %Arguments for specifying output of data
export OUTPUT_TYPE="width"
        
# %Arguments for creating and sampling real data
export ABNORMAL_DATA="2004_screening/abnormals"
export NORMAL_DATA="2004_screening/normals"
export SAMPLING_METHOD_MAMMO="sample_mammo_training_data"
export SAMPLING_METHOD_RADIAL="load_mass_training_data"
export VIEW=""
export FLIP=0
export NUM_FOLDS=10
export FOLD_ID=1
export DO_TEMPLATE=1
export DO_SCALE=1


# %Arguments for generating synthetic line data
export BG_TYPE="flat"
export BG_SIZE="[64 64]"
export BAR_TYPE="ellipse"
export WIDTH_RANGE="[1 8]"
export CONTRAST_RANGE="[1 8]"
export DECAY_RATE=0
export BG_STEM="bg"
export NUM_BGS=1
export BG_ZEROS=5
export BG_DIR="synthetic_lines/increasing_noise_exp/bg"
export BG_FMT="mat"
export NOISE_TYPE=""
export NOISE_PARAMS=0

# %Arguments for saving/resampling data
export SAVE_TRAINING_DATA="false"
export TRAINING_DATA_DIR="saved_training_data"
export TRAINING_DATA=""

# %classify_image arguments
export MAKE_SAMPLED_MAPS=1
export USE_SAMPLED_MAPS=0
export NUM_TREES_C=""
export MAX_SIZE=128
export USE_PROBS=0
export CLASSIFY_NORMALS=1
export CLASSIFY_ABNORMALS=1
export CLASSIFY_CONTRA=1

# %radial map arguments
export ANGULAR_RESOLUTION=24
export ANGULAR_BANDS=1
export DISTANCE_RANGE="[16 32 64 128 256]"
export RESIZE=""
export R_MIN=4
export R_KARSSEMEIJER=2
export SPACING=2
export OFFSET_X=0
export OFFSET_Y=0

#%Obsolete parameters (keep for back compatibility)
export SPLIT_CRITERION="dabs"
export VAR_CRITERION="mabs"
export DETECTION_TYPE="detection"
export FOREST_DIR="line_detection_rfs"
export MAMMO_RF_DIR="mammo_rfs"
export SAMPLING_METHOD_D="generate_dt_training_data"
export SAMPLING_METHOD_L="generate_linop_training_data"
export SAMPLING_METHOD_M="generate_monogenic_training_data"
export SAMPLING_METHOD_P="generate_pixel_training_data"
export FOREST_TYPE="isbe"