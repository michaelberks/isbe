% Qsub commands for the CSF

%Compute orientation and line maps for analytic Gaussian methods
NUM_JOBS=47 IMAGE_DIR="normals" qsub -V -l twoday -t 1-47 matlab_code/trunk/hydra/cuc/k_maps_script.sh
NUM_JOBS=73 IMAGE_DIR="abnormals" qsub -V -l twoday -t 1-73 matlab_code/trunk/hydra/cuc/k_maps_script.sh

%Build RF regressors for the Gaussian methods
NUM_TREES=10 NUM_BGS=840 BG_DIR="real512/train/" DETECTION_TYPE="orientation" DECOMP_TYPE="g2d" WIN_SIZE=3 SIGMA_RANGE="[1.1 1.9 3.2]" DO_UBOUND=0 END_CUT_MIN=0 WIDTH_RANGE="[2 16]" LINE_TYPE="sin" SPLIT_CRITERION="ssq" VAR_CRITERION="ssq" qsub -l twoday -t 1-20 -V matlab_code/trunk/hydra/cuc/build_rf_line_detector.sh
NUM_TREES=10 NUM_BGS=840 BG_DIR="real512/train/" DETECTION_TYPE="orientation" DECOMP_TYPE="g2d" WIN_SIZE=3 SIGMA_RANGE="[1.1]" DO_UBOUND=0 END_CUT_MIN=0 WIDTH_RANGE="[2 16]" LINE_TYPE="sin" SPLIT_CRITERION="ssq" VAR_CRITERION="ssq" qsub -l twoday -t 1-20 -V matlab_code/trunk/hydra/cuc/build_rf_line_detector.sh
NUM_TREES=10 NUM_BGS=840 BG_DIR="real512/train/" DETECTION_TYPE="orientation" DECOMP_TYPE="g2d" WIN_SIZE=3 SIGMA_RANGE="[1.9]" DO_UBOUND=0 END_CUT_MIN=0 WIDTH_RANGE="[2 16]" LINE_TYPE="sin" SPLIT_CRITERION="ssq" VAR_CRITERION="ssq" qsub -l twoday -t 1-20 -V matlab_code/trunk/hydra/cuc/build_rf_line_detector.sh
NUM_TREES=10 NUM_BGS=840 BG_DIR="real512/train/" DETECTION_TYPE="orientation" DECOMP_TYPE="g2d" WIN_SIZE=3 SIGMA_RANGE="[3.2]" DO_UBOUND=0 END_CUT_MIN=0 WIDTH_RANGE="[2 16]" LINE_TYPE="sin" SPLIT_CRITERION="ssq" VAR_CRITERION="ssq" qsub -l twoday -t 1-20 -V matlab_code/trunk/hydra/cuc/build_rf_line_detector.sh

FOREST_JOB="'63707'" FOREST_DIR="line_orientation_rfs" qsub -V -hold_jid 63707 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
FOREST_JOB="'64260'" FOREST_DIR="line_orientation_rfs" qsub -V -hold_jid 64260 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
FOREST_JOB="'64261'" FOREST_DIR="line_orientation_rfs" qsub -V -hold_jid 64261 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
FOREST_JOB="'3XXXX'" FOREST_DIR="line_orientation_rfs" qsub -V -hold_jid 3XXXX -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh

%Run RF regressors over mammograms
NUM_JOBS=73 FOREST_JOB="'63707'" TEST_IMAGE_DIR="'mammograms/2004_screening_processed/abnormals'" MASK_DIR="masks/2004_screening_processed/abnormals" FOREST_DIR="line_orientation_rfs" qsub -hold_jid 63708 -t 1-73 -V -l twoday matlab_code/trunk/hydra/cuc/classify_image_set.sh
NUM_JOBS=47 FOREST_JOB="'63707'" TEST_IMAGE_DIR="'mammograms/2004_screening_processed/normals'" MASK_DIR="masks/2004_screening_processed/normals" FOREST_DIR="line_orientation_rfs" qsub -hold_jid 63708 -t 1-47 -V -l twoday matlab_code/trunk/hydra/cuc/classify_image_set.sh

NUM_JOBS=73 FOREST_JOB="'64260'" TEST_IMAGE_DIR="'mammograms/2004_screening_processed/abnormals'" MASK_DIR="masks/2004_screening_processed/abnormals" FOREST_DIR="line_orientation_rfs" qsub -hold_jid 64262 -t 1-73 -V -l twoday matlab_code/trunk/hydra/cuc/classify_image_set.sh
NUM_JOBS=47 FOREST_JOB="'64260'" TEST_IMAGE_DIR="'mammograms/2004_screening_processed/normals'" MASK_DIR="masks/2004_screening_processed/normals" FOREST_DIR="line_orientation_rfs" qsub -hold_jid 64262 -t 1-47 -V -l twoday matlab_code/trunk/hydra/cuc/classify_image_set.sh

NUM_JOBS=73 FOREST_JOB="'64261'" TEST_IMAGE_DIR="'mammograms/2004_screening_processed/abnormals'" MASK_DIR="masks/2004_screening_processed/abnormals" FOREST_DIR="line_orientation_rfs" qsub -hold_jid 64263 -t 1-73 -V -l twoday matlab_code/trunk/hydra/cuc/classify_image_set.sh
NUM_JOBS=47 FOREST_JOB="'64261'" TEST_IMAGE_DIR="'mammograms/2004_screening_processed/normals'" MASK_DIR="masks/2004_screening_processed/normals" FOREST_DIR="line_orientation_rfs" qsub -hold_jid 64263 -t 1-47 -V -l twoday matlab_code/trunk/hydra/cuc/classify_image_set.sh

NUM_JOBS=73 FOREST_JOB="'3XXXX'" TEST_IMAGE_DIR="'mammograms/2004_screening_processed/abnormals'" MASK_DIR="masks/2004_screening_processed/abnormals" FOREST_DIR="line_orientation_rfs" qsub -t 1-73 -V -l twoday matlab_code/trunk/hydra/cuc/classify_image_set.sh
NUM_JOBS=47 FOREST_JOB="'3XXXX'" TEST_IMAGE_DIR="'mammograms/2004_screening_processed/normals'" MASK_DIR="masks/2004_screening_processed/normals" FOREST_DIR="line_orientation_rfs" qsub -t 1-47 -V -l twoday matlab_code/trunk/hydra/cuc/classify_image_set.sh

%Move the results from the mammograms/results folders to the orientation
%maps folders
mkdir scratch/asym/data/orientation_maps/g2d_rf_all/2004_screening_processed/normals
mkdir scratch/asym/data/orientation_maps/g2d_rf_all/2004_screening_processed/abnormals

mkdir scratch/asym/data/orientation_maps/g2d_rf_1
mkdir scratch/asym/data/orientation_maps/g2d_rf_1/2004_screening_processed
mkdir scratch/asym/data/orientation_maps/g2d_rf_1/2004_screening_processed/normals
mkdir scratch/asym/data/orientation_maps/g2d_rf_1/2004_screening_processed/abnormals

mkdir scratch/asym/data/orientation_maps/g2d_rf_2
mkdir scratch/asym/data/orientation_maps/g2d_rf_2/2004_screening_processed
mkdir scratch/asym/data/orientation_maps/g2d_rf_2/2004_screening_processed/normals
mkdir scratch/asym/data/orientation_maps/g2d_rf_2/2004_screening_processed/abnormals

mkdir scratch/asym/data/orientation_maps/g2d_rf_3/2004_screening_processed/normals
mkdir scratch/asym/data/orientation_maps/g2d_rf_3/2004_screening_processed/abnormals

cp scratch/asym/data/mammograms/2004_screening_processed/normals/results/63707/* scratch/asym/data/orientation_maps/g2d_rf_all/2004_screening_processed/normals
cp scratch/asym/data/mammograms/2004_screening_processed/abnormals/results/63707/* scratch/asym/data/orientation_maps/g2d_rf_all/2004_screening_processed/abnormals

cp scratch/asym/data/mammograms/2004_screening_processed/normals/results/64260/* scratch/asym/data/orientation_maps/g2d_rf_1/2004_screening_processed/normals
cp scratch/asym/data/mammograms/2004_screening_processed/abnormals/results/64260/* scratch/asym/data/orientation_maps/g2d_rf_1/2004_screening_processed/abnormals

cp scratch/asym/data/mammograms/2004_screening_processed/normals/results/64261/* scratch/asym/data/orientation_maps/g2d_rf_2/2004_screening_processed/normals
cp scratch/asym/data/mammograms/2004_screening_processed/abnormals/results/64261/* scratch/asym/data/orientation_maps/g2d_rf_2/2004_screening_processed/abnormals

cp scratch/asym/data/mammograms/2004_screening_processed/normals/results/3XXXX/* scratch/asym/data/orientation_maps/g2d_rf_3/2004_screening_processed/normals
cp scratch/asym/data/mammograms/2004_screening_processed/abnormals/results/3XXXX/* scratch/asym/data/orientation_maps/g2d_rf_3/2004_screening_processed/abnormals

%Build k-maps for all the methods...
NUM_JOBS=73 DATA_TYPE="'2004_screening_processed/abnormals'" METHOD_TYPE="'g2d_all'" MASK_DIR="masks" PECTORAL_DIR="masks_pectoral" LINE_DIR="line_maps" ALL_PIXELS_K=0 SPACING=4 qsub -l twoday -V -t 1-73 matlab_code/trunk/hydra/cuc/compute_k_maps_batch.sh
NUM_JOBS=47 DATA_TYPE="'2004_screening_processed/normals'" METHOD_TYPE="'g2d_all'" MASK_DIR="masks" PECTORAL_DIR="masks_pectoral" LINE_DIR="line_maps" ALL_PIXELS_K=0 SPACING=4 qsub -l twoday -V -t 1-47 matlab_code/trunk/hydra/cuc/compute_k_maps_batch.sh

NUM_JOBS=73 DATA_TYPE="'2004_screening_processed/abnormals'" METHOD_TYPE="'g2d_1'" MASK_DIR="masks" PECTORAL_DIR="masks_pectoral" LINE_DIR="line_maps" ALL_PIXELS_K=0 SPACING=4 qsub -l twoday -V -t 1-73 matlab_code/trunk/hydra/cuc/compute_k_maps_batch.sh
NUM_JOBS=47 DATA_TYPE="'2004_screening_processed/normals'" METHOD_TYPE="'g2d_1'" MASK_DIR="masks" PECTORAL_DIR="masks_pectoral" LINE_DIR="line_maps" ALL_PIXELS_K=0 SPACING=4 qsub -l twoday -V -t 1-47 matlab_code/trunk/hydra/cuc/compute_k_maps_batch.sh

NUM_JOBS=73 DATA_TYPE="'2004_screening_processed/abnormals'" METHOD_TYPE="'g2d_2'" MASK_DIR="masks" PECTORAL_DIR="masks_pectoral" LINE_DIR="line_maps" ALL_PIXELS_K=0 SPACING=4 qsub -l twoday -V -t 1-73 matlab_code/trunk/hydra/cuc/compute_k_maps_batch.sh
NUM_JOBS=47 DATA_TYPE="'2004_screening_processed/normals'" METHOD_TYPE="'g2d_2'" MASK_DIR="masks" PECTORAL_DIR="masks_pectoral" LINE_DIR="line_maps" ALL_PIXELS_K=0 SPACING=4 qsub -l twoday -V -t 1-47 matlab_code/trunk/hydra/cuc/compute_k_maps_batch.sh

NUM_JOBS=73 DATA_TYPE="'2004_screening_processed/abnormals'" METHOD_TYPE="'g2d_3'" MASK_DIR="masks" PECTORAL_DIR="masks_pectoral" LINE_DIR="line_maps" ALL_PIXELS_K=0 SPACING=4 qsub -l twoday -V -t 1-73 matlab_code/trunk/hydra/cuc/compute_k_maps_batch.sh
NUM_JOBS=47 DATA_TYPE="'2004_screening_processed/normals'" METHOD_TYPE="'g2d_3'" MASK_DIR="masks" PECTORAL_DIR="masks_pectoral" LINE_DIR="line_maps" ALL_PIXELS_K=0 SPACING=4 qsub -l twoday -V -t 1-47 matlab_code/trunk/hydra/cuc/compute_k_maps_batch.sh

NUM_JOBS=73 DATA_TYPE="'2004_screening_processed/abnormals'" METHOD_TYPE="'g2d_rf_all'" MASK_DIR="masks" PECTORAL_DIR="masks_pectoral" LINE_DIR="" ALL_PIXELS_K=0 SPACING=4 qsub -l twoday -V -t 1-73 matlab_code/trunk/hydra/cuc/compute_k_maps_batch.sh
NUM_JOBS=47 DATA_TYPE="'2004_screening_processed/normals'" METHOD_TYPE="'g2d_rf_all'" MASK_DIR="masks" PECTORAL_DIR="masks_pectoral" LINE_DIR="" ALL_PIXELS_K=0 SPACING=4 qsub -l twoday -V -t 1-47 matlab_code/trunk/hydra/cuc/compute_k_maps_batch.sh

NUM_JOBS=73 DATA_TYPE="'2004_screening_processed/abnormals'" METHOD_TYPE="'g2d_rf_1'" MASK_DIR="masks" PECTORAL_DIR="masks_pectoral" LINE_DIR="" ALL_PIXELS_K=0 SPACING=4 qsub -l twoday -V -t 1-73 matlab_code/trunk/hydra/cuc/compute_k_maps_batch.sh
NUM_JOBS=47 DATA_TYPE="'2004_screening_processed/normals'" METHOD_TYPE="'g2d_rf_1'" MASK_DIR="masks" PECTORAL_DIR="masks_pectoral" LINE_DIR="" ALL_PIXELS_K=0 SPACING=4 qsub -l twoday -V -t 1-47 matlab_code/trunk/hydra/cuc/compute_k_maps_batch.sh

NUM_JOBS=73 DATA_TYPE="'2004_screening_processed/abnormals'" METHOD_TYPE="'g2d_rf_2'" MASK_DIR="masks" PECTORAL_DIR="masks_pectoral" LINE_DIR="" ALL_PIXELS_K=0 SPACING=4 qsub -l twoday -V -t 1-73 matlab_code/trunk/hydra/cuc/compute_k_maps_batch.sh
NUM_JOBS=47 DATA_TYPE="'2004_screening_processed/normals'" METHOD_TYPE="'g2d_rf_2'" MASK_DIR="masks" PECTORAL_DIR="masks_pectoral" LINE_DIR="" ALL_PIXELS_K=0 SPACING=4 qsub -l twoday -V -t 1-47 matlab_code/trunk/hydra/cuc/compute_k_maps_batch.sh


