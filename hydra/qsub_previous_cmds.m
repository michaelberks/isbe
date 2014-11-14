%{

NUM_JOBS=10 FOREST_JOB="'191630'" TEST_IMAGE_DIR="'synthetic_lines/lines512'" FOREST_DIR="line_detection_rfs" qsub -t 1-10 -V matlab_code/trunk/hydra/classify_image_set.sh
NUM_JOBS=10 FOREST_JOB="'191933'" TEST_IMAGE_DIR="'synthetic_lines/lines512'" FOREST_DIR="line_orientation_rfs" FOREST_TYPE="regression" qsub -t 1-10 -V matlab_code/trunk/hydra/classify_image_set.sh
NUM_JOBS=10 FOREST_JOB="'191934'" TEST_IMAGE_DIR="'synthetic_lines/lines512'" FOREST_DIR="line_orientation_rfs" FOREST_TYPE="regression" qsub -t 1-10 -V matlab_code/trunk/hydra/classify_image_set.sh

NUM_TREES=20 NUM_BGS=1000 DETECTION_TYPE="orientation" DECOMP_TYPE="mono" WIN_SIZE=1 NUM_LEVELS=4 qsub -N rf_line -t 1:10 -V matlab_code/trunk/hydra/build_rf_line_detector.sh
NUM_TREES=20 NUM_BGS=1000 DETECTION_TYPE="orientation" DECOMP_TYPE="linop" WIN_SIZE=1 qsub -N rf_line -t 1:10 -V matlab_code/trunk/hydra/build_rf_line_detector.sh

NUM_TREES=20 NUM_BGS=1000 DECOMP_TYPE="mono" WIN_SIZE=1 NUM_LEVELS=4 qsub -N rf_line -t 1:10 -V matlab_code/trunk/hydra/build_rf_line_detector.sh
NUM_TREES=20 NUM_BGS=1000 DECOMP_TYPE="linop" WIN_SIZE=1 qsub -N rf_line -t 1:10 -V matlab_code/trunk/hydra/build_rf_line_detector.sh

FOREST_JOB="'191958'" FOREST_DIR="line_orientation_rfs" qsub -V matlab_code/trunk/hydra/combine_hydra_rfs.sh
FOREST_JOB="'191959'" FOREST_DIR="line_orientation_rfs" qsub -V matlab_code/trunk/hydra/combine_hydra_rfs.sh

FOREST_JOB="'191960'" qsub -V matlab_code/trunk/hydra/combine_hydra_rfs.sh
FOREST_JOB="'191961'" qsub -V matlab_code/trunk/hydra/combine_hydra_rfs.sh

NUM_JOBS=10 FOREST_JOB="'191958'" TEST_IMAGE_DIR="'synthetic_lines/lines512'" FOREST_DIR="line_orientation_rfs" FOREST_TYPE="regression" qsub -t 1-10 -V matlab_code/trunk/hydra/classify_image_set.sh
NUM_JOBS=10 FOREST_JOB="'191959'" TEST_IMAGE_DIR="'synthetic_lines/lines512'" FOREST_DIR="line_orientation_rfs" FOREST_TYPE="regression" qsub -t 1-10 -V matlab_code/trunk/hydra/classify_image_set.sh

NUM_JOBS=10 FOREST_JOB="'191960'" TEST_IMAGE_DIR="'synthetic_lines/lines512'" FOREST_DIR="line_detection_rfs" qsub -t 1-10 -V matlab_code/trunk/hydra/classify_image_set.sh
NUM_JOBS=10 FOREST_JOB="'191961'" TEST_IMAGE_DIR="'synthetic_lines/lines512'" FOREST_DIR="line_detection_rfs" qsub -t 1-10 -V matlab_code/trunk/hydra/classify_image_set.sh

VIEW="CC" FOLD_ID=2 NUM_TREES=20 WIN_SIZE=1 FEATURE_TYPE="conj" USE_NAG=1 qsub -N mammo_c -V -t 1-10 matlab_code/trunk/hydra/build_rf_mammo_class.sh

DATA_TYPE="'2004_screening/normals'" qsub -V -t 1-20 matlab_code/trunk/hydra/compute_rad_maps_batch.sh

FOREST_JOB="'lines_01'" FOREST_DIR="mammo_rfs" qsub -V matlab_code/trunk/hydra/combine_hydra_rfs.sh

NUM_TREES=20 NUM_BGS=1000 WIN_SIZE=3 FEATURE_TYPE="all" qsub -t 1:10 -V matlab_code/trunk/hydra/build_rf_line_detector.sh

NUM_JOBS=8 VIEW="CC" CLASSIFY_NORMALS=0 CLASSIFY_ABNORMALS=1 USE_NAG=1 qsub -t 1-8 -V matlab_code/trunk/hydra/classify_mammo_fold.sh
NUM_JOBS=8 VIEW="CC" CLASSIFY_NORMALS=1 CLASSIFY_ABNORMALS=0 USE_NAG=1 qsub -t 1-8 -V matlab_code/trunk/hydra/classify_mammo_fold.sh

NUM_TREES=20 NUM_BGS=1000 DETECTION_TYPE="orientation" DECOMP_TYPE="linop" WIN_SIZE=1 NUM_ANGLES=6 qsub -N rf_line -t 1:10 -V matlab_code/trunk/hydra/build_rf_line_detector.sh
NUM_TREES=20 NUM_BGS=1000 DECOMP_TYPE="pixel" WIN_SIZE=15 qsub -N rf_line -t 1:10 -V matlab_code/trunk/hydra/build_rf_line_detector.sh

CUSTOM_ID="roi" ABNORMAL_DIR="2004_screening/contralateral_roi/abnormals/" NORMAL_DIR="2004_screening/contralateral_roi/normals/" ABNORMAL_MASK_DIR="2004_screening/contralateral_roi/abnormals/" NORMAL_MASK_DIR="2004_screening/contralateral_roi/normals/" FOLD_ID=1 NUM_TREES=20 WIN_SIZE=1 FEATURE_TYPE="conj" USE_NAG=1 qsub -N mammo_c -V -t 1-10 matlab_code/trunk/hydra/build_rf_mammo_class.sh
NUM_JOBS=15 CUSTOM_ID="roi" ABNORMAL_DIR="2004_screening/contralateral_roi/abnormals/" NORMAL_DIR="2004_screening/contralateral_roi/normals/" ABNORMAL_MASK_DIR="2004_screening/contralateral_roi/abnormals/" NORMAL_MASK_DIR="2004_screening/contralateral_roi/normals/" CLASSIFY_NORMALS=1 CLASSIFY_ABNORMALS=1 FOLD_ID=1 qsub -t 1-15 -V matlab_code/trunk/hydra/classify_mammo_fold.sh

DATA_TYPE="'2004_screening/abnormals'" RELEVANCE_DIR="relevance_maps" RADIAL_DIR="weighted_radial_maps" DISTANCE_RANGE="[64 128 256 inf]" qsub -V -t 1-6 matlab_code/trunk/hydra/compute_rad_maps_batch.sh

DATA_TYPE="'2004_screening/abnormals'" RELEVANCE_DIR="relevance_maps" RADIAL_DIR="weighted_radial_maps" DISTANCE_RANGE="[64 128 256 inf]" VIEW="LCC" qsub -V -t 7-20 matlab_code/trunk/hydra/compute_rad_maps_batch.sh
DATA_TYPE="'2004_screening/abnormals'" RELEVANCE_DIR="relevance_maps" RADIAL_DIR="weighted_radial_maps" DISTANCE_RANGE="[64 128 256 inf]" VIEW="LML" qsub -V -t 7-20 matlab_code/trunk/hydra/compute_rad_maps_batch.sh
DATA_TYPE="'2004_screening/abnormals'" RELEVANCE_DIR="relevance_maps" RADIAL_DIR="weighted_radial_maps" DISTANCE_RANGE="[64 128 256 inf]" VIEW="R" qsub -V -t 1-20 matlab_code/trunk/hydra/compute_rad_maps_batch.sh



NUM_JOBS=10 FOREST_JOB="'191905'" TEST_IMAGE_DIR="'mammograms/2004_screening/abnormals/mat'" MASK_DIR="masks/2004_screening/abnormals" qsub -t 1-10 -V matlab_code/trunk/hydra/classify_image_set.sh

NUM_JOBS=94 FOREST_JOB="'191905'" TEST_IMAGE_DIR="'mammograms/2004_screening/normals'" MASK_DIR="masks/2004_screening/normals" FOREST_DIR="line_detection_rfs" qsub -t 21-25:2 -V matlab_code/trunk/hydra/classify_image_set.sh
NUM_JOBS=94 FOREST_JOB="'191905'" TEST_IMAGE_DIR="'mammograms/2004_screening/normals'" MASK_DIR="masks/2004_screening/normals" FOREST_DIR="line_detection_rfs" qsub -t 31,39 -V matlab_code/trunk/hydra/classify_image_set.sh
NUM_JOBS=188 FOREST_JOB="'191905'" TEST_IMAGE_DIR="'mammograms/2004_screening/normals'" MASK_DIR="masks/2004_screening/normals" FOREST_DIR="line_detection_rfs" qsub -t 82,85 -V matlab_code/trunk/hydra/classify_image_set.sh

NUM_JOBS=94 FOREST_JOB="'191934'" TEST_IMAGE_DIR="'mammograms/2004_screening/normals'" MASK_DIR="masks/2004_screening/normals" FOREST_DIR="line_orientation_rfs" FOREST_TYPE="regression" qsub -t 21-25:2 -V matlab_code/trunk/hydra/classify_image_set.sh
NUM_JOBS=94 FOREST_JOB="'191934'" TEST_IMAGE_DIR="'mammograms/2004_screening/normals'" MASK_DIR="masks/2004_screening/normals" FOREST_DIR="line_orientation_rfs" FOREST_TYPE="regression" qsub -t 31,39 -V matlab_code/trunk/hydra/classify_image_set.sh

DATA_TYPE="'2004_screening/abnormals'" LINE_DIR="linop_maps_line" ORI_DIR="linop_maps_ori" RELEVANCE_DIR="relevance_maps" RADIAL_DIR="linop_radial_maps" qsub -V -t 1-20 matlab_code/trunk/hydra/compute_rad_maps_batch.sh

NUM_JOBS=188 FOREST_JOB="'191934'" TEST_IMAGE_DIR="'mammograms/2004_screening/normals'" MASK_DIR="masks/2004_screening/normals" FOREST_DIR="line_orientation_rfs" FOREST_TYPE="regression" qsub -t 82,85 -V matlab_code/trunk/hydra/classify_image_set.sh

NUM_JOBS=10 FOREST_JOB="'191934'" TEST_IMAGE_DIR="'mammograms/2004_screening/abnormals/mat'" MASK_DIR="masks/2004_screening/abnormals" FOREST_DIR="line_orientation_rfs" FOREST_TYPE="regression" qsub -t 1-10 -V matlab_code/trunk/hydra/classify_image_set.sh

CUSTOM_ID="lines" ABNORMAL_DIR="2004_screening/contralateral_roi/abnormals/" NORMAL_DIR="2004_screening/contralateral_roi/normals/" ABNORMAL_MASK_DIR="2004_screening/contralateral_roi/abnormal_lines/" NORMAL_MASK_DIR="2004_screening/contralateral_roi/normal_lines/" FOLD_ID=1 NUM_TREES=20 WIN_SIZE=1 FEATURE_TYPE="conj" USE_NAG=1 qsub -N mammo_c -V -t 1-10 matlab_code/trunk/hydra/build_rf_mammo_class.sh

NUM_JOBS=15 CUSTOM_ID="lines" ABNORMAL_DIR="2004_screening/contralateral_roi/abnormals/" NORMAL_DIR="2004_screening/contralateral_roi/normals/" ABNORMAL_MASK_DIR="2004_screening/contralateral_roi/abnormals/" NORMAL_MASK_DIR="2004_screening/contralateral_roi/normals/" CLASSIFY_NORMALS=1 CLASSIFY_ABNORMALS=1 FOLD_ID=1 qsub -t 1-15 -V matlab_code/trunk/hydra/classify_mammo_fold.sh

NUM_JOBS=94 FOREST_JOB="'lines_01'" TEST_IMAGE_DIR="'mammograms/2004_screening/normals'" MASK_DIR="masks/2004_screening/normals" FOREST_DIR="mammo_rfs" qsub -t 1 -V matlab_code/trunk/hydra/classify_image_set.sh

DATA_TYPE="'2004_screening/abnormals'" RELEVANCE_DIR="relevance_maps" RADIAL_DIR="weighted_radial_maps" DISTANCE_RANGE="[8 16 32 64 128 256 inf]" qsub -V -t 1-20 matlab_code/trunk/hydra/compute_rad_maps_batch.sh
DATA_TYPE="'2004_screening/normals'" RELEVANCE_DIR="relevance_maps" RADIAL_DIR="weighted_radial_maps" DISTANCE_RANGE="[8 16 32 64 128 256 inf]" qsub -V -t 1-20 matlab_code/trunk/hydra/compute_rad_maps_batch.sh

DATA_TYPE="'2004_screening/normals'" RELEVANCE_DIR="relevance_maps" RADIAL_DIR="weighted_radial_maps" DISTANCE_RANGE="[8 16 32 64 128 256 inf]" NUM_JOBS=94  qsub -V -t 21-25:2 matlab_code/trunk/hydra/compute_rad_maps_batch.sh
DATA_TYPE="'2004_screening/normals'" RELEVANCE_DIR="relevance_maps" RADIAL_DIR="weighted_radial_maps" DISTANCE_RANGE="[8 16 32 64 128 256 inf]" NUM_JOBS=94  qsub -V -t 31,39 matlab_code/trunk/hydra/compute_rad_maps_batch.sh
DATA_TYPE="'2004_screening/normals'" RELEVANCE_DIR="relevance_maps" RADIAL_DIR="weighted_radial_maps" DISTANCE_RANGE="[8 16 32 64 128 256 inf]" NUM_JOBS=188  qsub -V -t 82,85 matlab_code/trunk/hydra/compute_rad_maps_batch.sh

NUM_TREES=1 NUM_BGS=1000 WIN_SIZE=1 FEATURE_TYPE="conj" SAMPLING_METHOD="sample_blank_dt_data" qsub -t 1 -V matlab_code/trunk/hydra/build_rf_line_detector.sh

CUSTOM_ID="wr1" MASK_DIR="mass_masks" FOLD_ID=1 NUM_TREES=20 qsub -N mammo_c -V -t 1-10 matlab_code/trunk/hydra/build_rf_radial_class.sh

NUM_TREES=20 NUM_BGS=1000 WIN_SIZE=3 FEATURE_TYPE="conj" NUM_LEVELS=3 qsub -t 1:10 -V matlab_code/trunk/hydra/build_rf_line_detector.sh
NUM_TREES=20 NUM_BGS=1000 WIN_SIZE=3 FEATURE_TYPE="conj" NUM_LEVELS=4 qsub -t 1:10 -V matlab_code/trunk/hydra/build_rf_line_detector.sh

NUM_JOBS=112 FOREST_JOB="'191905'" TEST_IMAGE_DIR="'mammograms/2004_screening_processing/abnormals'" MASK_DIR="masks/2004_screening/abnormals" qsub -t 1-10 -V matlab_code/trunk/hydra/classify_image_set.sh

DATA_TYPE="'2004_screening/abnormals'" MAM_NAMES="mam_names/2004_screening_abnormals"  LINE_DIR="line_maps" ORI_DIR="orientation_maps" RADIAL_DIR="k_stellate_maps_wrf" RELEVANCE_DIR="relevance_maps" RESIZE=0.5 DISTANCE_RANGE="[60 90 120 150 180]" ANGULAR_RESOLUTION="[48 24]" NUM_JOBS=146 qsub -V -t 1-146 matlab_code/trunk/hydra/compute_k_maps_batch.sh
DATA_TYPE="'2004_screening_processed/abnormals'" MAM_NAMES="mam_names/2004_screening_abnormals"  LINE_DIR="line_maps" ORI_DIR="orientation_maps" RADIAL_DIR="k_stellate_maps_rf" RESIZE=0.5 DISTANCE_RANGE="[60 90 120 150 180]" ANGULAR_RESOLUTION="[48 24]" NUM_JOBS=146 qsub -V -t 1-146 matlab_code/trunk/hydra/compute_k_maps_batch.sh
DATA_TYPE="'2004_screening/normals'"  LINE_DIR="line_maps" ORI_DIR="orientation_maps" RADIAL_DIR="k_stellate_maps_wrf" RELEVANCE_DIR="relevance_maps" RESIZE=0.5 DISTANCE_RANGE="[60 90 120 150 180]" ANGULAR_RESOLUTION="[48 24]" NUM_JOBS=188 qsub -V -t 1-188 matlab_code/trunk/hydra/compute_k_maps_batch.sh
DATA_TYPE="'2004_screening/normals'"  LINE_DIR="line_maps" ORI_DIR="orientation_maps" RADIAL_DIR="k_stellate_maps_rf" RESIZE=0.5 DISTANCE_RANGE="[60 90 120 150 180]" ANGULAR_RESOLUTION="[48 24]" NUM_JOBS=188 qsub -V -t 1-188 matlab_code/trunk/hydra/compute_k_maps_batch.sh

DATA_TYPE="'2004_screening_processed/abnormals'" MAM_NAMES="mam_names/2004_screening_abnormals"  LINE_DIR="k_line_maps" ORI_DIR="k_ori_maps" RADIAL_DIR="k_stellate_maps_wg2" RELEVANCE_DIR="relevance_maps" RESIZE=0 DISTANCE_RANGE="[60 90 120 150 180]" ANGULAR_RESOLUTION="[48 24]" NUM_JOBS=146 qsub -V -t 1-146 matlab_code/trunk/hydra/compute_k_maps_batch.sh
DATA_TYPE="'2004_screening_processed/abnormals'" MAM_NAMES="mam_names/2004_screening_abnormals"  LINE_DIR="k_line_maps" ORI_DIR="k_ori_maps" RADIAL_DIR="k_stellate_maps_g2" RESIZE=0 DISTANCE_RANGE="[60 90 120 150 180]" ANGULAR_RESOLUTION="[48 24]" NUM_JOBS=146 qsub -V -t 1-146 matlab_code/trunk/hydra/compute_k_maps_batch.sh
DATA_TYPE="'2004_screening_processed/normals'"  LINE_DIR="k_line_maps" ORI_DIR="k_ori_maps" RADIAL_DIR="k_stellate_maps_wg2" RELEVANCE_DIR="relevance_maps" RESIZE=0 DISTANCE_RANGE="[60 90 120 150 180]" ANGULAR_RESOLUTION="[48 24]" NUM_JOBS=188 qsub -V -t 1-188 matlab_code/trunk/hydra/compute_k_maps_batch.sh
DATA_TYPE="'2004_screening_processed/normals'"  LINE_DIR="k_line_maps" ORI_DIR="k_ori_maps" RADIAL_DIR="k_stellate_maps_g2" RESIZE=0 DISTANCE_RANGE="[60 90 120 150 180]" ANGULAR_RESOLUTION="[48 24]" NUM_JOBS=188 qsub -V -t 1-188 matlab_code/trunk/hydra/compute_k_maps_batch.sh

DATA_TYPE="'2004_screening_processed/abnormals'" LINE_DIR="line_maps_g2" ORI_DIR="orientation_maps_g2" RADIAL_DIR="radial_maps_g2" NUM_JOBS=292 RESIZE=0 qsub -V -t 1-292 matlab_code/trunk/hydra/compute_rad_maps_batch.sh
DATA_TYPE="'2004_screening_processed/normals'" LINE_DIR="line_maps_g2" ORI_DIR="orientation_maps_g2" RADIAL_DIR="radial_maps_g2" NUM_JOBS=188 RESIZE=0 qsub -V -t 1-292 matlab_code/trunk/hydra/compute_rad_maps_batch.sh
DATA_TYPE="'2004_screening_processed/abnormals'" LINE_DIR="line_maps_g2" ORI_DIR="orientation_maps_g2" RELEVANCE_DIR="relevance_maps" RADIAL_DIR="radial_maps_wg2" NUM_JOBS=292 RESIZE=0 qsub -V -t 1-292 matlab_code/trunk/hydra/compute_rad_maps_batch.sh
DATA_TYPE="'2004_screening_processed/normals'" LINE_DIR="line_maps_g2" ORI_DIR="orientation_maps_g2" RELEVANCE_DIR="relevance_maps" RADIAL_DIR="radial_maps_wg2" NUM_JOBS=188 RESIZE=0 qsub -V -t 1-292 matlab_code/trunk/hydra/compute_rad_maps_batch.sh

NUM_TREES=20 NUM_BGS=1000 DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -N rf_line -t 1:10 -V matlab_code/trunk/hydra/build_rf_line_detector.sh
NUM_TREES=20 NUM_BGS=1000 DECOMP_TYPE="g2d" DETECTION_TYPE="orientation" WIN_SIZE=1 qsub -N rf_line -t 1:10 -V matlab_code/trunk/hydra/build_rf_line_detector.sh

NUM_JOBS=20 FOREST_JOB="'233142'" TEST_IMAGE_DIR="'synthetic_lines/lines512/g2d_responses'" FOREST_DIR="line_orientation_rfs" FOREST_TYPE="regression" qsub -t 1-20 -V matlab_

NUM_TREES=10 NUM_BGS=840 BG_DIR="real512/train/" DETECTION_TYPE="orientation" WIN_SIZE=3 NUM_LEVELS=4 FEATURE_TYPE="conj" LINE_TYPE="curve" qsub -t 1:20 -V matlab_code/trunk/hydra/build_rf_line_detector.sh
NUM_TREES=20 NUM_BGS=940 BG_DIR="real512/" DETECTION_TYPE="orientation" WIN_SIZE=3 DECOMP_TYPE="g2d" LINE_TYPE="curve" qsub -t 1:10 -V matlab_code/trunk/hydra/build_rf_line_detector.sh
NUM_TREES=20 NUM_BGS=940 BG_DIR="smooth512/train/" DETECTION_TYPE="orientation" WIN_SIZE=3 DECOMP_TYPE="g2d" LINE_TYPE="curve" qsub -t 1:10 -V matlab_code/trunk/hydra/build_rf_line_detector.sh

NUM_JOBS=292 FOREST_JOB="'306079'" TEST_IMAGE_DIR="'mammograms/2004_screening/abnormals'" MASK_DIR="masks/2004_screening/abnormals" FOREST_DIR="line_orientation_rfs" qsub -t 1-292 -V matlab_code/trunk/hydra/classify_image_set.sh
NUM_JOBS=188 FOREST_JOB="'306079'" TEST_IMAGE_DIR="'mammograms/2004_screening/normals'" MASK_DIR="masks/2004_screening/normals" FOREST_DIR="line_orientation_rfs" qsub -t 1-188 -V matlab_code/trunk/hydra/classify_image_set.sh

%}
NUM_JOBS=292 DATA_TYPE="'2004_screening_processed/abnormals'" METHOD_TYPE="'rf'" MASK_DIR="masks" PECTORAL_DIR="masks_pectoral" qsub -V -t 1 matlab_code/trunk/hydra/compute_k_maps_batch.sh