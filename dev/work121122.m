DATA_NAME="STARE" MAKE_DATA=1 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=1 qsub -V -t 21-30 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh

DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=1 NUM_TREES=100 qsub -V -hold_jid 57819 -t 21-30 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=1 NUM_TREES=100 qsub -V -hold_jid 57819 -t 21-30 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=1 NUM_TREES=100 qsub -V -hold_jid 57819 -t 21-30 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=1 NUM_TREES=100 qsub -V -hold_jid 57819 -t 21-30 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=2 NUM_TREES=100 qsub -V -hold_jid 57819 -t 21-30 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=2 NUM_TREES=100 qsub -V -hold_jid 57819 -t 21-30 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=3 NUM_TREES=100 qsub -V -hold_jid 57819 -t 21-30 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=3 NUM_TREES=100 qsub -V -hold_jid 57819 -t 21-30 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=4 NUM_TREES=100 qsub -V -hold_jid 57819 -t 21-30 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh

DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=2 NUM_TREES=100 qsub -V -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=2 NUM_TREES=100 qsub -V -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=3 NUM_TREES=100 qsub -V -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=3 NUM_TREES=100 qsub -V -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=4 NUM_TREES=100 qsub -V -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh

MAKE_DATA=1 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=1 qsub -V -t 61 -l short matlab_code/trunk/hydra/cuc/compare_fibre_experiment.csf
%%
comparing_fibre_experiment_csf(1, ...
    'num_pts',      unixenv('NUM_SAMPLES', 2000), ...
    'num_trees',    unixenv('NUM_TREES', 2), ...
    'make_data',    unixenv('MAKE_DATA',1), ...
    'do_orientation', unixenv('DO_ORIENTATION',1), ...
    'do_detection', unixenv('DO_DETECTION',1), ...
    'do_tests', unixenv('DO_TESTS',1), ...
    'win_sizes', unixenv('WIN_SIZES',1) ...
);
%%
fibre_list = dir([asymmetryroot 'data\fibre\training\fibre_masks\*.mat']);
for ii = 1:200
    load([asymmetryroot 'data\fibre\training\fibre_masks\' fibre_list(ii).name], 'fibre_mask');
    num_pts(ii,1) = sum(fibre_mask(:));
end
    
%%
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=1 NUM_TREES=100 qsub -V -hold_jid 58075 -t 101-110 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=1 NUM_TREES=100 qsub -V -hold_jid 58075 -t 101-110 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=1 NUM_TREES=100 qsub -V -hold_jid 58075 -t 101-110 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=1 NUM_TREES=100 qsub -V -hold_jid 58075 -t 101-110 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=7 NUM_TREES=100 qsub -V -hold_jid 58075 -t 101-110 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=7 NUM_TREES=100 qsub -V -hold_jid 58075 -t 101-110 -l twoday matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
%%
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=3 NUM_TREES=100 qsub -V -t 51-60 -l twoday matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=3 NUM_TREES=100 qsub -V -t 51-60 -l twoday matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=2 NUM_TREES=100 qsub -V -t 1-10 -l highmem matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=2 NUM_TREES=100 qsub -V -t 1-10 -l highmem matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=2 NUM_TREES=100 qsub -V -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=2 NUM_TREES=100 qsub -V -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=2 NUM_TREES=100 qsub -V -t 21-30 -l highmem matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=2 NUM_TREES=100 qsub -V -t 21-30 -l highmem matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=2 NUM_TREES=100 qsub -V -t 21-30 -l twoday matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
DATA_NAME="STARE" MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=2 NUM_TREES=100 qsub -V -t 21-30 -l twoday matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
DATA_NAME="STARE" MAKE_DATA=1 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=3 NUM_TREES=100 qsub -V -t 101-110 -l twoday matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh

%%
MAKE_DATA=1 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=1 qsub -V -t 1-60 -l twoday matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=1 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=1 qsub -V -t 71-100 -l twoday matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh

MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=1 NUM_TREES=100 qsub -V -t 1-100 matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=1 NUM_TREES=100 qsub -V -t 1-100 matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=1 NUM_TREES=100 qsub -V -t 1-100 matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=1 NUM_TREES=100 qsub -V -t 1-100 matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh

MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=2 NUM_TREES=100 qsub -V -t 1-100 -l highmem matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=2 NUM_TREES=100 qsub -V -t 1-100 -l highmem matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=3 NUM_TREES=100 qsub -V -t 1-100 -l twoday matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=3 NUM_TREES=100 qsub -V -t 1-100 -l twoday matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=4 NUM_TREES=100 qsub -V -t 1-100 -l twoday matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=5 NUM_TREES=100 qsub -V -t 1-100 -l twoday matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=5 NUM_TREES=100 qsub -V -t 1-100 -l twoday matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=6 NUM_TREES=100 qsub -V -t 1-100 -l vhighmem matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=6 NUM_TREES=100 qsub -V -t 1-100 -l vhighmem matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=7 NUM_TREES=100 qsub -V -t 1-100 -l twoday matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=7 NUM_TREES=100 qsub -V -t 1-100 -l twoday matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
%%
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="[3]" LEVELS="[1:16]" NUM_ANGLES=18 EXP_NAME="all" qsub -V -t 2-10 -l highmem matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="[3]" LEVELS="[1:16]" NUM_ANGLES=18 EXP_NAME="all" qsub -V -t 2-10 -l highmem matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="[1]" LEVELS="[1:16]" NUM_ANGLES=18 EXP_NAME="all" qsub -V -t 2-10 -l highmem matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="[1]" LEVELS="[1:16]" NUM_ANGLES=18 EXP_NAME="all" qsub -V -t 2-10 -l highmem matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh

MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="[3]" LEVELS="[1 2 4 8 16]" NUM_ANGLES=18 EXP_NAME="angles" qsub -V -t 2-10 -l highmem matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="[3]" LEVELS="[1 2 4 8 16]" NUM_ANGLES=18 EXP_NAME="angles" qsub -V -t 2-10 -l highmem matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="[1]" LEVELS="[1 2 4 8 16]" NUM_ANGLES=18 EXP_NAME="angles" qsub -V -t 2-10 -l highmem matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="[1]" LEVELS="[1 2 4 8 16]" NUM_ANGLES=18 EXP_NAME="angles" qsub -V -t 2-10 -l highmem matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh

MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="[3]" LEVELS="[1:16]" NUM_ANGLES=6 EXP_NAME="scales" qsub -V -t 2-10 -l highmem matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="[3]" LEVELS="[1:16]" NUM_ANGLES=6 EXP_NAME="scales" qsub -V -t 2-10 -l highmem matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="[1]" LEVELS="[1:16]" NUM_ANGLES=6 EXP_NAME="scales" qsub -V -t 2-10 -l highmem matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="[1]" LEVELS="[1:16]" NUM_ANGLES=6 EXP_NAME="scales" qsub -V -t 2-10 -l highmem matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh
%%
comparing_DRIVE_g2d_csf(1, ... % non-strict mode
    'exp_name',     'jim', ...
    'image_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/images/'], ...
    'vessel_mask_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/vessel_masks/'], ...
    'fov_mask_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/fov_masks/'], ...
    'ori_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/orientations/'], ...
    'width_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/width_maps/'], ...
    'selected_images', [],...
    'num_pts',      unixenv('NUM_SAMPLES', 200), ...
    'num_trees',    unixenv('NUM_TREES', 2), ...
    'make_data',    unixenv('MAKE_DATA',1), ...
    'num_angles',   unixenv('NUM_ANGLES', 6),...
    'levels',       unixenv('LEVELS', [1 2 4]),...
    'do_orientation', unixenv('DO_ORIENTATION',0), ...
    'do_detection', unixenv('DO_DETECTION',0), ...
    'do_width', unixenv('DO_WIDTH',0), ...
    'win_sizes', unixenv('WIN_SIZES',[1 3]) ...
);
%%
rm -r scratch/asym/experiments/fibre/comparing_representations/results/1/detection/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/2/detection/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/3/detection/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/4/detection/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/5/detection/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/6/detection/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/7/detection/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/8/detection/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/9/detection/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/10/detection/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/1/orientation/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/2/orientation/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/3/orientation/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/4/orientation/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/5/orientation/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/6/orientation/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/7/orientation/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/8/orientation/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/9/orientation/gabor
rm -r scratch/asym/experiments/fibre/comparing_representations/results/10/orientation/gabor
%%
MAKE_DATA=1 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=1 NUM_TREES=100 qsub -V -t 21-30 matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=1 NUM_TREES=100 qsub -V -l twoday -t 21-30 matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=1 NUM_TREES=100 qsub -V -l twoday -t 21-30 matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=1 NUM_TREES=100 qsub -V -l highmem -t 21-30 matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=1 NUM_TREES=100 qsub -V -l highmem -t 21-30 matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh

MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=2 NUM_TREES=100 qsub -V -t 21-30 -l highmem matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=2 NUM_TREES=100 qsub -V -t 21-30 -l highmem matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=3 NUM_TREES=100 qsub -V -t 21-30 -l twoday matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=3 NUM_TREES=100 qsub -V -t 21-30 -l twoday matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=4 NUM_TREES=100 qsub -V -t 21-30 -l twoday matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=4 NUM_TREES=100 qsub -V -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh

MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=4 NUM_TREES=100 qsub -V -t 81-90 -l twoday matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=6 NUM_TREES=100 qsub -V -t 61-70 -l vhighmem matlab_code/trunk/hydra/cuc/

MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=7 NUM_TREES=100 qsub -V -l twoday -t 11-20 matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=7 NUM_TREES=100 qsub -V -l twoday -t 11-20 matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=2 NUM_TREES=100 qsub -V -l highmem -t 1-10 matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=2 NUM_TREES=100 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=2 NUM_TREES=100 qsub -V -l highmem -t 21-30 matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=2 NUM_TREES=100 qsub -V -l twoday -t 21-30 matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh

MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=1 NUM_TREES=100 qsub -V -l twoday -t 101-110 matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=1 NUM_TREES=100 qsub -V -l twoday -t 101-110 matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=1 NUM_TREES=100 qsub -V -l twoday -t 101-110 matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh

MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=7 NUM_TREES=100 qsub -V -l twoday -t 101-110 matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=7 NUM_TREES=100 qsub -V -l twoday -t 101-110 matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh

MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=5 NUM_TREES=100 qsub -V -l twoday -t 41-50 matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=5 NUM_TREES=100 qsub -V -l twoday -t 41-50 matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=5 NUM_TREES=100 qsub -V -l twoday -t 41-50 matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=5 NUM_TREES=100 qsub -V -l twoday -t 41-50 matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh

MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=8 NUM_TREES=50 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh

MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=5 NUM_TREES=100 qsub -V -l twoday -t 41-50 matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=5 NUM_TREES=100 qsub -V -l twoday -t 41-50 matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=5 NUM_TREES=100 qsub -V -l twoday -t 41-50 matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="1" DO_TESTS=5 NUM_TREES=100 qsub -V -l twoday -t 41-50 matlab_code/trunk/hydra/cuc/compare_fibre_experiment_csf.sh

MAKE_DATA=1 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="3" DO_TESTS=1 NUM_TREES=100 qsub -V -l twoday -t 1 matlab_code/trunk/hydra/cuc/compare_DRIVE_experiment_csf.sh

%%

base_dir = 'C:\isbe\asymmetry_project\experiments\DRIVE\comparing_representations\results\';

d1 = {'orig' 'feature_types' 'levels' 'do_max' 'rotate'};
d2{1} = {'1' '3'};
d2{2} = {'conj' 'all' 'imag' 'real' 'real_imag' 'mag' 'phase' 'real_abs_imag'};
d2{3} = {'1\3' '2\3' '3\3' '4\3' '5\3'};
d2{4} = {'3'};
d2{5} = {'3'};

for repeat = 1:10
    g_dir = [base_dir num2str(repeat) '\orientation\gabor\'];

    for i1 = 1:5 
        for i2 = 1:length(d2{i1})
            
            r_dir = [g_dir d1{i1} '\' d2{i1}{i2}];
            r_list = dir([r_dir '\*.mat']);
            
            for i_r = 1:length(r_list)
                
                if str2num(r_list(i_r).date(1:2)) < 20
                    movefile(...
                        [r_dir '\' r_list(i_r).name],...
                        [r_dir '\old_' r_list(i_r).name]);
                end
            end
        end
    end
end
%%
idx = 1:1e4;
for ii = 0:19
    figure; hist(predicted_lines(idx), 100);
    idx = idx + 1e4;
end
%%
comparing_fibre_experiment_csf(41, ... % non-strict mode
    'image_dir',      [asymmetryroot,'data/fibre/training/images/'], ...
    'fibre_mask_dir',      [asymmetryroot,'data/fibre/training/fibre_masks/'], ...
    'fov_mask_dir',      [asymmetryroot,'data/fibre/training/fov_masks/'], ...
    'ori_dir',      [asymmetryroot,'data/fibre/training/orientation_maps/'], ...
    'num_images', [],...
    'num_pts',      unixenv('NUM_SAMPLES', 10000), ...
    'num_trees',    unixenv('NUM_TREES', 2), ...
    'make_data',    unixenv('MAKE_DATA',1), ...
    'do_orientation', unixenv('DO_ORIENTATION',0), ...
    'do_detection', unixenv('DO_DETECTION',0), ...
    'do_tests', unixenv('DO_TESTS',0), ...
    'win_sizes', unixenv('WIN_SIZES',[1 3]) ...
);
%%
MAKE_DATA=1 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="[3]" LEVELS="[1]" qsub -V -t 1 -l short matlab_code/trunk/hydra/cuc/compare_DRIVE_g2d_csf.sh
%%
for ii = 1:20
    mask2 = imread(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\1st_manual\' zerostr(ii,2) '_manual1.gif']);
    mask1 = imread(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\2nd_manual\' zerostr(ii,2) '_manual2.gif']);
    %mask = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\vessel_masks\' zerostr(ii,2) '_test_v_mask.mat']);
    fov_mask = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\fov_masks\' zerostr(ii,2) '_test_f_mask.mat']);
    
    mask = mask1 > 0;
    mask2 = mask2 > 0;
    
    t_pos(ii) = sum(mask2(mask & fov_mask));
    f_pos(ii) = sum(mask2(~mask & fov_mask));
    pos(ii) = sum(mask(:) & fov_mask(:));
    neg(ii) = sum(~mask(:) & fov_mask(:));
    
    figure; 
    subplot(2,2,1); imgray(mask);
    subplot(2,2,2); imgray(mask2);
    subplot(2,2,3); imgray(mask & ~mask2);
    subplot(2,2,4); imgray(~mask & mask2);
end

display(sum(t_pos) / sum(pos))
display(sum(f_pos) / sum(neg))
%%
mkdir C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\vessel_masks2\
for ii = 1:20
    v_mask = imread(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\2nd_manual\' zerostr(ii,2) '_manual2.gif']); 
    v_mask = v_mask > 0;
    save(['C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\vessel_masks2\' zerostr(ii,2) '_test_v_mask.mat'], 'v_mask');
end
%%
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'24303'" NUM_JOBS=10 IMAGE_ROOT="retinograms/STARE/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
%%
rf_codes( 1,:) = {'9162', 'dt', '3', 'orig'};
rf_codes( 2,:) = {'9160', 'g2', '3', 'orig'};
rf_codes( 3,:) = {'12393', 'g12d', '3', 'orig'};
rf_codes( 4,:) = {'12394', 'mono', '3', 'orig'};
rf_codes( 5,:) = {'24303', 'gabor_a', '3', 'orig'}; %Old version of filters '11480'
rf_codes( 6,:) = {'42381', 'g2di', '3', 'orig'};
rf_codes( 7,:) = {'42382', 'gabor_i', '3', 'orig'};
rf_codes( 8,:) = {'50981', 'gabor_big_scale', '3', 'orig'};
rf_codes( 9,:) = {'51417', 'gabor_big_angle', '3', 'orig'};

rf_codes( 1,:) = {'42813', 'dt', '3', '0.75'};
rf_codes( 2,:) = {'13293', 'g2', '3', '0.75'};
rf_codes( 3,:) = {'13294', 'g12d', '3', '0.75'};
rf_codes( 4,:) = {'13295', 'mono', '3', '0.75'};
rf_codes( 5,:) = {'24308', 'gabor_a', '3', '0.75'}; %Old version of filters '11480'
rf_codes( 6,:) = {'42373', 'g2di', '3', '0.75'};
rf_codes( 7,:) = {'42374', 'gabor_i', '3', '0.75'};
MODEL_ROOT="models/vessel/orientation/rf_regression" MODEL_PATH="'42813'" NUM_JOBS=20 IMAGE_ROOT="retinograms/STARE/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -t 1-20 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/vessel/orientation/rf_regression" MODEL_PATH="'13295'" NUM_JOBS=20 IMAGE_ROOT="retinograms/STARE/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -t 1-20 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/vessel/orientation/rf_regression" MODEL_PATH="'24308'" NUM_JOBS=20 IMAGE_ROOT="retinograms/STARE/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -t 1-20 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/vessel/orientation/rf_regression" MODEL_PATH="'88171'" NUM_JOBS=20 IMAGE_ROOT="retinograms/STARE/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -t 1-20 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
%%
for ii = 1:20
    load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\predictions\detection\rf_classification\24303\' zerostr(ii,2) '_test_pred.mat'])
    figure; imgray(prediction_image);
end
%%
for ii = 1:20
    load(['C:\isbe\asymmetry_project\data\retinograms\STARE\training\predictions\detection\rf_classification\24303\' zerostr(ii,2) '_training_pred.mat'])
    figure; imgray(prediction_image);
end
%%
for ii = 2%1:20
    p1 = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\predictions\detection\rf_classification\9160\' zerostr(ii,2) '_test_pred.mat']);
    p2 = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\predictions\detection\rf_classification\12393\' zerostr(ii,2) '_test_pred.mat']);
    figure; 
    subplot(1,2,1); imgray(p1);
    subplot(1,2,2); imgray(p2);
end
%%
for ii = 1:9
    load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\predictions\detection\rf_classification\' rf_codes{ii,1} '\14_test_pred.mat'])
    figure; imgray(prediction_image);
end
%%
im = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\images\14_test.mat']);      
%f_mask = u_load('C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\fov_masks\14_test_f_mask.mat');
f_mask = false(size(im));
f_mask(200:250, 200:250) = 1;
%%
im = u_load(['C:\isbe\asymmetry_project\data\retinograms\STARE\training\images\01_training.mat']);      
f_mask = false(size(im));
f_mask(200:250, 200:250) = 1;
%%
im = 100 - create_gauss_bar(1, 10, 22.5, 256, 256, 128, 128);
f_mask = false(size(im));
f_mask(97:160, 97:160) = 1;
%%
rf = u_load('C:\isbe\asymmetry_project\data\models\vessel\detection\rf_classification\9162\predictor.mat');
job_args = u_load('C:\isbe\asymmetry_project\data\models\vessel\detection\rf_classification\9162\job_args.mat');
job_args.decomposition_args.decomp_type = {job_args.decomposition_args.decomp_type};
rf.tree_root = 'C:\isbe\asymmetry_project\data\models\vessel\detection\rf_classification\';
[p_pad] = predict_image(...
    'image_in', im,...
    'decomposition_args', job_args.decomposition_args,...
    'predictor', rf, ...
    'prediction_type', 'rf_classification',...
    'output_type', 'detection',...
    'use_probs', 0,...
    'mask', f_mask,...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 128,...
    'incremental_results', 0);
%%
[roc_pts, auc, auc_individual, tp_counts, fp_counts, t_counts, f_counts] = compute_roc_image_set([pred_dir rf_codes{9,1} '\'], label_dir, fov_mask_dir);
%%
DECOMP_TYPE="gabor" WIN_SIZE=3 PROBABILITY_DIR="24303/" MAKE_RESAMPLING_MAPS=0.75 MAX_N_IMAGES=10 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="retinograms/STARE/training" NUM_SAMPLES=100000 MODEL_ROOT="models/vessel" qsub -V -l twoday -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88155'" qsub -V -hold_jid 88155 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88155'" NUM_JOBS=20 IMAGE_ROOT="retinograms/STARE/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 88156 -t 11-20 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh

MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'9162'" NUM_JOBS=20 IMAGE_ROOT="retinograms/STARE/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -t 1-20 -l highmem matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'12394'" NUM_JOBS=20 IMAGE_ROOT="retinograms/STARE/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -t 1-20 -l highmem matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88126'" NUM_JOBS=20 IMAGE_ROOT="retinograms/STARE/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -t 1-20 -l highmem matlab_code/trunk/hydra/cuc/predict_image_set.sh

DECOMP_TYPE="gabor" WIN_SIZE=3 PROBABILITY_DIR="24303/" MAKE_RESAMPLING_MAPS=0.75 MAX_N_IMAGES=10 SHIFT_IMAGES=10 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="retinograms/STARE/training" NUM_SAMPLES=100000 MODEL_ROOT="models/vessel" qsub -V -l twoday -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="dt" WIN_SIZE=3 PROBABILITY_DIR="9162/" MAKE_RESAMPLING_MAPS=0.75 MAX_N_IMAGES=10 SHIFT_IMAGES=10 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="retinograms/STARE/training" NUM_SAMPLES=100000 MODEL_ROOT="models/vessel" qsub -V -l highmem -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="mono" WIN_SIZE=3 PROBABILITY_DIR="12394/" MAKE_RESAMPLING_MAPS=0.75 MAX_N_IMAGES=10 SHIFT_IMAGES=10 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="retinograms/STARE/training" NUM_SAMPLES=100000 MODEL_ROOT="models/vessel" qsub -V -l highmem -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="{'g2da','h2da'}" WIN_SIZE=3 PROBABILITY_DIR="88126/" MAKE_RESAMPLING_MAPS=0.75 MAX_N_IMAGES=10 SHIFT_IMAGES=10 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="retinograms/STARE/training" NUM_SAMPLES=100000 MODEL_ROOT="models/vessel" qsub -V -l highmem -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh

MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88359'" qsub -V -hold_jid 88359 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88360'" qsub -V -hold_jid 88360 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88361'" qsub -V -hold_jid 88361 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88362'" qsub -V -hold_jid 88362 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh

MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88335'" NUM_JOBS=20 IMAGE_ROOT="retinograms/STARE/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -t 11-20 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88336'" NUM_JOBS=20 IMAGE_ROOT="retinograms/STARE/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -t 11-20 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88337'" NUM_JOBS=20 IMAGE_ROOT="retinograms/STARE/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -t 11-20 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh

MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88359'" NUM_JOBS=20 IMAGE_ROOT="retinograms/STARE/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 88363 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88360'" NUM_JOBS=20 IMAGE_ROOT="retinograms/STARE/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 88364 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88361'" NUM_JOBS=20 IMAGE_ROOT="retinograms/STARE/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 88365 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88362'" NUM_JOBS=20 IMAGE_ROOT="retinograms/STARE/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 88366 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh

%%
DECOMP_TYPE="gabor" WIN_SIZE=3 FEATURE_TYPE="real_imag" PROBABILITY_DIR="24066/" MAKE_RESAMPLING_MAPS=0.75 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="retinograms/DRIVE/training" NUM_SAMPLES=100000 MODEL_ROOT="models/vessel" qsub -V -l highmem -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88122'" qsub -V -hold_jid 88122 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'88122'" NUM_JOBS=10 IMAGE_ROOT="retinograms/DRIVE/test" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 88123 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh

DECOMP_TYPE="{'g2da','h2da'}" BG_RATIO=2 WIN_SIZE=3 PROBABILITY_DIR="11476/" MAKE_RESAMPLING_MAPS=0.75 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="retinograms/DRIVE/training" NUM_SAMPLES=100000 MODEL_ROOT="models/vessel" qsub -V -l highmem -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'93293'" qsub -V -hold_jid 93293 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'93293'" NUM_JOBS=10 IMAGE_ROOT="retinograms/DRIVE/test" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 93294 -t 1-10 -l highmem matlab_code/trunk/hydra/cuc/predict_image_set.sh
%%
DECOMP_TYPE="gabor" WIN_SIZE=3 FEATURE_TYPE="real_imag" OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="fibre/training" FG_MASK_DIR="fibre_masks" NUM_SAMPLES=100000 MAX_N_IMAGES=100 MODEL_ROOT="models/fibre" qsub -V -l twoday -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
MODEL_ROOT="models/fibre/detection/rf_classification" MODEL_PATH="'88129'" qsub -V -hold_jid 88129 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/fibre/detection/rf_classification" MODEL_PATH="'88129'" NUM_JOBS=100 IMAGE_ROOT="fibre/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 88130 -t 1-50 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh

DECOMP_TYPE="gabor" WIN_SIZE=3 FEATURE_TYPE="conj" OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="fibre/training" FG_MASK_DIR="fibre_masks" NUM_SAMPLES=100000 MAX_N_IMAGES=100 MODEL_ROOT="models/fibre" qsub -V -l twoday -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="dt" WIN_SIZE=3 FEATURE_TYPE="conj" OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="fibre/training" FG_MASK_DIR="fibre_masks" NUM_SAMPLES=100000 MAX_N_IMAGES=100 MODEL_ROOT="models/fibre" qsub -V -l twoday -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="{'g2da','h2da'}" NUM_ANGLES=6 WIN_SIZE=3 OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="fibre/training" FG_MASK_DIR="fibre_masks" NUM_SAMPLES=100000 MAX_N_IMAGES=100 MODEL_ROOT="models/fibre" qsub -V -l twoday -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="mono" WIN_SIZE=3 OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="fibre/training" FG_MASK_DIR="fibre_masks" NUM_SAMPLES=100000 MAX_N_IMAGES=100 MODEL_ROOT="models/fibre" qsub -V -l twoday -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh

MODEL_ROOT="models/fibre/orientation/rf_regression" MODEL_PATH="'96675'" qsub -V -hold_jid 96675 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/fibre/orientation/rf_regression" MODEL_PATH="'96676'" qsub -V -hold_jid 96676 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/fibre/orientation/rf_regression" MODEL_PATH="'96677'" qsub -V -hold_jid 96677 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/fibre/orientation/rf_regression" MODEL_PATH="'96678'" qsub -V -hold_jid 96678 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh

MODEL_ROOT="models/fibre/detection/rf_classification" MODEL_PATH="'88348'" qsub -V -hold_jid 88348 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/fibre/detection/rf_classification" MODEL_PATH="'88301'" qsub -V -hold_jid 88301 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/fibre/detection/rf_classification" MODEL_PATH="'88302'" qsub -V -hold_jid 88302 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh

MODEL_ROOT="models/fibre/detection/rf_classification" MODEL_PATH="'88348'" NUM_JOBS=100 IMAGE_ROOT="fibre/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 88349 -t 1-50 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/fibre/detection/rf_classification" MODEL_PATH="'88301'" NUM_JOBS=100 IMAGE_ROOT="fibre/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 88304 -t 1-50 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/fibre/detection/rf_classification" MODEL_PATH="'88302'" NUM_JOBS=100 IMAGE_ROOT="fibre/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 88305 -t 1-50 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh

MODEL_ROOT="models/fibre/orientation/rf_regression" MODEL_PATH="'96675'" NUM_JOBS=20 IMAGE_ROOT="fibre/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=1 qsub -V -hold_jid 96679 -t 1-20 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/fibre/orientation/rf_regression" MODEL_PATH="'96676'" NUM_JOBS=20 IMAGE_ROOT="fibre/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=1 qsub -V -hold_jid 96680 -t 1-20 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/fibre/orientation/rf_regression" MODEL_PATH="'96677'" NUM_JOBS=20 IMAGE_ROOT="fibre/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=1 qsub -V -hold_jid 96681 -t 1-20 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/fibre/orientation/rf_regression" MODEL_PATH="'96678'" NUM_JOBS=20 IMAGE_ROOT="fibre/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=1 qsub -V -hold_jid 96682 -t 1-20 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh

DECOMP_TYPE="gabor" PROBABILITY_DIR="88129/" MAKE_RESAMPLING_MAPS=0.75 WIN_SIZE=3 FEATURE_TYPE="real_imag" OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="fibre/training" FG_MASK_DIR="fibre_masks" NUM_SAMPLES=100000 MAX_N_IMAGES=100 MODEL_ROOT="models/fibre" qsub -V -l highmem -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="dt" PROBABILITY_DIR="88300/" MAKE_RESAMPLING_MAPS=0.75 WIN_SIZE=3 FEATURE_TYPE="real_imag" OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="fibre/training" FG_MASK_DIR="fibre_masks" NUM_SAMPLES=100000 MAX_N_IMAGES=100 MODEL_ROOT="models/fibre" qsub -V -l highmem -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="{'g2da','h2da'}" PROBABILITY_DIR="88301/" MAKE_RESAMPLING_MAPS=0.75 NUM_ANGLES=6 WIN_SIZE=3 FEATURE_TYPE="real_imag" OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="fibre/training" FG_MASK_DIR="fibre_masks" NUM_SAMPLES=100000 MAX_N_IMAGES=100 MODEL_ROOT="models/fibre" qsub -V -l highmem -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="mono" PROBABILITY_DIR="88302/" MAKE_RESAMPLING_MAPS=0.75 WIN_SIZE=3 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="fibre/training" FG_MASK_DIR="fibre_masks" NUM_SAMPLES=100000 MAX_N_IMAGES=100 MODEL_ROOT="models/fibre" qsub -V -l highmem -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh

MODEL_ROOT="models/fibre/detection/rf_classification" MODEL_PATH="'88343'" qsub -V -hold_jid 88343 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/fibre/detection/rf_classification" MODEL_PATH="'88344'" qsub -V -hold_jid 88344 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/fibre/detection/rf_classification" MODEL_PATH="'88345'" qsub -V -hold_jid 88345 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh

MODEL_ROOT="models/fibre/detection/rf_classification" MODEL_PATH="'88311'" NUM_JOBS=100 IMAGE_ROOT="fibre/test" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 88312 -t 1-50 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/fibre/detection/rf_classification" MODEL_PATH="'88343'" NUM_JOBS=100 IMAGE_ROOT="fibre/test" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 88346 -t 1-50 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/fibre/detection/rf_classification" MODEL_PATH="'88344'" NUM_JOBS=100 IMAGE_ROOT="fibre/test" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 88347 -t 1-50 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/fibre/detection/rf_classification" MODEL_PATH="'88345'" NUM_JOBS=100 IMAGE_ROOT="fibre/test" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 88348 -t 1-50 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh

DECOMP_TYPE="gabor" PROBABILITY_DIR="96675/" MAKE_RESAMPLING_MAPS=0.75 WIN_SIZE=3 FEATURE_TYPE="real_imag" OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="fibre/training" FG_MASK_DIR="fibre_masks" NUM_SAMPLES=100000 MAX_N_IMAGES=100 MODEL_ROOT="models/fibre" qsub -V -l highmem -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="dt" PROBABILITY_DIR="96676/" MAKE_RESAMPLING_MAPS=0.75 WIN_SIZE=3 FEATURE_TYPE="real_imag" OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="fibre/training" FG_MASK_DIR="fibre_masks" NUM_SAMPLES=100000 MAX_N_IMAGES=100 MODEL_ROOT="models/fibre" qsub -V -l highmem -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="{'g2da','h2da'}" PROBABILITY_DIR="96677/" MAKE_RESAMPLING_MAPS=0.75 NUM_ANGLES=6 WIN_SIZE=3 FEATURE_TYPE="real_imag" OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="fibre/training" FG_MASK_DIR="fibre_masks" NUM_SAMPLES=100000 MAX_N_IMAGES=100 MODEL_ROOT="models/fibre" qsub -V -l highmem -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="mono" PROBABILITY_DIR="96678/" MAKE_RESAMPLING_MAPS=0.75 WIN_SIZE=3 OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="fibre/training" FG_MASK_DIR="fibre_masks" NUM_SAMPLES=100000 MAX_N_IMAGES=100 MODEL_ROOT="models/fibre" qsub -V -l highmem -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh

MODEL_ROOT="models/fibre/orientation/rf_regression" MODEL_PATH="'99189'" qsub -V -hold_jid 99189 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/fibre/orientation/rf_regression" MODEL_PATH="'99190'" qsub -V -hold_jid 99190 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/fibre/orientation/rf_regression" MODEL_PATH="'99191'" qsub -V -hold_jid 99191 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/fibre/orientation/rf_regression" MODEL_PATH="'99192'" qsub -V -hold_jid 99192 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh

MODEL_ROOT="models/fibre/orientation/rf_regression" MODEL_PATH="'99189'" NUM_JOBS=2 IMAGE_ROOT="fibre/test" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 MAX_SIZE=512 qsub -V -t 2 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/fibre/orientation/rf_regression" MODEL_PATH="'99190'" NUM_JOBS=2 IMAGE_ROOT="fibre/test" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 MAX_SIZE=512 qsub -V -t 2 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/fibre/orientation/rf_regression" MODEL_PATH="'99191'" NUM_JOBS=2 IMAGE_ROOT="fibre/test" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 MAX_SIZE=512 qsub -V -t 2 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/fibre/orientation/rf_regression" MODEL_PATH="'99192'" NUM_JOBS=2 IMAGE_ROOT="fibre/test" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 MAX_SIZE=512 qsub -V -t 2 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh

%%
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'24303'" NUM_JOBS=20 IMAGE_ROOT="retinograms/STARE/training" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -t 3 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
%%
DECOMP_TYPE="{'g2da','h2da'}" NUM_ANGLES=6 WIN_SIZE=3 PROBABILITY_DIR="12244/" MAKE_RESAMPLING_MAPS=0.75 OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="retinograms/DRIVE/training" NUM_SAMPLES=100000 MODEL_ROOT="models/vessel" qsub -V -l twoday -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
MODEL_ROOT="models/vessel/orientation/rf_regression" MODEL_PATH="'88171'" qsub -V -hold_jid 88171 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/vessel/orientation/rf_regression" MODEL_PATH="'88171'" NUM_JOBS=10 IMAGE_ROOT="retinograms/DRIVE/test" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 88172 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
%%
for ii = 12:15
     p1 = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\predictions\detection\rf_classification\9162\' zerostr(ii,2) '_test_pred.mat']);
     p2 = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\predictions\detection\rf_classification\12394\' zerostr(ii,2) '_test_pred.mat']);
     p3 = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\predictions\detection\rf_classification\88126\' zerostr(ii,2) '_test_pred.mat']);
     p4 = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\predictions\detection\rf_classification\88122\' zerostr(ii,2) '_test_pred.mat']);

     p1 = [p1(2:end,:); p1(end,:)];
     p1 = [p1(:,2:end), p1(:,end)];
     
     figure; 
     subplot(2,3,1); imgray(p1);
     subplot(2,3,2); imgray(p2);
     subplot(2,3,3); imgray(p1-p3);
     subplot(2,3,4); imgray(p3);
     subplot(2,3,5); imgray(p4);
     subplot(2,3,6); imgray(p2-p3);
end
%%
for ii = 14
     p1 = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\predictions\detection\rf_classification\9162\' zerostr(ii,2) '_test_pred.mat']);
     p2 = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\predictions\detection\rf_classification\88126\' zerostr(ii,2) '_test_pred.mat']);
     
     figure; imgray(p1-p2);
end
%%
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'9162'" NUM_JOBS=20 IMAGE_ROOT="retinograms/DRIVE/test" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -t 14 -l highmem matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'12394'" NUM_JOBS=20 IMAGE_ROOT="retinograms/DRIVE/test" MASK_DIR="fov_masks" USE_SAMPLED_MAPS=0 qsub -V -t 3 -l highmem matlab_code/trunk/hydra/cuc/predict_image_set.sh
%%
mkdir C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\predictions\detection\staal
for ii = 1:20
    vessel_prob = double(rgb2gray(imread(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\predictions\staal\' zerostr(ii-1, 2) '.bmp'])))/255;
    save(['C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\predictions\detection\staal\' zerostr(ii, 2) '_test_pred.mat'], 'vessel_prob');
end