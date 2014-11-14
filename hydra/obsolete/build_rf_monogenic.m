function [] = build_rf_monogenic(job_idx, n_samples, n_trees, d)
%BUILD_RF *Insert a one line summary here*
%   [] = build_rf(job_idx,image_dir)
%
% Inputs:
%      job_idx- *Insert description of input variable here*
%
%      image_dir- *Insert description of input variable here*
%
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 03-Feb-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%randomise Matlab seed (otheriwse each parallel job will build the same
%trees!
rand('twister', sum(100*clock));

%Set arguments to sample data
sampling_method_args.num_samples = n_samples;
sampling_method_args.width_range = [4 16];
sampling_method_args.contrast_range = [1 8];
sampling_method_args.win_size = 3;
sampling_method_args.num_levels = 5;
sampling_method_args.num_angles = 8;
sampling_method_args.bg_dir = [mberksroot, 'classification/data/smooth512x512_patches/train/'];

%set up arguments for main forest
forest_args.sampling_method = 'sample_training_data_linop';
forest_args.sampling_method_args = sampling_method_args;
forest_args.d = d;
forest_args.n_trees = n_trees;
forest_args.tree_root = mberksroot;
forest_args.tree_dir = ['classification/rf/rf_linop_3_5_8_job', zerostr(job_idx,2), '_trees/'];
forest_args.save_path = [mberksroot, 'classification/rf/rf_linop_3_5_8_job', zerostr(job_idx,2), '.mat'];

%build forest
mb_random_forest_class_train(forest_args);

