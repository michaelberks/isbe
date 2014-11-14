function [] = build_rf_allsubbands_pixelset3(JOB_ID, TASK_ID, n_samples, n_trees, d)
%
% [] = build_rf_chen(JOB_ID, TASK_ID, n_samples, n_trees, d)
% USAGE:
% submit a .sh job on Hydra cluster machine
% The describtion is:
% #!/bin/bash
% #
% matlab -nodisplay -nojvm -r "build_rf($JOB_ID, $SGE_TASK_ID, 2e4, 20, 10)" > logs/build_rf$SGE_TASK_ID.log
%
% Inputs:
%       JOB_ID - the job ID from cluster machine
%       TASK_ID - the task ID from cluster machine
%                ------- * *  When sue qsub to submit multiple job, the machine return a JOB ID and a TASK ID. * * ------
%
%       n_samples - the total number of observation vectors in each JOB
%       n_trees -  the number of trees to grow in each JOB
%       d - number of variables randomly sampled as candidate at each split.
%           the default values are sqrt(p) where p is the number of variables.
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
% Created: 09-Feb-2010
% Author: Zezhi Chen and Michael Berks
% Email : zezhi.chen@manchester.ac.uk
% Phone : +44 (0)161 275 7669
% Copyright: (C) University of Manchester
%
%
% NOTE: randomise Matlab seed (otheriwse each parallel job will build the same
%trees!


rand('twister', sum(100*clock));

if strcmpi(computer,'PCWIN') |strcmpi(computer,'PCWIN64')
    %% windows path
    % the root of random forests
    rfroot_chen = 'E:\DTCWTmini\data\';
    % the path of backgound image patch
    sampling_method_args.bg_dir = [rfroot_chen, 'smooth512x512_patches\train\'];
      % the path of parameters of each tree in the random forest 
    sampling_method_args.save_path = [rfroot_chen, 'rf_trees_', int2str(JOB_ID), '\rf_centre_', zerostr(TASK_ID,2), '_parameters\'];
   % the path of random forests 
    forest_args.tree_dir = ['rf_trees_', int2str(JOB_ID), '\rf_centre_', zerostr(TASK_ID,2), '_trees\'];
    % the path of a forest of each TASK from Hydra
    forest_args.save_path = [rfroot_chen, 'rf_trees_', int2str(JOB_ID), '\rf_centre_', zerostr(TASK_ID,2), '_forest.mat'];
else
    %% Linux path  
    rfroot_chen = '/data/zchen/';
    sampling_method_args.bg_dir = [rfroot_chen, 'smooth512x512_patches/train/'];
    sampling_method_args.save_path = [rfroot_chen, 'rf_trees_', int2str(JOB_ID), '/rf_centre_', zerostr(TASK_ID,2), '_parameters/'];
    forest_args.tree_dir = ['rf_trees_', int2str(JOB_ID), '/rf_centre_', zerostr(TASK_ID,2), '_trees/'];
    forest_args.save_path = [rfroot_chen, 'rf_trees_', int2str(JOB_ID), '/rf_centre_', zerostr(TASK_ID,2), '_forest.mat'];
end


%Set arguments to sample data
sampling_method_args.num_samples = n_samples;
sampling_method_args.width_range = [4 16];
sampling_method_args.contrast_range = [4 16];
sampling_method_args.num_levels = 3;
% sampling_method_args.bar_type = 'sin';
% sampling_method_args.feature_type = 'phase';

%set up arguments for main forest
forest_args.sampling_method = 'sample_training_data_pixelset';
forest_args.sampling_method_args = sampling_method_args;
forest_args.d = d;
forest_args.n_trees = n_trees;
forest_args.tree_root = rfroot_chen;

%build forest
mb_random_forest_class_train(forest_args);

