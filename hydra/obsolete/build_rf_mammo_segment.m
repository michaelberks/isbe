function [] = build_rf_mammo_segment(varargin)
%BUILD_RF_LINE_DETECTOR wrapper function to build an random forest line
%detector on hydra
%   [] = build_rf_line_detector(win_size, num_levels, do_max, feature_type, num_trees)
%
% BUILD_RF_LINE_DETECTOR uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names.
%
%
% Outputs:
%
% Example: Usage on hydra with qsub
%   NUM_TREES=20 NUM_BGS=1000 NUM_SAMPLES=200000 WIN_SIZE=1 qsub -N rf_line -t 1:10 -V matlab_code/trunk/hydra/build_rf_line_detector.sh                              Your job-array 191656.1-10:1 ("rf_line") has been submitted
%   FOREST_JOB="'191656'" qsub -N comb -hold_jid 191656 -V matlab_code/trunk/hydra/combine_hydra_rfs.sh
%
% Notes:
%
% See also:
%
% Created: 03-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%OK, there's some jiggery-pokery here to allow for optional setting of
%default variables through a unix shell (i.e. when running on hydra)

%Use varargin to merge default argument values with user set values
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'task_id', unixenv('SGE_TASK_ID',1), ...
    'custom_id', unixenv('CUSTOM_ID',['pc' datestr(now, 'yyyymmddTHHMMSS')]), ...
    'num_mammos', unixenv('NUM_FOLDS',12), ...
    'image_dir', unixenv('IMAGE_DIR','mammograms/2004_screening/abnormals'),...
    'mask_dir', unixenv('MASK_DIR','masks/2004_screening/abnormals'),...
    'image_type', unixenv('IMAGE_TYPE','.mat'),...
    'view', unixenv('VIEW',[]),...
    'resize', unixenv('RESIZE', 0.25),...
    'num_samples', unixenv('NUM_SAMPLES',2e3), ...
    'num_levels', unixenv('NUM_LEVELS',5), ...
    'feature_shape', unixenv('FEATURE_SHAPE','rect'), ...
    'feature_type', unixenv('FEATURE_TYPE','conj'), ...
    'win_size', unixenv('WIN_SIZE',1), ...
    'do_max', unixenv('DO_MAX',0), ...
    'rotate', unixenv('ROTATE',0), ...
    'sampling_method', unixenv('SAMPLING_METHOD_MAMMO','sample_mammo_segment_data'), ...
    'n_trees', unixenv('NUM_TREES',2),...
    'split_min', unixenv('SPLIT_MIN',100),...
    'overwrite', unixenv('OVERWRITE',0),...
    'use_nag', unixenv('USE_NAG',0));

if isempty(args.custom_id)    
    forest_name = ['segment_' args.view];
else
    forest_name = args.custom_id;
end

%Set arguments to sample data
sampling_method_args.mammo_dir = [asymmetryroot, 'data/' args.image_dir '/'];
sampling_method_args.mask_dir = [asymmetryroot, 'data/' args.mask_dir '/'];
sampling_method_args.num_mammos = args.num_mammos;
sampling_method_args.view = args.view;
sampling_method_args.resize = args.resize;
sampling_method_args.image_type = args.image_type;

sampling_method_args.num_samples = args.num_samples;
sampling_method_args.use_nag = args.use_nag;
sampling_method_args.num_levels = args.num_levels;
sampling_method_args.feature_shape = args.feature_shape;
sampling_method_args.feature_type = args.feature_type;
sampling_method_args.do_max = args.do_max;
sampling_method_args.rotate = args.rotate;
sampling_method_args.win_size = args.win_size;

%set up arguments for main forest
forest_args.sampling_method = args.sampling_method;
forest_args.sampling_method_args = sampling_method_args;

forest_args.n_trees = args.n_trees;
forest_args.split_min = args.split_min;
forest_args.overwrite = args.overwrite;

forest_args.tree_dir = [forest_name '/' zerostr(args.task_id,2), '_trees/'];
forest_args.save_path = [asymmetryroot 'data/mammo_rfs/' forest_name '/random_forest' zerostr(args.task_id,2) '.mat'];
forest_args.tree_root = [asymmetryroot 'data/mammo_rfs/'];

%set random seed based on task id and clock to ensure unique trees created
rand('twister', sum(args.task_id*clock));
randn('state', sum(args.task_id*clock));

%Clear the args structure and display the sampling/forest args
clear args;
display(sampling_method_args);
display(forest_args);

%build forest
mb_random_forest_class_train(forest_args);

%For the first task of job, save the sampling args
args_name = [asymmetryroot 'data/mammo_rfs/' forest_name '/sampling_args.mat'];
if ~exist(args_name, 'file');
    sampling_args.detection_type = 'classification';
    save(args_name, 'sampling_method_args');
end
display('Forest successfully constructed!');

