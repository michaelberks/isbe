function [] = build_rf_mammo_class(varargin)
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

%Set defaults for windows system
if ispc
    
    %task and job identifiers
    task_id = 1;
    custom_id = ['pc' datestr(now, 'yyyymmddTHHMMSS')];
    fold_id = 1;
    num_folds = 10;
    
    %arguments controlling data sampling
    sampling_method = 'sample_mammo_training_data';
    num_samples = 2e5;
    use_nag = 1;
    abnormal_data = '2004_screening\abnormals';
    normal_data = '2004_screening\normals';
    image_dir = 'mammograms';
    mask_dir = 'masks';
    
    image_type = '.mat';
    view = 'CC';
    
    num_levels = 5; %args for patch represententation
    feature_shape = 'rect';
    feature_type = 'all';
    do_max = 0;
    rotate = 0;
    win_size = 3;
    
    %arguments controlling forest construction
    n_trees = 200;
    split_min = 100;
    overwrite = 0;
else
%Set defaults for unix system - environment variables listed in .bashrc (and utils/mb_bash.m)
    %task and job identifiers
    [z task_id] = unix('echo $SGE_TASK_ID'); task_id = str2num(task_id); %#ok
    [z custom_id] = unix('echo $CUSTOM_ID'); custom_id(end) = [];
    [z fold_id] = unix('echo $FOLD_ID'); fold_id = str2num(fold_id); %#ok
    [z num_folds] = unix('echo $NUM_FOLDS'); num_folds = str2num(num_folds); %#ok
    
    %arguments controliing data sampling
    [z sampling_method] = unix('echo $SAMPLING_METHOD_MAMMO'); sampling_method(end) = [];
    [z num_samples] = unix('echo $NUM_SAMPLES'); num_samples = str2num(num_samples); %#ok
    [z pts_per_image] = unix('echo $PTS_PER_IMAGE'); pts_per_image = str2num(pts_per_image); %#ok
    [z use_nag] = unix('echo $USE_NAG'); use_nag = str2num(use_nag); %#ok
    [z abnormal_data] = unix('echo $ABNORMAL_DATA'); abnormal_data(end) = [];
    [z normal_data] = unix('echo $NORMAL_DATA'); normal_data(end) = [];
    [z image_dir] = unix('echo $IMAGE_DIR'); image_dir(end) = [];
    [z mask_dir] = unix('echo $MASK_DIR'); mask_dir(end) = [];
    [z image_type] = unix('echo $IMAGE_TYPE'); image_type(end) = [];
    [z view] = unix('echo $VIEW'); view(end) = [];
    
    %args for patch represententation
    [z num_levels] = unix('echo $NUM_LEVELS'); num_levels = str2num(num_levels); %#ok 
    [z feature_shape] = unix('echo $FEATURE_SHAPE'); feature_shape(end) = []; 
    [z feature_type] = unix('echo $FEATURE_TYPE'); feature_type(end) = []; %args for patch represententation
    [z win_size] = unix('echo $WIN_SIZE'); win_size = str2num(win_size); %#ok
    [z do_max] = unix('echo $DO_MAX'); do_max = str2num(do_max); %#ok
    [z rotate] = unix('echo $ROTATE'); rotate = str2num(rotate); %#ok
    
    %arguments controlling forest construction
    [z n_trees] = unix('echo $NUM_TREES'); n_trees = str2num(n_trees); %#ok
    [z split_min] = unix('echo $SPLIT_MIN'); split_min = str2num(split_min); %#ok
    [z overwrite] = unix('echo $OVERWRITE'); overwrite = str2num(overwrite); %#ok
    clear z;
end

%Now use varargin to merge default argument values with user set values
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'task_id', task_id, ...
    'custom_id', custom_id, ...
    'fold_id', fold_id, ...
    'num_folds', num_folds, ...
    'abnormal_data', abnormal_data,...
    'normal_data', normal_data,...
    'image_dir', image_dir,...
    'mask_dir', mask_dir,...
    'image_type', image_type,...
    'view', view,...
    'num_samples', num_samples, ...
    'num_levels', num_levels, ...
    'feature_shape', feature_shape, ...
    'feature_type', feature_type, ...
    'win_size', win_size, ...
    'do_max', do_max, ...
    'rotate', rotate, ...
    'sampling_method', sampling_method, ...
    'n_trees', n_trees,...
    'split_min', split_min,...
    'overwrite', overwrite,...
    'use_nag', use_nag);
    
fold_id = args.fold_id;
view = args.view;
if isempty(args.custom_id)    
    forest_name = [view '_' zerostr(fold_id,2)];
else
    forest_name = [custom_id '_' zerostr(fold_id,2)];
end

%Set arguments to sample data
sampling_method_args.abnormal_dir = [asymmetryroot, 'data/' args.image_dir '/' args.abnormal_data '/'];
sampling_method_args.normal_dir = [asymmetryroot, 'data/' args.image_dir '/' args.normal_data '/'];
sampling_method_args.abnormal_mask_dir = [asymmetryroot, 'data/' args.mask_dir '/' args.abnormal_data '/'];
sampling_method_args.normal_mask_dir = [asymmetryroot, 'data/' args.mask_dir '/' args.normal_data '/'];
sampling_method_args.fold_id = args.fold_id;
sampling_method_args.num_folds = args.num_folds;
sampling_method_args.view = args.view;
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
    save(args_name, 'sampling_method_args');
end
display('Forest successfully constructed!');

