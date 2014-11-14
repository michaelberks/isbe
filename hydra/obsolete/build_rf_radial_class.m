function [] = build_rf_radial_class(varargin)
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
    %sampling_method = 'sample_mass_training_data';
    sampling_method = 'load_mass_training_data';
    num_samples = 2e3;
    abnormal_data = '2004_screening/abnormals';
    normal_data = '2004_screening/normals';
    image_dir = 'mammograms';
    mask_dir = 'mass_masks';
    radial_dir = 'weighted_radial_maps';
    template_dir = 'template_maps';
    do_template = 1;
    do_scale = 1;
    
    image_type = '.mat';
    view = [];
    
    %arguments controlling feature vector composition
    dist_range = [16 32 64 128 256];
    sigma_range = [1 2 4 8];
    angular_res = 1;
    
    %arguments controlling forest construction
    n_trees = 2;
    split_min = 100;
    overwrite = 0;
    rand_seed = [];
else
%Set defaults for unix system - environment variables listed in .bashrc (and utils/mb_bash.m)
    %task and job identifiers
    [z task_id] = unix('echo $SGE_TASK_ID'); task_id = str2num(task_id); %#ok
    [z custom_id] = unix('echo $CUSTOM_ID'); custom_id(end) = [];
    [z fold_id] = unix('echo $FOLD_ID'); fold_id = str2num(fold_id); %#ok
    [z num_folds] = unix('echo $NUM_FOLDS'); num_folds = str2num(num_folds); %#ok
    
    %arguments controliing data sampling
    [z sampling_method] = unix('echo $SAMPLING_METHOD_RADIAL'); sampling_method(end) = [];
    [z num_samples] = unix('echo $NUM_SAMPLES'); num_samples = str2num(num_samples); %#ok
    [z pts_per_image] = unix('echo $PTS_PER_IMAGE'); pts_per_image = str2num(pts_per_image); %#ok
    [z abnormal_data] = unix('echo $ABNORMAL_DATA'); abnormal_data(end) = [];
    [z normal_data] = unix('echo $NORMAL_DATA'); normal_data(end) = [];
    [z image_dir] = unix('echo $IMAGE_DIR'); image_dir(end) = [];
    [z mask_dir] = unix('echo $MASK_DIR'); mask_dir(end) = [];
    [z radial_dir] = unix('echo $RADIAL_DIR'); radial_dir(end) = [];
    [z template_dir] = unix('echo $TEMPLATE_DIR'); template_dir(end) = [];
    [z image_type] = unix('echo $IMAGE_TYPE'); image_type(end) = [];
    [z view] = unix('echo $VIEW'); view(end) = [];
    [z do_template] = unix('echo $DO_TEMPLATE'); do_template = str2num(do_template); %#ok
    [z do_scale] = unix('echo $DO_SCALE'); do_scale = str2num(do_scale); %#ok
    
    %args for patch represententation
    [z dist_range] = unix('echo $DISTANCE_RANGE'); dist_range = str2num(dist_range); %#ok 
    [z sigma_range] = unix('echo $SIGMA_RANGE'); sigma_range = str2num(sigma_range); %#ok
    [z angular_res] = unix('echo $ANGULAR_BANDS'); angular_res = str2num(angular_res); %#ok
    
    %arguments controlling forest construction
    [z n_trees] = unix('echo $NUM_TREES'); n_trees = str2num(n_trees); %#ok
    [z split_min] = unix('echo $SPLIT_MIN'); split_min = str2num(split_min); %#ok
    [z overwrite] = unix('echo $OVERWRITE'); overwrite = str2num(overwrite); %#ok
    [z rand_seed] = unix('echo $RAND_SEED'); rand_seed = str2num(rand_seed); %#ok
    clear z;
end

%Now use varargin to merge default argument values with user set values
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'task_id', task_id, ...
    'custom_id', custom_id, ...
    'fold_id', fold_id, ...
    'num_folds', num_folds, ...
    'num_samples', num_samples,...
    'abnormal_data', abnormal_data,...
    'normal_data', normal_data,... 
    'image_dir', image_dir,...
    'mask_dir', mask_dir,...
    'radial_dir', radial_dir,...
    'template_dir', template_dir,...
    'view', view,...
    'do_template', do_template,...
    'do_scale', do_scale,...
    'dist_range', dist_range,...
    'sigma_range', sigma_range,...
    'angular_res', angular_res,...
    'image_type', image_type,...
    'sampling_method', sampling_method, ...
    'n_trees', n_trees,...
    'split_min', split_min,...
    'overwrite', overwrite,...
    'rand_seed', rand_seed);
    
task_id = args.task_id;
fold_id = args.fold_id;
view = args.view;
if isempty(args.custom_id)    
    forest_name = [view '_' zerostr(fold_id,2)];
else
    forest_name = [args.custom_id '_' zerostr(fold_id,2)];
end

%Set arguments to sample data
sampling_method_args.abnormal_data  = args.abnormal_data;  
sampling_method_args.normal_data  = args.normal_data; 
sampling_method_args.image_dir = [asymmetryroot, 'data/' args.image_dir '/'];
sampling_method_args.save_dir = [asymmetryroot, 'data/' args.radial_dir '/'];
sampling_method_args.save_name = args.custom_id;
% sampling_method_args.mask_dir = [asymmetryroot, 'data/' args.mask_dir '/'];
sampling_method_args.radial_dir = [asymmetryroot, 'data/' args.radial_dir '/'];
sampling_method_args.template_dir = [asymmetryroot, 'data/' args.template_dir '/'];
sampling_method_args.fold_id = args.fold_id;
sampling_method_args.num_folds = args.num_folds;
sampling_method_args.view = args.view;
sampling_method_args.do_template = args.do_template;
sampling_method_args.do_scale = args.do_scale;
sampling_method_args.dist_range = args.dist_range;
sampling_method_args.sigma_range = args.sigma_range;
sampling_method_args.angular_res = args.angular_res;
sampling_method_args.image_type = args.image_type;
sampling_method_args.num_samples = args.num_samples;

%set up arguments for main forest
forest_args.sampling_method = args.sampling_method;
forest_args.sampling_method_args = sampling_method_args;

forest_args.n_trees = args.n_trees;
forest_args.split_min = args.split_min;
forest_args.overwrite = args.overwrite;

forest_args.tree_dir = [forest_name '/' zerostr(task_id,2), '_trees/'];
forest_args.save_path = [asymmetryroot 'data/radial_rfs/' forest_name '/random_forest' zerostr(task_id,2) '.mat'];
forest_args.tree_root = [asymmetryroot 'data/radial_rfs/'];

%set random seed based on task id and clock to ensure unique trees created
if isempty(args.rand_seed)
    rand_seed = sum(task_id*clock);
else
    rand_seed = args.rand_seed;
end
rand('twister', rand_seed);
randn('state', rand_seed);

%Clear the args structure and display the sampling/forest args
clear args;
display(sampling_method_args);
display(forest_args);

%build forest
mb_random_forest_class_train(forest_args);

%For the first task of job, save the sampling args
args_name = [asymmetryroot 'data/radial_rfs/' forest_name '/sampling_args.mat'];
if ~exist(args_name, 'file');
    sampling_method_args.rand_seed = rand_seed;
    sampling_method_args.task_id = task_id;
    save(args_name, 'sampling_method_args');
end
display('Forest successfully constructed!');

