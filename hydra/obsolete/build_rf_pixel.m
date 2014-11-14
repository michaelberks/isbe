function [] = build_rf_pixel(varargin)
%BUILD_RF_PIXEL wrapper function to build an random forest line
%detector on hydra
%   [] = build_rf_pixel(win_size, num_levels, do_max, feature_type, num_trees)
%
% BUILD_RF_PIXEL uses the U_PACKARGS interface function
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
    job_id = ['pc' datestr(now, 'yyyymmddTHHMMSS')];
    
    %Choose between line classification and orientation regression
    detection_type = 'detection';
    
    %arguments controlling data sampling
    sampling_method = 'generate_pixel_training_data';
    num_samples = 2e5;
    pts_per_image = 500;
    
    bg_stem = 'bg'; %args for background patches
    num_bgs = 10460;
    bg_zeros = 5;
    bg_dir = 'smooth512/train/';
    bg_ratio = 1;
    
    line_type = 'sin'; %args for generating lines
    width_range = [4 16];
    contrast_range = [1 8];
    decay_rate = 4;

    %args for patch represententation
    win_size = 15;
    
    %arguments controlling forest construction
    n_trees = 200;
    split_min = 100;
    overwrite = 0;
else
%Set defaults for unix system - environment variables listed in .bashrc (and utils/mb_bash.m)
    %task and job identifiers
    [z task_id] = unix('echo $SGE_TASK_ID'); task_id = str2num(task_id); %#ok
    [z job_id] = unix('echo $CUSTOM_ID'); job_id = str2num(job_id); %#ok
    if isempty(job_id)
        [z job_id] = unix('echo $JOB_ID'); job_id = str2num(job_id); %#ok
    end
    
    %Choose between line classification and orientation regression
    [z detection_type] = unix('echo $DETECTION_TYPE'); detection_type(end) = [];
    
    %arguments controliing data sampling
    [z sampling_method] = unix('echo $SAMPLING_METHOD_P'); sampling_method(end) = [];
    [z num_samples] = unix('echo $NUM_SAMPLES'); num_samples = str2num(num_samples); %#ok
    [z pts_per_image] = unix('echo $PTS_PER_IMAGE'); pts_per_image = str2num(pts_per_image); %#ok
    
    [z bg_stem] = unix('echo $BG_STEM'); bg_stem(end) = []; %args for background patches
    [z num_bgs] = unix('echo $NUM_BGS'); num_bgs = str2num(num_bgs); %#ok
    [z bg_zeros] = unix('echo $BG_ZEROS'); bg_zeros = str2num(bg_zeros); %#ok
    [z bg_dir] = unix('echo $BG_DIR'); bg_dir(end) = [];
    [z bg_ratio] = unix('echo $BG_RATIO'); bg_ratio = str2num(bg_ratio); %#ok
    
    [z line_type] = unix('echo $BAR_TYPE'); line_type(end) = []; %args for generating lines
    [z width_range] = unix('echo $WIDTH_RANGE'); width_range = str2num(width_range); %#ok
    [z contrast_range] = unix('echo $CONTRAST_RANGE'); contrast_range = str2num(contrast_range); %#ok
    [z decay_rate] = unix('echo $DECAY_RATE'); decay_rate = str2num(decay_rate); %#ok
     
    [z win_size] = unix('echo $WIN_SIZE'); win_size = str2num(win_size); %#ok %args for patch represententation
    
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
    'job_id', job_id, ...
    'detection_type', detection_type, ...
    'num_samples', num_samples, ...
    'pts_per_image', pts_per_image, ...
    'win_size', win_size, ...
    'width_range', width_range, ...
    'contrast_range', contrast_range, ...
    'decay_rate', decay_rate, ...
    'line_type', line_type, ...
    'bg_stem', bg_stem, ...
    'num_bgs', num_bgs, ...
    'bg_zeros', bg_zeros, ...
    'bg_dir', bg_dir, ...
    'bg_ratio', bg_ratio, ...
    'sampling_method', sampling_method, ...
    'n_trees', n_trees,...
    'split_min', split_min,...
    'overwrite', overwrite);
    
job_id = args.job_id;

sampling_method_args.bg_dir = [asymmetryroot, 'data/synthetic_backgrounds/' args.bg_dir];
sampling_method_args.save_path = [asymmetryroot 'data/line_' args.detection_type '_rfs/' num2str(args.job_id) '/line_parameters/' zerostr(args.task_id,2) '/'];

%Set arguments to sample data
sampling_method_args.detection_type = args.detection_type;
sampling_method_args.num_samples = args.num_samples;
sampling_method_args.pts_per_image = args.pts_per_image;

sampling_method_args.bg_stem = args.bg_stem; %args for background patches
sampling_method_args.num_bgs = args.num_bgs;
sampling_method_args.bg_zeros = args.bg_zeros;
sampling_method_args.bg_ratio = args.bg_ratio;

sampling_method_args.line_type = args.line_type; %args for generating lines
sampling_method_args.width_range = args.width_range;
sampling_method_args.contrast_range = args.contrast_range;
sampling_method_args.decay_rate = args.decay_rate;
 
sampling_method_args.win_size = args.win_size; %args for patch represententation

%set up arguments for main forest
forest_args.sampling_method = args.sampling_method;
forest_args.sampling_method_args = sampling_method_args;

forest_args.n_trees = args.n_trees;
forest_args.split_min = args.split_min;
forest_args.overwrite = args.overwrite;
%forest_args.d = round(args.win_size * sqrt(2*args.num_levels * (6 -
%5*args.do_max))); %let this be selected automatically

forest_args.tree_dir = [num2str(job_id) '/' zerostr(args.task_id,2), '_trees/'];
forest_args.save_path = [asymmetryroot 'data/line_' args.detection_type '_rfs/' num2str(job_id) '/random_forest' zerostr(args.task_id,2) '.mat'];
forest_args.tree_root = [asymmetryroot 'data/line_' args.detection_type '_rfs/'];

%set random seed based on task id and clock to ensure unique trees created
rand('twister', sum(args.task_id*clock));

%Clear the args structure and display the sampling/forest args
clear args;
display(sampling_method_args);
display(forest_args);

%build forest
if strcmpi(sampling_method_args.detection_type, 'orientation')
    forest_args.mod = 180;
    mb_random_forest_reg_train(forest_args);
else   
    mb_random_forest_class_train(forest_args);
end

%For the first task of job, save the sampling args
args_name = [asymmetryroot 'data/line_' sampling_method_args.detection_type '_rfs/' num2str(job_id) '/sampling_args.mat'];
if ~exist(args_name, 'file');
    save(args_name, 'sampling_method_args');
end
display('Forest successfully constructed!');

