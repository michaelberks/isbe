function [predictor] = build_rf_line_detector(varargin)
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

warning('off', 'ASYM:unexpectedArgument');

%% use varargin to merge default argument values with user set values

% this now uses the utils/unixenv function that takes the key name of the
% environment variable and returns its corresponding value, returning the
% default value if the environment variable does not exist
%   e.g. envval = unixenv('ENVVARNAME',default_val);

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'task_id',				unixenv('SGE_TASK_ID',1), ...
    'job_id',				unixenv('JOB_ID',['pc',datestr(now,'yyyymmddTHHMMSS')]), ...
    'detection_type',		unixenv('DETECTION_TYPE','detection'), ...
	'num_samples',			unixenv('NUM_SAMPLES',2e3), ...
    'pts_per_image',		unixenv('PTS_PER_IMAGE',500), ...
    'normalise', 			unixenv('NORMALISE',0), ...
    'num_levels', 			unixenv('NUM_LEVELS',5), ...
    'decomp_type', 			unixenv('DECOMP_TYPE','dt'), ...
    'feature_shape', 		unixenv('FEATURE_SHAPE','rect'), ...
    'feature_type',			unixenv('FEATURE_TYPE','conj'), ...
    'win_size',				unixenv('WIN_SIZE',3), ...
    'do_max',				unixenv('DO_MAX',0), ...
    'rotate',				unixenv('ROTATE',0), ...
    'pca',					unixenv('PCA',[]), ...
    'num_angles', 			unixenv('NUM_ANGLES',8), ...
    'min_wavelength',		unixenv('MIN_WAVELENGTH',4), ...
    'onf',					unixenv('ONF',0.65), ...
    'sigma_range', 			unixenv('SIGMA_RANGE',[1 2 4 8]), ...
    'width_range', 			unixenv('WIDTH_RANGE',[4 16]), ...
    'contrast_range',       unixenv('CONTRAST_RANGE',[4 8]), ...
    'decay_rate', 			unixenv('DECAY_RATE',4), ...
    'line_type',			unixenv('BAR_TYPE','sin'), ...
    'bg_stem',				unixenv('BG_STEM','bg'), ...
    'num_bgs',				unixenv('NUM_BGS',100), ...
    'bg_zeros',				unixenv('BG_ZEROS',5), ...
    'bg_dir',				unixenv('BG_DIR','real512/train/'), ...
    'bg_mask_dir',          unixenv('MASK_DIR',[]), ...
	'bg_fmt',				unixenv('BG_FMT','mat'), ...
    'bg_ratio',				unixenv('BG_RATIO',1), ...
    'sampling_method',      unixenv('SAMPLING_METHOD','generate_line_training_data'), ...
    'image_type', 			unixenv('IMAGE_TYPE','line'), ...
    'rgb_channel',          unixenv('RGB_CHANNEL','all'), ...
    'n_trees',				unixenv('NUM_TREES',2), ...
    'split_criterion',      unixenv('SPLIT_CRITERION', 'dabs'),...
    'var_criterion',		unixenv('VAR_CRITERION', 'mabs'),...
    'split_min',			unixenv('SPLIT_MIN',10), ...
    'end_cut_min',			unixenv('END_CUT_MIN',1), ...
    'do_ubound',			unixenv('DO_UBOUND',1), ...
    'do_circular',			unixenv('DO_CIRCULAR',[]), ...
	'w_prior',				unixenv('W_PRIOR',0), ...
    'impure_thresh',		unixenv('IMPURE_THRESH',1e-4), ...
	'minimise_size',		unixenv('MINIMIZE_TREE',0), ...
    'd',                    unixenv('d',[]), ...
    'overwrite',			unixenv('OVERWRITE',0), ...
    'use_nag',				unixenv('USE_NAG',1),...
    'rand_seed',			unixenv('RAND_SEED',[]));

% if $CUSTOM_ID is defined and not empty then use that instead of $JOB_ID
custom_id = unixenv('CUSTOM_ID');
if ~isempty(custom_id)
	args.job_id = custom_id;
end

job_id = args.job_id;

%% copy arguments for sampling
switch args.sampling_method
    case {'generate_line_training_data','generate_training_data'}
        
        args.bg_dir = [asymmetryroot, 'data/synthetic_backgrounds/' args.bg_dir];
        args.bg_mask_dir = [asymmetryroot, 'data/synthetic_backgrounds/' args.bg_mask_dir];
        args.save_path = ...
            [asymmetryroot 'data/line_' args.detection_type '_rfs/' num2str(args.job_id) ...
            '/line_parameters/' zerostr(args.task_id,2) '/'];
    
        fields_to_copy = {...
            'bg_dir','bg_fmt','save_path','pca', ...
            'detection_type','num_samples','pts_per_image','normalise', ... % args to sample data
            'bg_stem','num_bgs','bg_zeros','bg_ratio','bg_mask_dir',... %args for background patches
            'image_type', ... % use line or grain images
            'line_type','width_range','contrast_range','decay_rate', ... %args for generating lines
            'decomp_type'}; %args for patch represententation
        
        % include arguments specific to sampling method
        switch args.decomp_type
            case 'dt',
                fields_to_copy = [fields_to_copy {...
                    'num_levels','feature_shape','feature_type',...
                    'do_max','rotate','win_size','use_nag'}];

            case 'mono'
                fields_to_copy = [fields_to_copy {...
                    'num_levels','win_size','min_wavelength','onf'}];

            case {'g2d', 'g2di', 'clover', 'haar', 'g1d'}
                fields_to_copy = [fields_to_copy {...
                    'sigma_range','win_size'}];

            case 'linop'
                fields_to_copy = [fields_to_copy {...
                    'win_size','num_levels','num_angles','do_max','rotate'}];

            case 'pixel'
                fields_to_copy = [fields_to_copy {...
                    'win_size'}];
        end
        
    case 'sample_saved_dt_line_data'
        args.saved_data_dir = [asymmetryroot, 'data/synthetic_data/' args.bg_dir];
        args.id_offset = (args.task_id-1)*args.n_trees; 
        fields_to_copy = {...
            'num_samples','saved_data_dir', 'detection_type', 'id_offset'...
            'num_levels','feature_shape','feature_type',...
            'do_max','rotate','win_size'};
		if strcmp(get_username,'ptresadern')
			fields_to_copy{end+1} = 'pts_per_image';
        end
        
    case 'sample_vessel_dt_data'
        args.saved_data_dir = [asymmetryroot, 'data/synthetic_data/' args.bg_dir];
        args.id_offset = (args.task_id-1)*args.n_trees; 
        fields_to_copy = {...
            'num_samples','saved_data_dir', 'detection_type', 'id_offset'...
            'num_levels','feature_shape','feature_type',...
            'do_max','rotate','win_size', 'rgb_channel'};
		if strcmp(get_username,'ptresadern')
			fields_to_copy{end+1} = 'pts_per_image';
		end
end

% copy marked fields from args to sampling_method_args
sampling_method_args = get_substructure(args,fields_to_copy);

% load in saved PCA data if necessary
if ~isempty(args.pca)
    sampling_method_args.pca = u_load([asymmetryroot, 'data/pca/' args.pca]);
else
    sampling_method_args.pca = [];
end

%% set up arguments for main forest
predictor_args = get_substructure(args,...
	{	'sampling_method','n_trees','d','w_prior','impure_thresh',...
		'split_criterion','var_criterion','split_min','end_cut_min',...
		'do_ubound', 'do_circular', 'overwrite','minimise_size' });

% manually set a few fields
predictor_args.sampling_method_args = sampling_method_args;
predictor_args.tree_dir = [num2str(job_id) '/' zerostr(args.task_id,2), '_trees/'];
predictor_args.save_path = [asymmetryroot 'data/line_' args.detection_type '_rfs/' num2str(job_id) '/random_forest' zerostr(args.task_id,2) '.mat'];
predictor_args.tree_root = [asymmetryroot 'data/line_' args.detection_type '_rfs/'];

% set random seed based on task id and clock to ensure unique trees created
if isempty(args.rand_seed);
    rand('twister', sum(args.task_id*clock)); %#ok
	randn('state', sum(args.task_id*clock));
else
    rand('twister', args.rand_seed);
	randn('state', args.rand_seed);
end

%Clear the args structure and display the sampling/forest args
clear args;
display(sampling_method_args);
display(predictor_args);

%build forest
switch sampling_method_args.detection_type
	case 'detection'
		predictor = mb_random_forest_class_train(predictor_args);
	case 'orientation'
		predictor = mb_random_forest_reg_train(predictor_args);    
	case 'width'
		predictor = mb_random_forest_reg_train(predictor_args);
	case 'linear_regression'
		predictor_args.tree_dir = [num2str(job_id) '/'];
		predictor.n_trees = 1;
		predictor = linear_regressor_train(predictor_args);
	case 'boosted_regression'
		predictor_args.tree_dir = [num2str(job_id) '/'];
		predictor.n_trees = 1;
		predictor = boosted_regressor_train(predictor_args);
	otherwise
		display(['Detection type: ' sampling_method_args.detection_type ' not recognised, using ''detection''']);
				% shouldn't this be forest_args.sampling_method_args... = 'detection'?
		sampling_method_args.detection_type = 'detection';
		predictor_args.sampling_method_args = sampling_method_args;
		predictor = mb_random_forest_class_train(predictor_args);
end
	
%For the first task of job, save the sampling args
args_name = [asymmetryroot 'data/line_' sampling_method_args.detection_type '_rfs/' num2str(job_id) '/sampling_args.mat'];
if ~exist(args_name, 'file');
    save(args_name, 'sampling_method_args');
end

%For the first task of job, dump the username, forest arguments and
%sampling arguments to a text file.
filename = [asymmetryroot 'data/line_' sampling_method_args.detection_type '_rfs/' num2str(job_id) '/args.txt'];
if ~exist(filename, 'file');
	fid = fopen(filename,'w');
	fprintf(fid,'User = %s\n',get_username());
	fprintf(fid,'%s\n',evalc('predictor_args'));
	fprintf(fid,'%s\n',evalc('sampling_method_args'));
	fclose(fid);
end

display('Predictor successfully constructed!');
warning('on', 'ASYM:unexpectedArgument');
