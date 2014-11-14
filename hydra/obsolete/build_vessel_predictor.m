function predictor = build_vessel_predictor(varargin)
%BUILD_VESSEL_PREDICTOR wrapper function to build a predictor for
%retinograms
%detector on hydra
%   [] = build_vessel_predictor(win_size, num_levels, do_max, feature_type, num_trees)
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
    .... % What we want to predict
    'prediction_type',		unixenv('PREDICTION_TYPE','detection'), ...
	... % Sampling parameters
    'num_samples',			unixenv('NUM_SAMPLES',2e3), ...
    'pts_per_image',		unixenv('PTS_PER_IMAGE',500), ...
    'sampling_method',      unixenv('SAMPLING_METHOD','generate_vessel_data'), ...
    'ret_dir',              unixenv('IMAGE_DIR', 'retinograms/DRIVE/training'),...
    'mask_dir',             unixenv('MASK_DIR', []),...
    'ori_dir',              unixenv('ORI_DIR', []),...
    'width_dir',            unixenv('WIDTH_DIR', 'width_maps'), ...
    'image_type', 			unixenv('IMAGE_TYPE','line'), ...
    'rgb_channel',          unixenv('RGB_CHANNEL','rgb'), ...
    'normalise', 			unixenv('NORMALISE',0), ...
	... % Image feature parameters
    'decomp_type', 			unixenv('DECOMP_TYPE','dt'), ...
    'num_levels', 			unixenv('NUM_LEVELS',5), ...
    'win_size',				unixenv('WIN_SIZE',3), ...
    'pca',					unixenv('PCA',[]), ...
    'do_max',				unixenv('DO_MAX',0), ...
    'rotate',				unixenv('ROTATE',0), ...
    ... % DTCWT parameters
    'feature_shape', 		unixenv('FEATURE_SHAPE','rect'), ...
    'feature_type',			unixenv('FEATURE_TYPE','conj'), ...
    ... % Gaussian derivative parameters
    'sigma_range', 			unixenv('SIGMA_RANGE',[1 2 4 8]), ...
    ... % LinOp parameters
    'num_angles', 			unixenv('NUM_ANGLES',8), ...
    ... % Monogenic parameters
    'min_wavelength',		unixenv('MIN_WAVELENGTH',4), ...
    'onf',					unixenv('ONF',0.65), ...
    ... % Tree parameters
    'n_trees',				unixenv('NUM_TREES',2), ...
    'split_criterion',      unixenv('SPLIT_CRITERION', 'gdi'),...
    'var_criterion',		unixenv('VAR_CRITERION', 'mabs'),...
    'split_min',			unixenv('SPLIT_MIN',10), ...
    'end_cut_min',			unixenv('END_CUT_MIN',1), ...
    'do_ubound',			unixenv('DO_UBOUND',1), ...
    'do_circular',			unixenv('DO_CIRCULAR',[]), ...
	'w_prior',				unixenv('W_PRIOR',0), ...
    'impure_thresh',		unixenv('IMPURE_THRESH',1e-4), ...
	'minimise_size',		unixenv('MINIMIZE_TREE',0), ...
    'd',                    unixenv('d',[]), ...
    ... % Miscellaneous parameters
    'overwrite',			unixenv('OVERWRITE',0), ...
    'use_nag',				unixenv('USE_NAG',1),...
    'rand_seed',			unixenv('RAND_SEED',[]));

% if $CUSTOM_ID is defined and not empty then use that instead of $JOB_ID
custom_id = unixenv('CUSTOM_ID');
if ~isempty(custom_id)
	args.job_id = custom_id;
end
display(args);

job_id = args.job_id;

if strcmp(get_username(),'ptresadern') && ispc
    args.use_nag = 0;
end

%% copy arguments for sampling

% Start with those common to every experiment
fields_to_copy = {...
    'image_type'; 'decomp_type'; 'prediction_type';
    'num_samples'; 'pca'; 'rgb_channel'; 'bg_ratio';
    'image_dir'; 
};

% Add those that depend on the image decomposition (inputs)
decomp_args = get_sampling_args_from(args);
fields_to_copy = [fields_to_copy(:); fieldnames(decomp_args)];

% Add those that depend on what we want to predict (outputs)
ret_dir = [asymmetryroot 'data/'  args.ret_dir];
switch args.prediction_type
    case {'detection'}
        args.bg_ratio = 1;

    case {'orientation', 'centre_orientation', ...
          'mixed_orientation', 'mixed_centre_orientation', ...
          'linear_regression', 'boosted_regression'}
        args.bg_ratio = 0;
        args.ori_dir = [ret_dir '/' args.ori_dir '/'];
        fields_to_copy = [fields_to_copy(:); 'ori_dir'];

    case {'width'}
        args.bg_ratio = 0;
        args.width_dir = [ret_dir '/' args.width_dir '/'];
        fields_to_copy = [fields_to_copy(:); 'width_dir'];

    otherwise
        error(['Unknown prediction type:', args.prediction_type]);
end

% Add those that depend on from where we sample data
switch args.sampling_method
    case {'generate_vessel_data', 'generate_vessel_data_pt'}
        args.image_dir = [ret_dir '/images_extended/'];
        args.foveal_mask_dir = [ret_dir '/foveal_masks/'];
        args.vessel_mask_dir = [ret_dir '/vessel_masks/'];
        
        fields_to_copy = [fields_to_copy(:);
            'foveal_mask_dir'; 'vessel_mask_dir'; 
        ];
        
    case {'resample_vessel_data'}
        args.image_dir = [ret_dir '/images/'];
        args.bg_sample_map_dir = [asymmetryroot 'data/'  args.mask_dir '/bg_maps/'];
        args.vessel_sample_map_dir = [asymmetryroot 'data/'  args.mask_dir '/vessel_maps/'];

        fields_to_copy = [fields_to_copy(:);
            'bg_sample_map_dir'; 'vessel_sample_map_dir'; 
        ];
        
    case {'resample_retinogram_data'}
        args.image_dir = [ret_dir '/images/'];
        args.sample_map_dir = [asymmetryroot 'data/'  args.mask_dir '/'];
        args.vessel_mask_dir = [ret_dir '/vessel_masks/'];
        
        fields_to_copy = [fields_to_copy(:);
            'sample_map_dir'; 'vessel_mask_dir';
        ];
    
    otherwise
        error(['Unknown sampling method: ', args.sampling_method]);
end

% copy marked fields from args to sampling_method_args
sampling_method_args = get_substructure(args, fields_to_copy);

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

predictor_dir = [asymmetryroot 'data/models/vessel/' args.prediction_type];
predictor_args.save_path = [predictor_dir '/' num2str(job_id) '/random_forest' zerostr(args.task_id,2) '.mat'];
predictor_args.tree_root = [predictor_dir '/'];

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
switch sampling_method_args.prediction_type
	case {'detection', 'centre_detection', ...
          'junction_detection', 'junction_centre_detection', ...
          'foveal_detection', 'foveal_centre_detection'}
		predictor = mb_random_forest_class_train(predictor_args);
        
	case {'orientation', 'centre_orientation', ...
          'mixed_orientation','mixed_centre_orientation',...
          'width'}
		predictor = mb_random_forest_reg_train(predictor_args);

    case {'linear_regression'}
		predictor_args.tree_dir = [num2str(job_id) '/'];
		predictor = linear_regressor_train(predictor_args);

    case 'boosted_regression'
		predictor_args.tree_dir = [num2str(job_id) '/'];
    	predictor = boosted_regressor_train(predictor_args);
        
    otherwise
        if strcmp(get_username(),'ptresadern')
            error(['Unknown predictor: ',sampling_method_args.prediction_type]);
        else
            display(['Detection type: ' sampling_method_args.prediction_type ' not recognised, using ''detection''']);
                    % shouldn't this be forest_args.sampling_method_args... = 'detection'?
            sampling_method_args.prediction_type = 'detection';
            predictor_args.sampling_method_args = sampling_method_args;
            predictor = mb_random_forest_class_train(predictor_args);
        end
end
	
% define the output path
outpath = [asymmetryroot 'data/models/vessel/' ...
           sampling_method_args.prediction_type '/' num2str(job_id)];

%For the first task of job, save the sampling args
args_name = [outpath '/sampling_args.mat'];
if ~exist(args_name, 'file');
    save(args_name, 'sampling_method_args');
end

%For the first task of job, dump the username, forest arguments and
%sampling arguments to a text file.
filename = [outpath '/args.txt'];
if ~exist(filename, 'file');
	fid = fopen(filename,'w');
	fprintf(fid,'User = %s\n',get_username());
	fprintf(fid,'%s\n',evalc('predictor_args'));
	fprintf(fid,'%s\n',evalc('sampling_method_args'));
	fclose(fid);
end

display('Predictor successfully constructed!');
warning('on', 'ASYM:unexpectedArgument');


