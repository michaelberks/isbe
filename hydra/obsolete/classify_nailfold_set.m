function [] = classify_nailfold_set(forest_job, image_dir, varargin)
%
%
%
% Inputs:
%
%
% Outputs:
%
% Example: for use on hydra
%   NUM_JOBS=20 FOREST_JOB="'191658'" TEST_IMAGE_DIR="'lines512'" qsub -N class_im -t 1-20 -V matlab_code/trunk/hydra/classify_image_set.sh
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

% this now uses the utils/unixenv function that takes the key name of the
% environment variable and returns its corresponding value, returning the
% default value if the environment variable does not exist
%   e.g. envval = unixenv('ENVVARNAME',default_val);

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'task_id',			unixenv('SGE_TASK_ID',1), ...
    'num_jobs',			unixenv('NUM_JOBS',1), ...
    'forest_dir',		unixenv('FOREST_DIR','line_detection_rfs'), ...
    'mask_dir',			unixenv('MASK_DIR',[]), ...
    'max_size',			unixenv('MAX_SIZE',128), ...
    'num_trees',		unixenv('NUM_TREES_C',[]), ...
    'use_probs',		unixenv('USE_PROBS',0), ...
    'use_nag',			unixenv('USE_NAG',1), ...
    'view',				unixenv('VIEW',[]), ...
    'resize',			unixenv('RESIZE',[]), ...
    'flip',				unixenv('FLIP',0));


%Construct paths to images and results directories
image_dir = [nailfoldroot image_dir '/'];
results_dir = strrep([image_dir 'results/' forest_job '/'],'\','/');
results_dir = strrep(results_dir,'//','/');
predictor_name = [nailfoldroot args.forest_dir '/' forest_job '/random_forest.mat'];

display(args);
display(image_dir);
display(results_dir);
display(predictor_name);

%Load in random forest and samplign arguments
predictor = u_load(predictor_name);
sampling_args = u_load([nailfoldroot args.forest_dir '/' forest_job '/sampling_args.mat']);

% set prediction_type for older predictors that use detection_type
if ~isfield(sampling_args,'prediction_type')
	sampling_args.prediction_type = sampling_args.detection_type;
end

%Check the type of forest we have
switch sampling_args.prediction_type
	case {'detection', 'classification', 'centre_detection', 'junction_detection',...
            'junction_centre_detection', 'foveal_detection', 'foveal_centre_detection'}
		forest_type = 'isbe';
        predictor.tree_root = [nailfoldroot args.forest_dir '/'];
	case {'orientation', 'centre_orientation','mixed_orientation','mixed_centre_orientation'}
		forest_type = 'orientation';
        predictor.tree_root = [nailfoldroot args.forest_dir '/'];
	case 'width'
		forest_type = 'regression';
        predictor.tree_root = [nailfoldroot args.forest_dir '/'];
	case 'linear_regression'
		forest_type = 'linear_regression';
	case 'boosted_regression'
		forest_type = 'boosted_regression';
	otherwise
		display('Detection type not recognised, forest_type set to isbe');
		forest_type = 'isbe';
        predictor.tree_root = [nailfoldroot args.forest_dir '/'];
end

%Get list of images
image_list = dir([image_dir '*' args.view '*.mat']);

%Check whether we need a list of masks
if ~isempty(args.mask_dir)
    mask_dir = [nailfoldroot args.mask_dir '/'];
    mask_list = dir([mask_dir '*' args.view '*.mat']);
end

%1. Workout number of images in job
job_size = ceil(length(image_list) / args.num_jobs);

%2. Workout start and end indices for job
start_i	= (args.task_id-1)*job_size + 1;
end_i	= min(args.task_id*job_size, length(image_list));

display(['classifying images ' num2str(start_i) ' to ' num2str(end_i)]);

%Check for decomposition type
if isfield(sampling_args, 'decomp_type')
    decomp_type = sampling_args.decomp_type;
else
    decomp_type = 'dt';
end
switch decomp_type
    case 'dt'
        sampling_args_c.num_levels = sampling_args.num_levels;
        sampling_args_c.feature_shape = sampling_args.feature_shape;
        sampling_args_c.feature_type = sampling_args.feature_type;
        sampling_args_c.do_max = sampling_args.do_max;
        sampling_args_c.rotate = sampling_args.rotate;
        sampling_args_c.win_size = sampling_args.win_size;
        if isfield(sampling_args, 'use_nag')
            sampling_args_c.use_nag = sampling_args.use_nag;
        else
            sampling_args_c.use_nag = args.use_nag;
        end
    case 'mono'
        sampling_args_c.num_levels = sampling_args.num_levels;
        sampling_args_c.win_size = sampling_args.win_size;
        sampling_args_c.min_wavelength = sampling_args.min_wavelength;
        sampling_args_c.onf = sampling_args.onf;
    case {'g2d', 'g2di', 'clover', 'haar', 'g1d', 'g12d'}
        sampling_args_c.win_size = sampling_args.win_size;
        sampling_args_c.sigma_range = sampling_args.sigma_range;
    case 'linop'
        sampling_args_c.win_size = sampling_args.win_size;
        sampling_args_c.num_levels = sampling_args.num_levels;
        sampling_args_c.num_angles = sampling_args.num_angles;
        sampling_args_c.do_max = sampling_args.do_max;
        sampling_args_c.rotate = sampling_args.rotate;
    case 'pixel'
        sampling_args_c.win_size = sampling_args.win_size;
    case 'dtg2'
        sampling_args_c.num_levels = sampling_args.num_levels;
        sampling_args_c.feature_shape = sampling_args.feature_shape;
        sampling_args_c.feature_type = sampling_args.feature_type;
        sampling_args_c.do_max = sampling_args.do_max;
        sampling_args_c.rotate = sampling_args.rotate;
        sampling_args_c.win_size = sampling_args.win_size;
        sampling_args_c.sigma_range = sampling_args.sigma_range;
        if isfield(sampling_args, 'use_nag')
            sampling_args_c.use_nag = sampling_args.use_nag;
        else
            sampling_args_c.use_nag = args.use_nag;
        end
end
if isfield(sampling_args, 'pca') && ~isempty(sampling_args.pca);
    %pca should have fields mean and modes
    sampling_args_c.pca = u_load(sampling_args.pca);
end
if isfield(sampling_args, 'normalise')
    normalise = sampling_args.normalise;
else
    normalise = 0;
end
classify_fcn = @classify_image;

display(sampling_args_c);
display(['classifying images ' num2str(start_i) ' to ' num2str(end_i)]);

%3. Classify each image using forest
if isempty(args.num_trees)
    num_num_trees = 1;
    num_trees = [];
else
    num_num_trees = length(args.num_trees);
end

for jj = 1:num_num_trees
	% Make top-level output folder
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
		if ~ispc,
			fileattrib(results_dir,'+w','g');
		end
    end
    if isempty(args.num_trees)
		% Put results in top-level folder
		results_dir_jj = results_dir;
    else
		% Create a subfolder and put results there
        num_trees = args.num_trees(jj);
        results_dir_jj = [results_dir zerostr(num_trees,3) '/'];
        if ~exist(results_dir_jj, 'dir')
			mkdir(results_dir_jj);
            if ~ispc
				fileattrib(results_dir_jj,'+w','g');
            end
        end
    end
	
    for ii = start_i:end_i
        display(['classifying image ' num2str(ii)]);

        %Load nailfold to classify
        test_image = u_load([image_dir image_list(ii).name]);
    
        %Check if we need to resize
        if ~isempty(args.resize)
            test_image = imresize(test_image, args.resize, 'bilinear');
        end
        
        %Load mask if necessary
        if isempty(args.mask_dir)
            mask = [];
        else
            mask = u_load([mask_dir mask_list(ii).name]);
            if any(size(mask) ~= [size(test_image,1) size(test_image,2)])
                mask = imresize(mask, size(test_image));
            end
        end
            
        %Check if we need to flip
        if args.flip
            if ~isempty(strfind(image_list(ii).name, 'R'))
                test_image = fliplr(test_image);
                mask = fliplr(mask);
            end
        end

        %Check if we need to normalise the image
        if normalise
            im_mean = mean(test_image(:));
            im_std = std(test_image(:));
            test_image = (test_image - im_mean) / im_std;
        end

        [prediction_image] = classify_fcn(...
            'image_in', test_image, ...
            'forest', predictor,...
            'sampling_args', sampling_args_c,...
            'forest_type', forest_type,...
            'decomp_type', decomp_type,...
            'mask', mask,...
            'num_trees', num_trees, ...
            'max_size', args.max_size,...
            'use_probs', args.use_probs);

        if args.flip
            if ~isempty(strfind(image_list(ii).name, 'R'))
                prediction_image = fliplr(prediction_image); %#ok
            end
        end

        %Save the image in the output directory
        save_uint8([results_dir_jj image_list(ii).name(1:end-4) '_class.mat'], prediction_image);
    end
end
display('All images successfully classified');
