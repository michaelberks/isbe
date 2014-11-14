function [] = predict_image_set(varargin)
%PREDICT_IMAGE_SET wrapper function to predict ouput for a set of images in
%a given directory
%
%
% Inputs:
%
%
% Outputs:
%
% Example: for use on hydra
%   NUM_JOBS=20 JOB_ID="'191658'" TEST_IMAGE_DIR="'lines512'" qsub -N class_im -t 1-20 -V matlab_code/trunk/hydra/predict_image_set.sh
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

warning('off', 'ASYM:unexpectedArgument');

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'model_id',         unixenv('MODEL_PATH', 'model'),...
    'image_dir',        unixenv('IMAGE_DIR', 'images'),...
    'prediction_dir',   unixenv('PREDICTION_DIR', 'predictions'),...
    'model_name',       unixenv('MODEL_NAME', 'predictor'),...
    'model_root',		[unixenv('DATA_ROOT',[]) unixenv('MODEL_ROOT',[asymmetryroot,'data/models/vessel'])], ...
    'image_root',       [unixenv('DATA_ROOT',[]) unixenv('IMAGE_ROOT', [asymmetryroot,'data/retinograms/DRIVE/test'])],...
    'task_id',			unixenv('SGE_TASK_ID',1), ...
    'num_jobs',			unixenv('NUM_JOBS',1), ...
    'mask_dir',			unixenv('MASK_DIR',[]), ...
    'max_size',			unixenv('MAX_SIZE',128), ...
    'num_trees',		unixenv('NUM_TREES_C',[]), ...
    'use_sampled_maps', unixenv('USE_SAMPLED_MAPS', 0),...
    'use_probs',		unixenv('USE_PROBS',0), ...
    'view',				unixenv('VIEW',[]), ...
    'resize',			unixenv('RESIZE',[]), ...
    'flip',				unixenv('FLIP',0),...
    'overwrite',		unixenv('OVERWRITE',false), ...
    'save_type',        unixenv('SAVE_TYPE','normal'));


%Construct paths to images and predictor
image_dir = [args.image_root '/' args.image_dir '/'];
image_dir = prettypath(image_dir);
model_dir = [args.model_root '/' args.model_id '/'];
model_dir = prettypath(model_dir);
predictor_name = [model_dir args.model_name '.mat'];


%Load in random forest and job arguments
predictor = u_load(predictor_name);
job_args = u_load([model_dir 'job_args.mat']);
prediction_type = job_args.predictor_args.prediction_type;
decomposition_args = job_args.decomposition_args;
output_type = job_args.sampling_args.output_type;
decomposition_args.decomp_type = check_decomp_type(decomposition_args.decomp_type);
clear job_args;

%Construct paths to results directory
results_dir = [args.image_root '/' args.prediction_dir '/' output_type '/' ...
                      prediction_type '/' args.model_id  '/'];
results_dir = prettypath(results_dir);
create_folder(results_dir);

display(image_dir);
display(results_dir);
display(predictor_name);
display(prediction_type);
display(decomposition_args);

%Get list of images
image_list = dir([image_dir '*' args.view '*.mat']);

%Check whether we need a list of masks
if ~isempty(args.mask_dir)
    mask_dir = [args.image_root '/' args.mask_dir '/'];
    mask_list = dir([mask_dir '*' args.view '*.mat']);
end

%Check whether we've been asked to use sampled maps (i.e. to predict
%training images)
if args.use_sampled_maps
    %We need to know the path to the lookup table that link image paths to
    %the sampled maps
    sampled_maps_lookup = u_load([predictor.sampled_data_dir 'sampled_maps_lookup.mat']);
end

%1. Workout number of images in job
job_size = ceil(length(image_list) / args.num_jobs);

%2. Workout start and end indices for job
start_i	= (args.task_id-1)*job_size + 1;
end_i	= min(args.task_id*job_size, length(image_list));
display(['Predicting images ' num2str(start_i) ' to ' num2str(end_i)]);

%3. Predict output for each each image

for i_im = start_i:end_i
    
    save_name = [results_dir image_list(i_im).name(1:end-4) '_pred.mat'];
    
    if ~args.overwrite && exist(save_name, 'file')
        display(['Image ' num2str(i_im) ...
            ' already exists. Skipping (switch overwrite to 1 if necessary)']);
        continue;
    end
    
    display(['Predicting image ' num2str(i_im)]);
    

    %Load image to predict
    test_image = u_load([image_dir image_list(i_im).name]);

    %Check if we need to resize
    if ~isempty(args.resize)
        test_image = imresize(test_image, args.resize, 'bilinear');
    end

    %Load mask if necessary
    if isempty(args.mask_dir)
        mask = [];
    else
        mask = u_load([mask_dir mask_list(i_im).name]);
        if any(size(mask) ~= [size(test_image,1) size(test_image,2)])
            mask = imresize(mask, size(test_image));
        end
    end
    
    %Load sampled map if necessary
    if ~args.use_sampled_maps
        tree_mask = [];
    else
        %lookup which sample map matches the image then load it
        idx = strcmp([image_dir image_list(i_im).name], sampled_maps_lookup(:,1));
        tree_mask = u_load(sampled_maps_lookup{idx,2});
    end

    %Check if we need to flip
    if args.flip
        if ~isempty(strfind(image_list(i_im).name, 'R'))
            test_image = fliplr(test_image);
            mask = fliplr(mask);
        end
    end

    %Do prediction
    [prediction_image] = predict_image(...
        'image_in', test_image, ...
        'predictor', predictor,...
        'decomposition_args', decomposition_args,...
        'prediction_type', prediction_type,...
        'output_type', output_type,...
        'mask', mask,...
        'tree_mask', tree_mask,...
        'num_trees', args.num_trees, ...
        'max_size', args.max_size,...
        'use_probs', args.use_probs);

    if args.flip
        if ~isempty(strfind(image_list(i_im).name, 'R'))
            prediction_image = fliplr(prediction_image);
        end
    end

    %Save the image in the output directory
    switch args.save_type
        case 'normal'
            save(save_name, 'prediction_image');
            
        case 'uint8'
            save_uint8(save_name, prediction_image);
            
        otherwise
            display(['Save method "' args.save_method '" not recognised, using normal save']);
            save(save_name, 'prediction_image');
    end
            
end

display('All images successfully predicted');

%Copy the job args txt file into the results dir
if ~exist([results_dir 'args.txt'], 'file')
    copyfile([model_dir 'args.txt'], [results_dir 'args.txt']);
end
