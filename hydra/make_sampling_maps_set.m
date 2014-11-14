function [] = make_sampling_maps_set(varargin)
%MAKE_SAMPLING_MAPS_SET wrapper function to compute resampling maps to
%weight the probability of sampling particular pixels in training data
%
%
% Inputs:
%
%
% Outputs:
%
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
    {'model_id'}, ...
    'error_ratio',           unixenv('MAKE_RESAMPLING_MAPS', 0.5), ...
    'output_type',      unixenv('OUTPUT_TYPE', 'width'), ...
	'prediction_type',	unixenv('PREDICTION_TYPE','rf_regression'), ...
    'fg_mask_dir',      unixenv('FG_MASK_DIR', 'vessel_masks'),...
    'prediction_dir',   unixenv('PREDICTION_DIR', 'predictions'),...
    'probability_dir',  unixenv('PROBABILITY_DIR', 'resampling_maps'),...    
    'image_root',       [unixenv('DATA_ROOT',[]) unixenv('IMAGE_ROOT', [asymmetryroot,'data/retinograms/DRIVE/training'])],...
    'fov_mask_dir',		unixenv('FOV_MASK_DIR',[]),...
    'ori_dir',          unixenv('ORI_DIR', 'orientations'),...
    'width_dir',        unixenv('WIDTH_DIR', 'width_maps'), ...
    'view',				unixenv('VIEW',[]) ...
);


%Construct paths to masks
fg_mask_dir = [args.image_root '/' args.fg_mask_dir '/'];
fg_mask_dir = prettypath(fg_mask_dir);

%Construct paths to prediction images
prediction_dir = [args.image_root '/' args.prediction_dir '/' args.output_type '/' ...
                      args.prediction_type '/' args.model_id  '/'];
prediction_dir = prettypath(prediction_dir);

%Construct paths to output resampling maps and create dir
resampling_maps_dir = [args.image_root '/' args.probability_dir '/' args.output_type '/' ...
                      args.prediction_type '/' args.model_id  '/' num2str(args.error_ratio) '/'];
resampling_maps_dir = prettypath(resampling_maps_dir);                 
create_folder(resampling_maps_dir);

display(fg_mask_dir);
display(prediction_dir);
display(resampling_maps_dir);

%Get list of prediction images
prediction_list = dir([prediction_dir '*' args.view '*.mat']);

%Get list of FG masks
fg_mask_list = dir([fg_mask_dir '*' args.view '*.mat']);

%Check whether we need a list of masks
if ~isempty(args.fov_mask_dir)
    fov_mask_dir = [args.image_root '/' args.fov_mask_dir '/'];
    fov_mask_list = dir([fov_mask_dir '*' args.view '*.mat']);
end

%Get lists of GT maps if necessary
switch args.output_type
    case {'orientation', 'centre_orientation'}
        gt_dir = [args.image_root '/' args.ori_dir '/'];
        gt_list = dir([gt_dir '*' args.view '*.mat']);
        
    case 'width'
        gt_dir = [args.image_root '/' args.width_dir '/'];
        gt_list = dir([gt_dir '*' args.view '*.mat']);
end

%3. Predict output for each each image

for ii = 1:length(prediction_list)
    display(['Generating map for image ' num2str(ii)]);
    
    %Load preiction image
    prediction_image = u_load([prediction_dir prediction_list(ii).name]);
    
    %Load fg_mask
    fg_mask = u_load([fg_mask_dir fg_mask_list(ii).name]);
    
    %Load mask if necessary
    if ~isempty(args.fov_mask_dir)
        fov_mask = u_load([fov_mask_dir fov_mask_list(ii).name]);
    end
    
    %Load in GT if necessary
    if exist('gt_list', 'var')
        gt_map = u_load([gt_dir gt_list(ii).name]);
    else
        gt_map = [];
    end

    %Do prediction
    [resampling_map] = make_resampling_map(...
        prediction_image, fg_mask, fov_mask, gt_map, args.error_ratio, args.output_type); %#ok

    save([resampling_maps_dir prediction_list(ii).name], 'resampling_map');
            
end
