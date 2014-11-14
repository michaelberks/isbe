function [orientation_errors predicted_orientations error_stats gt_orientations] =...
    compute_image_orientation_errors(prediction_dir, fg_mask_dir, varargin)
%COMPUTE_IMAGE_ORIENTATION_ERRORS Compute ROC curve for a set of test images that have
%been classified
%   [orientation_errors] = compute_image_orientation_errors(label_dir,prediction_dir,label_name)
%
% Inputs:
%      image_dir- directory containing test images with known orientations
%
%      prediction_dir- directory containing orientation maps of the test images
%
%
% Outputs:
%      mean_errors
%
%      std_errors
%
%      total_samples
%
%
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 08-Feb-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
warning('off', 'load_uint8:missing_variables');
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'label_type',	[], ...
    'gt_orientations',      [],...
    'centre_only',  0,...
    'centre_idx',   [],...
    'label_dir',	[], ...
    'fov_mask_dir',	[],...
    'save_dir', []);

%get directory of preiction images
prediction_dir = prettypath([prediction_dir '/']);
prediction_list = dir([prediction_dir, '*.mat']);

%get directory of preiction images
fg_mask_dir = prettypath([fg_mask_dir '/']);
fg_mask_list = dir([fg_mask_dir, '*.mat']);

%Get list of masks
if ~isempty(args.fov_mask_dir)
    args.fov_mask_dir = prettypath([args.fov_mask_dir filesep]);
    fov_mask_list = dir([args.fov_mask_dir '*.mat']);
end

%Checki if we need to sample the GT oris first
gt_orientations = args.gt_orientations; args = rmfield(args, 'gt_orientations');
if isempty(gt_orientations)
    if isempty(args.label_dir)
        error('Either ground truth orientations or a label directory must be provided');
    end
    args.label_dir = prettypath([args.label_dir '/']);
    label_list = dir([args.label_dir, '*.mat']);
    
    %Check lists are same length
    if length(prediction_list) ~= length(label_list)
        error('Number of labels and prediction images differ');
    end
    
    sample_gt = 1;
else
    sample_gt = 0;
    if args.centre_only
        gt_orientations = gt_orientations(args.centre_idx);
    end
end  

%Pre-allocate orientations
predicted_orientations = [];

%For each image
for ii = 1:length(prediction_list)
    
    %Load prediction image
    predicted_ori = load_uint8([prediction_dir, prediction_list(ii).name]);
    
    %Check if not complex - if not convert to cos(2x) + i*sin(2x)
    if isreal(predicted_ori)
        predicted_ori = cos(2*predicted_ori) + 1i*sin(2*predicted_ori);
    end
    
    %Load fg mask
    fg_mask = u_load([fg_mask_dir, fg_mask_list(ii).name]);
    
    if args.centre_only
        fg_mask = bwmorph(fg_mask, 'thin', inf);
    end
    
    %Load FoV mask
    if ~isempty(args.fov_mask_dir)
        fov_mask = u_load([args.fov_mask_dir, fov_mask_list(ii).name]);
        fg_mask = fg_mask & fov_mask;
    end
    
    %Sample predicted orientations
    predicted_orientations = [predicted_orientations; predicted_ori(fg_mask)]; %#ok
    
    if sample_gt
        %Load image label and sample GT orientations
        gt_ori = load_uint8([args.label_dir label_list(ii).name]);
        
        gt_orientations = [gt_orientations; gt_ori(fg_mask)]; %#ok
    end
end

%Compute orientation errors
[orientation_errors, error_stats] = ori_error(...
    predicted_orientations, gt_orientations);

if ~isempty(args.save_dir)
    args.save_dir = prettypath([args.save_dir '\']);
    create_folder(args.save_dir);
    save([args.save_dir 'error_stats.mat'],...
        'orientation_errors', 'error_stats', 'predicted_orientations');
end
