function [width_errors predicted_widths gt_widths] =...
    compute_image_width_errors(im_names, prediction_dir, fg_mask_dir, varargin)
%COMPUTE_IMAGE_ORIENTATION_ERRORS Compute ROC curve for a set of test images that have
%been classified
%   [orientation_errors] = compute_image_orientation_errors(label_dir,prediction_dir,label_name)
%
% Inputs:
%      image_dir- directory containing test images with known widths
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
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'label_type',	[], ...
    'gt_widths',      [],...
    'centre_idx',   [],...
    'label_dir',	[], ...
    'fov_mask_dir',	[],...
    'pred_ext',     '_pred.mat',...
    'label_ext',    '_width.mat',...
    'fg_mask_ext',  '_v_mask.mat',...
    'fov_mask_ext', '_f_mask.mat',...
    'do_log',       0);

%get directory of preiction images
prediction_dir = prettypath([prediction_dir '/']);

%get directory of preiction images
fg_mask_dir = prettypath([fg_mask_dir '/']);

%Get list of masks
if ~isempty(args.fov_mask_dir)
    args.fov_mask_dir = prettypath([args.fov_mask_dir filesep]);
end

%Checki if we need to sample the GT oris first
gt_widths = args.gt_widths; args = rmfield(args, 'gt_widths');
if isempty(gt_widths)
    if isempty(args.label_dir)
        error('Either ground truth widths or a label directory must be provided');
    end
    args.label_dir = prettypath([args.label_dir '/']);    
    sample_gt = 1;
else
    sample_gt = 0;
end  

%Pre-allocate widths
predicted_widths = [];

%For each image
for i_im = 1:length(im_names)
    
    %Load prediction image
    predicted_width = load_uint8([prediction_dir, im_names{i_im} args.pred_ext]);
    
    %Load fg mask
    fg_mask = load_uint8([fg_mask_dir, im_names{i_im} args.fg_mask_ext]);
    fg_mask = bwmorph(fg_mask, 'thin', inf);
    
    %Load FoV mask
    if ~isempty(args.fov_mask_dir)
        fov_mask = load_uint8([args.fov_mask_dir, im_names{i_im} args.fov_mask_ext]);
        fg_mask = fg_mask & fov_mask;
    end
    
    %Sample predicted widths
    predicted_widths = [predicted_widths; predicted_width(fg_mask)]; %#ok
    
    if sample_gt
        %Load image label and sample GT widths
        gt_width = load_uint8([args.label_dir im_names{i_im} args.label_ext]);
        
        gt_widths = [gt_widths; gt_width(fg_mask)]; %#ok
    end
end

if args.do_log
    predicted_widths = exp(predicted_widths);
end

%Compute orientation errors
width_errors = predicted_widths - gt_widths; 
      
