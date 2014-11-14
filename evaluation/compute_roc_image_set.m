function [roc_pts, auc, auc_individual, tp_counts, fp_counts, t_counts, f_counts] = ...
    compute_roc_image_set(prob_dir, label_dir, varargin)
%COMPUTE_ROC_IMAGE_SET Compute ROC curve for a set of test images that have
%been classified
%   [roc_pts,auc,tp_counts,fp_counts] = compute_roc_image_set(image_dir,prob_dir,label_name)
%
% Inputs:
%      image_dir- directory containing (pre-classified) test images
%
%      prob_dir- directory containing probability maps of the test images
%
%      label_name- name of the variable containing truth labels for the
%       images (default: 'label_centre')
%
%
% Outputs:
%      roc_pts- N*2 array defining pts on ROC curve for data
%
%      auc- Area under the ROC curve
%
%      tp_counts- Total number of true positives for the data
%
%      fp_counts- Total number of false positives for the data
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

args = u_packargs(varargin, '0', ...
    'fov_mask_dir', [],...
    'centre_only', [],...
    'tol_mask_dir', [],...
    'dilate', 0,...
    'thin', 0,...
    'min_length', 0,...
    'plot', 0);


%Set up args for comput_roc_image
c_args.dilate = args.dilate;
c_args.thin = args.thin;
c_args.min_length = args.min_length;
c_args.plot = args.plot;

%get directory of probability images
prob_dir = prettypath([prob_dir filesep]);
prob_list = dir([prob_dir '*.mat']);

%get directory of image labels
label_dir = prettypath([label_dir filesep]);
label_list = dir([label_dir '*.mat']);

%Check lists are same length
if length(prob_list) ~= length(label_list)
    error('Number of image labels and probability images differ');
end

%Get list of masks
if ~isempty(args.fov_mask_dir)
    args.fov_mask_dir = prettypath([args.fov_mask_dir filesep]);
    fov_mask_list = dir([args.fov_mask_dir '*.mat']);
end

if ~isempty(args.tol_mask_dir)
    args.tol_mask_dir = prettypath([args.tol_mask_dir filesep]);
    tol_mask_list = dir([args.tol_mask_dir '*.mat']);
end

tp_counts = zeros(length(prob_list), 102);
fp_counts = zeros(length(prob_list), 102);
t_counts = zeros(length(prob_list), 1);
f_counts = zeros(length(prob_list), 1);
auc_individual = zeros(length(prob_list), 1);
%For each image
for ii = 1:length(prob_list)
    
    %Load probability image
    probability_image = load_uint8([prob_dir, prob_list(ii).name]);
    if ~isreal(probability_image)
        probability_image = abs(probability_image);
    end
        
    %Load image label
    label = load_uint8([label_dir label_list(ii).name]);
    
    
    
    %Load FoV mask
    if ~isempty(args.fov_mask_dir)
        fov_mask = load_uint8([args.fov_mask_dir fov_mask_list(ii).name]);
    else
        fov_mask = true(size(label));
    end
    
    if ~isempty(args.tol_mask_dir)
        c_args.tol_mask = load_uint8([args.tol_mask_dir tol_mask_list(ii).name]);
    end
    
    if args.centre_only
        label = bwmorph(label, 'thin', 'inf');
    end
    
    %Compute ROC counts for image
    [d d tp_count fp_count d t_count f_count] = calculate_roc_image(...
        probability_image, label, (-1:100)/100, fov_mask,...
        c_args);
    
    %Increment total counts
    tp_counts(ii,:) = tp_count;
    fp_counts(ii,:) = fp_count;
    t_counts(ii) = t_count;
    f_counts(ii) = f_count;
    
    [~,auc_individual(ii)] = ...
            do_roc(tp_count, fp_count, t_counts(ii), f_counts(ii));
    
end

%Compute ROC points for complete set of data
[roc_pts auc] = do_roc(sum(tp_counts), sum(fp_counts), sum(t_counts), sum(f_counts));

function [roc_pts auc] = do_roc(tp_count, fp_count, t_count, f_count)
    roc_pts = [fp_count(:) / f_count tp_count(:) / t_count];
    auc = sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* roc_pts(1:end-1,2)) + ...
        0.5*sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* (roc_pts(2:end,2)-roc_pts(1:end-1,2)) );
    

        
