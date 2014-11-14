function [roc_pts, auc, tp_counts, fp_counts, t_counts, f_counts] = compute_roc_image_set_lines(label_dir,prob_dir,label_type)
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
% Notes: This is a very specific file to compute ROC curves for our
% synthetic lines experiment. The label names it assumes are specific to
% this task and may not be suitable for more general image sets
%
% See also:
%
% Created: 08-April-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

if nargin < 3
    label_type = 'centre_line';
end

%get directory listing of input images
if ~isempty(label_dir) && ~strcmp(label_dir(end), filesep)
    label_dir = [label_dir filesep];
end
label_list = dir([label_dir, 'label*.mat']);

%get directory of probability images
if ~isempty(prob_dir) && ~strcmp(prob_dir(end), filesep)
    prob_dir = [prob_dir filesep];
end
prob_list = dir([prob_dir, '*.mat']);

%Check lists are same length
if length(prob_list) ~= length(label_list)
    error('Number of labels and probability images differ');
end

tp_counts = zeros(length(prob_list), 102);
fp_counts = zeros(length(prob_list), 102);
t_counts = zeros(length(prob_list), 1);
f_counts = zeros(length(prob_list), 1);

%For each image
for ii = 1:length(prob_list)
    
    %Load probability image
    probability_image = load_uint8([prob_dir, prob_list(ii).name]);
    
    %Load image label
    s = load([label_dir label_list(ii).name]); %#ok
    
    switch label_type
        case 'centre_line'
            label = s.label_centre(:);
            probability_image = probability_image(:);
            
        case 'all_line'
            label = s.label(:);
            probability_image = probability_image(:);
            
        case 'centre_not_all'
            keep = ~s.label | s.label_centre;
            label = s.label_centre(keep);
            probability_image = probability_image(keep);
        otherwise
            warning('Label type not recognised, using centre line'); %#ok
            label = s.label_centre(:);
            probability_image = probability_image(:);
            
    end
    
    %Compute ROC counts for image
    [d d tp_count fp_count] = calculate_roc_curve(probability_image, label,(-1:100)/100);
    
    %Increment total counts
    tp_counts(ii,:) = tp_count;
    fp_counts(ii,:) = fp_count;
    t_counts(ii) = sum(label(:));
    f_counts(ii) = sum(~label(:));
    
end

%Compute ROC points for complete set of data
roc_pts = [sum(fp_counts)' / sum(f_counts) sum(tp_counts)' / sum(t_counts)];

%Compute AUC for ROC curve
auc = sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* roc_pts(1:end-1,2)) + ...
        0.5*sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* (roc_pts(2:end,2)-roc_pts(1:end-1,2)) );
