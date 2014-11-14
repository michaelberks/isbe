function [roc_pts auc tp_count fp_count auc_se] = calculate_roc_nms(prob_map, prob_ori, pos_label, neg_label, operating_pts, if_plot)
%CALCULATE_ROC_CURVE computes the points on an ROC curve for a set of input
% values and true class labels
%   [roc_pts] = calculate_roc_curve(class_vals,class_labels,operating_points)
%
% Inputs:
%      class_vals- 1D vector of values assigned to each point (usually a probability, but
%      not necessarily so)
%
%      class_labels - 1D vector of class labels, either {1,0} 
%
%      operating_pts- Threshold values at which to calculate ROC points,
%      for probabilities these will vary from 0 to 1, for other input
%      values they should vary from min(class_vals) to max(class_vals);
%
%
% Outputs:
%      roc_pts- nx2 vector of points that make up a the ROC curve
%      auc - area under the curve
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 06-Oct-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
if nargin < 5
    if_plot = 0;
end

%compute number of positive/negative labels
n_pos = sum(pos_label(:));
n_neg = sum(neg_label(:));

%now dilate the positive label

%pre-allocate space for roc_pts
n_pts = length(operating_pts);
tp_count = zeros(n_pts,1);
fp_count = zeros(n_pts,1);

%Sort thresholds from high to low
operating_pts = sort(operating_pts, 'descend');

% Calcultae ROC points at thresholds from 0 to 5
for ii = 1:n_pts
    
    th = operating_pts(ii);
    
    %threshold probability
    thresh_map = prob_map > th;
    
    if if_dilate
        %Now skeletonise then dilate this map
        thresh_map = bwmorph(thresh_map, 'thin', 2);
        %thresh_map = imdilate(thresh_map, strel('disk', 1));
    end
    
    if if_plot && (~mod(ii,5) || ii == n_pts-1)
        figure; image(make_tp_fp_image(pos_label, neg_label, thresh_map)); axis image;
    end
    
    %Compute true positives, and true positive fraction
    %TP = sum(class_vals(class_labels) > th); %also FN = n_pos - TP;
    %TPF = TP / n_pos; %since TPF = TP/(TP+FN) = TP/n_pos;
    
    tp_count(ii) = sum(thresh_map(:) & pos_label(:)); %also FN = n_pos - TP;


    %Compute false positives, and false positive fraction
    %FP = sum(class_vals(~class_labels) > th); %also TN = n_neg - FP; 
    %FPF = FP / n_neg; %since FPF = FP/(FP+TN) = FP/n_neg;

    fp_count(ii) = sum(thresh_map(:) & neg_label(:)); %also TN = n_neg - FP;
    
    %roc_pts(ii,1) = FPF;
    %roc_pts(ii,2) = TPF;

end

%Set of points on ROC curve are the true positive fractions vs flase
%positive fractions (see commented out code above)
roc_pts = [fp_count/n_neg tp_count/n_pos];

auc = sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* roc_pts(1:end-1,2)) + ...
        0.5*sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* (roc_pts(2:end,2)-roc_pts(1:end-1,2)) );
    
q2 = sum( ...
    (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* ...
    (roc_pts(1:end-1,2).^2 + roc_pts(1:end-1,2).*(roc_pts(2:end,2)-roc_pts(1:end-1,2)) +...
    (roc_pts(2:end,2)-roc_pts(1:end-1,2)).*(roc_pts(2:end,2)-roc_pts(1:end-1,2))/3) );

q1 = sum( ...
    (roc_pts(2:end,2)-roc_pts(1:end-1,2)) .* ...
    ((1-roc_pts(2:end,1)).^2 + (1-roc_pts(2:end,1)).*(roc_pts(2:end,1)-roc_pts(1:end-1,1)) +...
    (roc_pts(2:end,1)-roc_pts(1:end-1,1)).*(roc_pts(2:end,1)-roc_pts(1:end-1,1))/3) );

auc_se = sqrt( (auc*(1-auc) + (n_pos-1)*(q1-auc^2) + (n_neg-1)*(q2-auc^2)) / (n_neg*n_pos) );