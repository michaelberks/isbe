function [roc_pts auc tp_count fp_count auc_se] = calculate_roc_connected(prob_map,label_map,operating_pts, ignore_mask, varargin)
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
args = u_packargs(varargin, '0', ...
    'dilate', 0,...
    'thin', 0,...
    'plot', 0);

%compute number of positive/negative labels
pos_label = label_map & ignore_mask;
neg_label = ~label_map & ignore_mask;
n_pos = sum(pos_label(:));
n_neg = sum(neg_label(:));

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
    
    %Get largest connected component
    CC = bwconncomp(thresh_map);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,max_comp] = max(numPixels);
    
    connected_map = false(size(thresh_map));
    connected_map(CC.PixelIdxList{max_comp}) = 1;
    
    if args.thin
        %Thin the map
        connected_map = bwmorph(connected_map, 'thin', args.thin);
    end
    if args.dilate
        %Dilate the map
        connected_map = imdilate(connected_map, strel('disk', args.dilate));
    end
    
    if args.plot && (~mod(ii,5) || ii == n_pts-1)
        figure; image(make_tp_fp_image(pos_label, neg_label, connected_map)); axis image;
    end
    
    %Compute true positives, and true positive fraction
    %TP = sum(class_vals(class_labels) > th); %also FN = n_pos - TP;
    %TPF = TP / n_pos; %since TPF = TP/(TP+FN) = TP/n_pos;
    
    tp_count(ii) = sum(connected_map(:) & pos_label(:)); %also FN = n_pos - TP;


    %Compute false positives, and false positive fraction
    %FP = sum(class_vals(~class_labels) > th); %also TN = n_neg - FP; 
    %FPF = FP / n_neg; %since FPF = FP/(FP+TN) = FP/n_neg;

    fp_count(ii) = sum(connected_map(:) & neg_label(:)); %also TN = n_neg - FP;
    
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