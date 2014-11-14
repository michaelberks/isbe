function [roc_pts auc tp_count fp_count] = calculate_roc_curve3D(class_vals1, class_vals2,class_labels,operating_pts)
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

%make sure class_labels are logicals
if ~islogical(class_labels)
    class_labels = logical(class_labels);
end

%compute number of positive/negative labels
n_pos = sum(class_labels);
n_neg = length(class_labels) - n_pos;

%pre-allocate space for roc_pts
n_pts = length(operating_pts);
tp_count = zeros(n_pts,n_pts);
fp_count = zeros(n_pts,n_pts);
roc_pts = zeros(n_pts, 2, n_pts);
auc = zeros(n_pts, 1);

%Sort thresholds from high to low
operating_pts = sort(operating_pts, 'descend');

% Calcultae ROC points at thresholds from 0 to 5
for ii = 1:n_pts
    
    th1 = operating_pts(ii);
    for jj = 1:n_pts
    
        th2 = operating_pts(ii);

        %Compute true positives, and true positive fraction
        tp_count(jj,ii) = sum(class_vals1(class_labels) > th1 & class_vals2(class_labels) > th2); %also FN = n_pos - TP;


        %Compute false positives, and false positive fraction
        fp_count(jj,ii) = sum(class_vals1(~class_labels) > th1 & class_vals2(~class_labels) > th2); %also TN = n_neg - FP;
    
    end
    
    %Set of points on ROC curve are the true positive fractions vs flase
    %positive fractions (see commented out code above)
    roc_pts(:, :, ii) = [fp_count(:,ii)/n_neg tp_count(:,ii)/n_pos];
    
    auc(ii) = sum( (roc_pts(2:end,1,ii)-roc_pts(1:end-1,1,ii)) .* roc_pts(1:end-1,2,ii)) + ...
            0.5*sum( (roc_pts(2:end,1,ii)-roc_pts(1:end-1,1,ii)) .* (roc_pts(2:end,2,ii)-roc_pts(1:end-1,2,ii)) );

end
  
% q2 = sum( ...
%     (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* ...
%     (roc_pts(1:end-1,2).^2 + roc_pts(1:end-1,2).*(roc_pts(2:end,2)-roc_pts(1:end-1,2)) +...
%     (roc_pts(2:end,2)-roc_pts(1:end-1,2)).*(roc_pts(2:end,2)-roc_pts(1:end-1,2))/3) );
% 
% q1 = sum( ...
%     (roc_pts(2:end,2)-roc_pts(1:end-1,2)) .* ...
%     ((1-roc_pts(2:end,1)).^2 + (1-roc_pts(2:end,1)).*(roc_pts(2:end,1)-roc_pts(1:end-1,1)) +...
%     (roc_pts(2:end,1)-roc_pts(1:end-1,1)).*(roc_pts(2:end,1)-roc_pts(1:end-1,1))/3) );
% 
% auc_se = sqrt( (auc*(1-auc) + (n_pos-1)*(q1-auc^2) + (n_neg-1)*(q2-auc^2)) / (n_neg*n_pos) );