function [froc_pts] = compute_froc_curve(tp_counts, fp_counts, fp_areas, prompt_area)
%COMPUTE_FROC_CURVE *Insert a one line summary here*
%   [froc_pts] = compute_froc_curve(tp_counts, fp_counts, fp_areas)
%
% Inputs:
%      tp_counts - *Insert description of input variable here*
%
%      fp_counts - *Insert description of input variable here*
%
%      fp_areas - *Insert description of input variable here*
%
%
% Outputs:
%      froc_pts - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 07-Jan-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

[num_images num_pts] =  size(fp_counts);

%For the FP counts, workout the threshold index with maximum count for each
%image
[max_fp max_fp_idx] = max(fp_counts, [], 2);

%Workout FP areas at the maximum indices
areas_at_max = fp_areas(num_pts*(0:num_images-1)' + max_fp_idx);

%For each image, swap the fp counts to the left of the maximum count with
%the increase in promt_area added to the maximum count (that's not really
%in english is it?)
fp_combined = fp_counts;
for ii = 1:num_images
    cols = 1:max_fp_idx(ii)-1;
    fp_combined(ii, cols) = max_fp(ii) + (fp_areas(ii,cols) - areas_at_max(ii)) / prompt_area;
end

%Now compute the FROC curve
froc_pts(:,1) = mean(fp_combined)'; %False positive per image
froc_pts(:,2) = mean(tp_counts)'; %Sensitivity
        
