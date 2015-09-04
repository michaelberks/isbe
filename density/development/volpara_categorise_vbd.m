function [vbd_cat4, vbd_cat5, vbd_mean, vbd_lr_max, outliers] = volpara_categorise_vbd(vbd_scores, view_flags, varargin)
%VOLPARA_CATEGORISE_VBD *Insert a one line summary here*
%   [vbd_cat4, vbd_cat5, vbd_mean, vbd_lr_max] = volpara_categorise_vbd(vbd_scores, view_flags, varargin)
%
% Inputs:
%      vbd_scores - 4 element vector of VBD scores, *MUST* be in order LCC,
%      LMLO, RCC, RMLO
%
%      view_flags - 4 element logical vector, specifying which views are
%      present (1) or missing (0)
%
%      varargin - Options supplied in name/value pairs:     
%       cat4_thresh: Thresholds for BIRADS 4th edition categories - default [4.5 7.5 15.5 inf]     
%     	cat5_thresh: Thresholds for BIRADS 5th edition categories - default [3.5 7.5 15.5 inf]
%     	do_outlier_rejection: default true;
%     	min_to_max_ssd_thresh: Outlier rejection: Min SSD to max SSD ratio threshold - default = 0.4^2;
%     	median_vbd_scale: Outlier rejection, scale factor applied to median to compute min SSD threshold - default 3
%     	median_vbd_offset: Outlier rejection, offset applied to median to compute min SSD threshold - default 0
%
% Outputs:
%      vbd_cat4 - VBD category group for BIRADS 4 thresholds (uses vbd_mean)
%
%      vbd_cat5 - VBD category group for BIRADS 5 thresholds (uses vbd_lr_max) 
%
%      vbd_mean - Mean average of all present views
%
%      vbd_lr_max - Max of the left/right mean averages
%
%
% Example:
%
% Notes: This function replicates the behaviour of the MS Excel macro
% 
%
% See also: VOLPARA_CATEGORISE_VBD_BATCH
%
% Created: 04-Sep-2015
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 5285 
% Copyright: (C) University of Manchester 

%Check required input types so we can bug out early if incorrect format
if length(vbd_scores) ~= 4 || ~isnumeric(vbd_scores)
    error('vbd_scores must be a numeric vector with 4 elements');
end
%Allow view flags to be missing so can call function with just vbd_scores
if ~exist('view_flags', 'var') || isempty(view_flags)
    view_flags = true(4,1);
    
elseif length(view_flags) ~= 4 || ~islogical(view_flags)
    error('view flags must be a logical vector with 4 elements');
end
vbd_scores = vbd_scores(:);
view_flags = view_flags(:);


%Check if any optional arguments have been set
if isempty(varargin)
    %Use default category thresholds and outlier rejection parameters
    args = set_default_args();
else
    args = process_args(varargin);
end

%Crack on with the function...

%Outlier rejection
if args.do_outlier_rejection && sum(view_flags)==4
    outliers = detect_outliers(vbd_scores,...
        args.min_to_max_ssd_thresh, args.median_vbd_scale, args.median_vbd_offset);
    view_flags(outliers) = 0;
else
    outliers = false(4,1);
end

%Compute vbd mean and max of left/right means
view_flags_l = [view_flags([1 2]); false; false];
view_flags_r = [false; false; view_flags([3 4])];
vbd_mean = mean(vbd_scores(view_flags));
vbd_lr_max = max(mean(vbd_scores(view_flags_l)), mean(vbd_scores(view_flags_r)));

%Get BIRADS 4th edition cat using vbd_mean
vbd_cat4 = 1;
while (vbd_mean > args.cat4_thresh(vbd_cat4))
    vbd_cat4 = vbd_cat4 + 1;
end

%Get BIRADS 5th edition cat using vbd_lr_max
vbd_cat5 = 1;
while (vbd_lr_max > args.cat5_thresh(vbd_cat5))
    vbd_cat5 = vbd_cat5 + 1;
end

%All done!
return;


%--------------------------------------------------------------------------
function args = set_default_args()
    args.cat4_thresh = [4.5 7.5 15.5 inf];
    args.cat5_thresh = [3.5 7.5 15.5 inf];
    args.do_outlier_rejection = true;
    args.min_to_max_ssd_thresh = 0.4^2;  %Outlier rejection: Min SSD to max SSD ratio threshold
    args.median_vbd_scale = 3; %Outlier rejection: Median density scale factor
    args.median_vbd_offset = 0; %Outlier rejection: Median density offset

%--------------------------------------------------------------------------
function args = process_args(args_in)
%Note I'd usually just include a call to u_parkargs at the start of the
%main function, but assuming the majority use of this function will be to
%run 10,000's times using defaults, so more efficient to only call if the
%user has actually set some options. The downside is it means manually
%maintaining defaults in set_default_args and below.
args = u_packargs(args_in,... % the user's input
     '0', ... % non-strict mode
     'cat4_thresh', [4.5 7.5 15.5 inf],... % the optional arguments
     'cat5_thresh', [3.5 7.5 15.5 inf],...
     'do_outlier_rejection', true,...
     'min_to_max_ssd_thresh', 0.4^2,...
     'median_vbd_scale', 3, ...
     'median_vbd_offset', 0);

%--------------------------------------------------------------------------
function outliers = detect_outliers(vbd_scores, min_to_max_ssd_thresh, median_vbd_scale, median_vbd_offset)

%Compute sum of square differences between remaining 3 views, after leaving
%each view out it turn
sum_sq_diffs = zeros(4,1);
for i_view = 1:4  
    vbd_scores_i = vbd_scores([1:(i_view-1) (i_view+1):4]);
    sum_sq_diffs(i_view) = (...
        (vbd_scores_i(1) - vbd_scores_i(2))^2 +...
        (vbd_scores_i(3) - vbd_scores_i(2))^2 +...
        (vbd_scores_i(1) - vbd_scores_i(3))^2 ) / 3;
end

%Use median VBD to compute threshold for this case
threshold_ssd = median(vbd_scores)*median_vbd_scale + median_vbd_offset;

%Get min/max SSDs
min_ssd = min(sum_sq_diffs);
max_ssd = max(sum_sq_diffs);

%Check if outlier according Volpara's criteria
if (max_ssd > threshold_ssd) && ... %outlier_far_enough
    (min_ssd / max_ssd) < min_to_max_ssd_thresh %min_max_spread
        
    %Outlier is the view that returned the min SSD - note that computing it
    %as below permits the possibility of 2 views being listed as outliers
    %as per the code in the MACRO, although I suspect an proper analysis of
    %the rejection criteria woould show such a scenario is possible
    outliers = min_ssd == sum_sq_diffs;
else
    outliers = false(4,1);
end
    
        
