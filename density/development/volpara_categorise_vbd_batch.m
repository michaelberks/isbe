function [unique_case_ids, vbd_cat4, vbd_cat5, vbd_mean, vbd_lr_max] = ...
    volpara_categorise_vbd_batch(vbd_scores, case_ids, view_labels, varargin)
%VOLPARA_CATEGORISE_VBD_BATCH *Insert a one line summary here*
%   [unique_case_ids, vbd_cat4, vbd_cat5, vbd_mean, vbd_lr_max] = ...
%       volpara_categorise_vbd_batch(vbd_scores, case_ids, view_labels, varargin)
%
% Inputs:
%      vbd_scores - Nx1 array of VBD scores
%
%      case_ids - Nx1 cell array of case IDs
%
%      view_labels - Nx1 cell array of view labels (each entry must be one
%      of {LCC, LML, RCC, RML} or since I'm being nice, lower-case variants
%      or or L(R)MLO
%
%      varargin - Options supplied in name/value pairs:     
%       cat4_thresh: Thresholds for BIRADS 4th edition categories - default [4.5 7.5 15.5 inf]     
%     	cat5_thresh: Thresholds for BIRADS 5th edition categories - default [3.5 7.5 15.5 inf]
%     	do_outlier_rejection: default true;
%     	min_to_max_ssd_thresh: Outlier rejection: Min SSD to max SSD ratio threshold - default = 0.4^2;
%     	median_vbd_scale: Outlier rejection, scale factor applied to median to compute min SSD threshold - default 3
%     	median_vbd_offset: Outlier rejection, offset applied to median to compute min SSD threshold - default 0
%
%
% Outputs:
%      unique_case_ids - Mx1 array where M = length(unique(case_ids))
%
%      vbd_cat4 - Mx1 vector of BIRADS 4th edition categories for each case
%
%      vbd_cat5 - Mx1 vector of BIRADS 5th edition categories for each case
%
%      vbd_mean - Mx1 vector of mean VBD scores for ecah case
%
%      vbd_lr_max - Mx1 vector max left/right VBD scores for each case
%
%
% Example:
%
% Notes:
%
% See also: VOLPARA_CATEGORISE_VBD_BATCH
%
% Created: 04-Sep-2015
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

%Check inputs are same length
if length(vbd_scores) ~= length(case_ids) ||...
        length(vbd_scores) ~= length(view_labels)
    error('vbd_scores, case_ids and view labels must all be the same size');
end

use_default_args = isempty(varargin);

%Get individual case ids
unique_case_ids = unique(case_ids);
num_cases = length(unique_case_ids);

%Pre-allocate outputs
vbd_cat4 = zeros(num_cases, 1);
vbd_cat5 = zeros(num_cases, 1);
vbd_mean = zeros(num_cases, 1);
vbd_lr_max = zeros(num_cases, 1);

%Loop through each case
for i_case = 1:num_cases
    case_idx = strcmp(case_ids, unique_case_ids(i_case));
    vbd_scores_i = vbd_scores(case_idx);
    view_labels_i = view_labels(case_idx);
    
    view_flags = true(4,1);
    vbd_scores_full = zeros(4,1);
    
    %Check which views are present for this case
    lcc_idx = find(strcmpi(view_labels_i, 'LCC')); %case-insensitive
    if isempty(lcc_idx)
        view_flags(1) = 0;
    elseif length(lcc_idx) > 1
        display(['Duplicate LCC views found for case ' unique_case_ids(i_case) ', skipping case']);
    else
        vbd_scores_full(1) = vbd_scores_i(lcc_idx);
    end
    lml_idx = find(strcmpi(view_labels_i, 'LML') | strcmpi(view_labels_i, 'LMLO')); %case-insensitive
    if isempty(lml_idx)
        view_flags(2) = 0;
    elseif length(lml_idx) > 1
        display(['Duplicate LML views found for case ' unique_case_ids(i_case) ', skipping case']);
    else
        vbd_scores_full(2) = vbd_scores_i(lml_idx);
    end
    rcc_idx = find(strcmpi(view_labels_i, 'RCC')); %case-insensitive
    if isempty(rcc_idx)
        view_flags(3) = 0;
    elseif length(rcc_idx) > 1
        display(['Duplicate RCC views found for case ' unique_case_ids(i_case) ', skipping case']);
    else
        vbd_scores_full(3) = vbd_scores_i(lcc_idx);
    end
    rml_idx = find(strcmpi(view_labels_i, 'RML') | strcmpi(view_labels_i, 'RMLO')); %case-insensitive
    if isempty(rml_idx)
        view_flags(4) = 0;
    elseif length(rml_idx) > 1
        display(['Duplicate RML views found for case ' unique_case_ids(i_case) ', skipping case']);
    else
        vbd_scores_full(4) = vbd_scores_i(rml_idx);
    end
    
    %Compute categories for this case and save in the output containers
    if use_default_args
        [vbd_cat4(i_case), vbd_cat5(i_case), vbd_mean(i_case), vbd_lr_max(i_case)] = ...
            volpara_categorise_vbd(vbd_scores_full, view_flags);
    else
        [vbd_cat4(i_case), vbd_cat5(i_case), vbd_mean(i_case), vbd_lr_max(i_case)] = ...
            volpara_categorise_vbd(vbd_scores_full, view_flags, varargin);
    end
end


