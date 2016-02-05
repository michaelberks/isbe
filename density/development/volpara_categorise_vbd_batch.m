function [unique_case_ids, vdg4, vdg5, vbd_mean, vbd_lr_max,...
    vbd_scores_by_case, views_present, outliers] = ...
    volpara_categorise_vbd_batch(vbd_scores, case_ids, view_labels, varargin)
%VOLPARA_CATEGORISE_VBD_BATCH Apply Volpara categorisation into VDG groups
% to a set of input images 
%   [unique_case_ids, vdg4, vdg5, vbd_mean, vbd_lr_max] = ...
%       volpara_categorise_vbd_batch(vbd_scores, case_ids, view_labels, varargin)
%
% Inputs:
%      vbd_scores - Nx1 array of VBD scores
%
%      case_ids - Nxk cell array of case IDs. Each row should be a unique case
%      k can be as many data fields as necessary to identify a unique case
%      In the excel wrapper I use 3: patientID, detectorID and study date.
%      Note the case IDs do not need to be contiguous, as long each case
%      ID matches up with the image VBD score and view label they will be
%      grouped together
%
%      view_labels - Nx1 cell array of view labels (each entry must be one
%      of {LCC, LML, RCC, RML} or since I'm being nice, lower-case variants
%      or or L(R)MLO.
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
%      vdg4 - Mx1 vector of BIRADS 4th edition categories for each case
%
%      vdg5 - Mx1 vector of BIRADS 5th edition categories for each case
%
%      vbd_mean - Mx1 vector of mean VBD scores for ecah case
%
%      vbd_lr_max - Mx1 vector max left/right VBD scores for each case
%
%       vbd_scores_by_case - Mx4 array of the input VBD scores for each
%       case, in order LCC, LML, RCC, RML
%
%       views_present - Mx4 logical array flagging which views present
%
%       outliers - Mx4 logical array flagging which views marked as outliers
%
%
% Example:
%
% Notes: This function directly mimics the Volpara macro. One strange
% feature of this is when a case has multiple images of the same view, the
% Volpara macro overwrites each earlier image it detects, effectively only
% using the last (in terms of row number in the excel input) image for each
% view. This doesn't seem terribly satisfactory to me - so I warn when this
% happens, but for consistency I follow this protocol. The function can be
% called directly on an excel spreadsheet of image-by-image data produced
% by Volpara using the excel wrapper below
%
% See also: VOLPARA_CATEGORISE_VBD VOLPARA_CATEGORISE_VBD_EXCEL_WRAPPER
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

%Discard any missing values (will be NaNs in the vbd_scores) now - note if
%we do this before taking the unique case, we'll discard any cases with no
%views from the analysis, this matches the volpara macro output, but may be
%something we want to flag?
missing_idx = isnan(vbd_scores);
vbd_scores(missing_idx) = [];
case_ids(missing_idx,:) = [];
view_labels(missing_idx) = [];
[num_images num_features] = size(case_ids);

%Get individual case ids
unique_case_ids = unique_cell_rows(case_ids);
num_cases = size(unique_case_ids,1);

%Pre-allocate outputs
vdg4 = zeros(num_cases, 1);
vdg5 = zeros(num_cases, 1);
vbd_mean = zeros(num_cases, 1);
vbd_lr_max = zeros(num_cases, 1);
vbd_scores_by_case = zeros(num_cases, 4);
views_present = false(num_cases, 4);
outliers = false(num_cases, 4);

%Loop through each case
for i_case = 1:num_cases
    
    %Get indices that match this unique case
    case_idx = true(num_images,1);    
    for i_feature = 1:num_features
        case_idx = case_idx & strcmp(case_ids(:,i_feature), unique_case_ids(i_case,i_feature));
    end
    
    vbd_scores_i = vbd_scores(case_idx);
    view_labels_i = view_labels(case_idx);
    
    view_flags = true(4,1);
    vbd_scores_full = zeros(4,1);
    
    %Check which views are present for this case
    lcc_idx = find(strcmpi(view_labels_i, 'LCC')); %case-insensitive
    if isempty(lcc_idx)
        view_flags(1) = 0;
    else
        if length(lcc_idx) > 1
            display(['Duplicate LCC views found for case ' unique_case_ids{i_case}]);
        end
        vbd_scores_full(1) = vbd_scores_i(lcc_idx(end));
    end
    lml_idx = find(strcmpi(view_labels_i, 'LML') | strcmpi(view_labels_i, 'LMLO')); %case-insensitive
    if isempty(lml_idx)
        view_flags(2) = 0;
    else
        if length(lml_idx) > 1
            display(['Duplicate LML views found for case ' unique_case_ids{i_case}]);
        end
        vbd_scores_full(2) = vbd_scores_i(lml_idx(end));
    end
    rcc_idx = find(strcmpi(view_labels_i, 'RCC')); %case-insensitive
    if isempty(rcc_idx)
        view_flags(3) = 0;
    else
        if length(rcc_idx) > 1
            display(['Duplicate RCC views found for case ' unique_case_ids{i_case}]);
        end
        vbd_scores_full(3) = vbd_scores_i(rcc_idx(end));
    end
    rml_idx = find(strcmpi(view_labels_i, 'RML') | strcmpi(view_labels_i, 'RMLO')); %case-insensitive
    if isempty(rml_idx)
        view_flags(4) = 0;
    else
        if length(rml_idx) > 1
            display(['Duplicate RML views found for case ' unique_case_ids{i_case}]);
        end
        vbd_scores_full(4) = vbd_scores_i(rml_idx(end));
    end
    
    %Compute categories for this case and save in the output containers
    if use_default_args
        [vdg4(i_case), vdg5(i_case), vbd_mean(i_case), vbd_lr_max(i_case), outliers(i_case,:)] = ...
            volpara_categorise_vbd(vbd_scores_full, view_flags);
    else
        [vdg4(i_case), vdg5(i_case), vbd_mean(i_case), vbd_lr_max(i_case), outliers(i_case,:)] = ...
            volpara_categorise_vbd(vbd_scores_full, view_flags, varargin);
    end
    views_present(i_case,:) = view_flags;
    vbd_scores_by_case(i_case,:) = vbd_scores_full;
end

%--------------------------------------------------------------------------
%For some annoying reason Matlab doesn't support unique with 'rows' option
%for cell arrays
function unique_cell = unique_cell_rows(cell_in)

[num_rows num_cols] = size(cell_in);
if num_cols==1
    unique_cell = unique(cell_in);
    return;
end
unique_idx = zeros(num_rows, num_cols);
for i_col = 1:num_cols
    [~,~,unique_idx(:,i_col)] = unique(cell_in(:,i_col));
end
[~,ui] = unique(unique_idx, 'rows');
unique_cell = cell_in(ui,:);




