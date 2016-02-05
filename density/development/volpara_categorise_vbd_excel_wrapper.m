function [unique_case_ids, vdg4, vdg5, vbd_mean, vbd_lr_max,...
    vbd_scores_by_case, views_present, outliers] = ...
    volpara_categorise_vbd_excel_wrapper(xls_file_in, xls_file_out)
%VOLPARA_CATEGORISE_VBD_EXCEL_WRAPPER Apply VOLPARA_CATEGORISE_VBD_BATCH to input
% contained in Microsoft EXcel spreadsheet containing Volpara output
%   [unique_case_ids, vdg4, vdg5, vbd_mean, vbd_lr_max] = ...
%       volpara_categorise_vbd_excel_wrapper(xls_file_in, xls_file_out)
%
% Inputs:
%      xls_file_in - filepath to input spreadsheet
%
%      xls_file_out - filepath to spreadsheet in which case-by-case output
%      will be saved
%
% Outputs: These are the same as VOLPARA_CATEGORISE_VBD_BATCH
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
% Notes:
%
% See also: VOLPARA_CATEGORISE_VBD_BATCH VOLPARA_CATEGORISE_VBD
%
% Created: 04-Sep-2015
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

[~,~,by_image_data] = xlsread(xls_file_in, 'VolparaOutput');

%Get rid of any blank rows - although usually the headings are on row 1
while isnan(by_image_data{1})
    by_image_data(1,:) = [];
end
image_headings = by_image_data(1,:);
by_image_data(1,:) = [];

%Get columns we need
case_id_col = strcmpi(image_headings, 'PatientID');
detector_id_col = strcmpi(image_headings, 'DetectorID');
study_date_col = strcmpi(image_headings, 'StudyDate');
breast_label_col = strcmpi(image_headings, 'BreastSide');
view_label_col = strcmpi(image_headings, 'MammoView');
vbd_score_col = strcmpi(image_headings, 'VolumetricBreastDensity');

%There shouldn't be anything else below the rows of image data, so we'll
%assume all the remaining rows are valid

%For the case IDs, use a combination of patient ID, detector and study date
%to hopefully end up with unique cases
case_ids = [by_image_data(:, case_id_col)...
            by_image_data(:, detector_id_col)...
            cellstr(num2str(cell2mat(by_image_data(:, study_date_col))))];
breast_labels = by_image_data(:, breast_label_col);
view_labels = by_image_data(:, view_label_col);
vbd_scores = cell2mat(by_image_data(:, vbd_score_col));

%Merge breast side and view label to get LCC, LMLO, RCC, RMLO labels
l_breast_idx = strcmpi(breast_labels, 'left');
r_breast_idx = strcmpi(breast_labels, 'right');
view_labels(l_breast_idx) = strcat('L', view_labels(l_breast_idx));
view_labels(r_breast_idx) = strcat('R', view_labels(r_breast_idx));

[unique_case_ids, vdg4, vdg5, vbd_mean, vbd_lr_max,...
    vbd_scores_by_case, views_present, outliers] = ...
    volpara_categorise_vbd_batch(vbd_scores, case_ids, view_labels);

%If we've been a file to save the output to, write out in Excel format
if exist('xls_file_out', 'var') && ~isempty(xls_file_out)
    
    headings = {...
        'Patient ID',...
        'Detector ID',...
        'Study date',...
        'VBD% LCC',...
        'VBD% LML',...
        'VBD% RCC',...
        'VBD% RML',...
        'Mean VBD%',...
        'Max VBD%',...
        'VDG 4th Edition',...
        'VDG 5th Edition',...
        'Outlier LCC',...
        'Outlier LML',...
        'Outlier RCC',...
        'Outlier RML'};
    
    xls_data = [vbd_scores_by_case vbd_mean vbd_lr_max vdg4 vdg5 double(outliers)];
    
    xlswrite(xls_file_out, headings, 1, 'A1');
    xlswrite(xls_file_out, unique_case_ids, 1, 'A2');
    xlswrite(xls_file_out, xls_data, 1, 'D2');
end
    
        