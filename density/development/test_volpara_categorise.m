%--------------------------------------------------------------------------
% Script name: test_volpara_categorise
%--------------------------------------------------------------------------
%
% Script to test the functions VOLPARA_CATEGORISE_VBD and
% VOLPARA_CATEGORISE_VBD_BATCH produce the same output as the Excel macro
% provided by Volpara.
%
% Takes as input the a spreadsheet containing image-by-image Volpara data.
% Assumes the image-by-image data is in a worksheet named VolparaOuput and
% that the Volpara macro has been applied to this worksheet generating
% a worksheet (amongst others) with case-by-case data called ByCase
%
% Procudes comparisons of the case VBD% scores computed from the individual
% views and the resulting categorisations into density groups, for both the
% BIRADS 4th and 5th editions thresholds (as defined by Volpara)
%
%--------------------------------------------------------------------------
%xls_name = 'C:\isbe\density\assure\CANCERS_PROCAS_VOLPARA_COMBO.xls';
xls_name = 'C:\isbe\density\assure\CONTROLS_PROCAS_VOLPARA_COMBO.xlsx';
%%
%Get image-by-image input from volpara spreadsheet and run volpara
%categorise
[unique_case_ids, vdg4, vdg5, vbd_mean, vbd_lr_max,...
    vbd_scores_by_case, views_present, outliers] = ...
    volpara_categorise_vbd_excel_wrapper(xls_name);
%%
%Get macro output from volpara spreadsheet
[~,~,by_case_data] = xlsread(xls_name, 'ByCase');

%Get rid of any blank rows - usually the headings are on row 2
while isnan(by_case_data{1})
    by_case_data(1,:) = [];
end
case_headings = by_case_data(1,:);
by_case_data(1,:) = [];

end_row = 1;
while ~isnan(by_case_data{end_row,1})
    end_row = end_row + 1;
end
num_cases = end_row-1;

%Get columns we need
case_id_col = find(strcmpi(case_headings, 'PatientID'));
vdg4_col = find(strcmpi(case_headings, 'VDG 4thEd'));
vdg5_col = find(strcmpi(case_headings, 'VDG 5thEd'));
vbd_max_col = find(strcmpi(case_headings, 'Max VBD%'));
vbd_mean_col = find(strcmpi(case_headings, 'VBD%'));

case_ids_xls = by_case_data(1:num_cases, case_id_col);
vdg4_xls = cell2mat(by_case_data(1:num_cases, vdg4_col));
vdg5_xls = cell2mat(by_case_data(1:num_cases, vdg5_col));
vbd_max_xls = cell2mat(by_case_data(1:num_cases, vbd_max_col));
vbd_mean_xls = cell2mat(by_case_data(1:num_cases, vbd_mean_col));
    
%We know the cases returned from my matlab function will be sorted in case
%order, check the ones we've pulled from excel are
if ~issorted(case_ids_xls)
    [case_ids_xls, sorted_idx] = sort(case_ids_xls);
    vdg4_xls = vdg4_xls(sorted_idx);
    vdg5_xls = vdg5_xls(sorted_idx);
    vbd_max_xls = vbd_max_xls(sorted_idx);
    vbd_mean_xls = vbd_mean_xls(sorted_idx);
end
%%
%Now check we're getting the same result

%First of all, are the case IDs the same, if not, we're stuffed already
if size(unique_case_ids,1) == num_cases
    matching_case_ids = strcmp(case_ids_xls, unique_case_ids(:,1));
    if all(matching_case_ids)
        display('All case IDs match');
        %
        mismatched_mean_vbd = abs(vbd_mean_xls - vbd_mean) > 1e-4;
        mismatched_max_vbd = abs(vbd_max_xls - vbd_lr_max) > 1e-4;

        mismatched_vdg4 = vdg4_xls ~= vdg4;
        mismatched_vdg5 = vdg5_xls ~= vdg5;

        %Use sparse trick to fill co-occurrence matrices
        co_occurrence_vdg4 = full(sparse(vdg4_xls, vdg4, 1, 4, 4));
        co_occurrence_vdg5 = full(sparse(vdg5_xls, vdg5, 1, 4, 4));

        display(['# Cases with mis-matched mean VBD%: ' num2str(sum(mismatched_mean_vbd)) ' of ' num2str(num_cases)]);
        display(['# Cases with mis-matched max VBD%: ' num2str(sum(mismatched_max_vbd)) ' of ' num2str(num_cases)]);

        display(['# Cases with mis-matched BIRADS 4th edition VDG: ' num2str(sum(mismatched_vdg4)) ' of ' num2str(num_cases)]);
        display(['# Cases with mis-matched BIRADS 5th edition VDG: ' num2str(sum(mismatched_vdg5)) ' of ' num2str(num_cases)]);

        display('Co-occurrence matrix for  BIRADS 4th edition VDG:')
        display(co_occurrence_vdg4);

        display('Co-occurrence matrix for  BIRADS 5th edition VDG:')
        display(co_occurrence_vdg5);
    else
        display('Case IDs do not match');
        display('In Excel, but not in Matlab:')
        display(setdiff(case_ids_xls, unique_case_ids(:,1)));
        display('In Matlab, but not in Excel:')
        display(setdiff(unique_case_ids(:,1), case_ids_xls));
    end
else
    display('Number of cases do not match');
    display('In Excel, but not in Matlab:')
    display(setdiff(case_ids_xls, unique_case_ids(:,1)));
    display('In Matlab, but not in Excel:')
    display(setdiff(unique_case_ids(:,1), case_ids_xls));
end

