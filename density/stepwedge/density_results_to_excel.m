function [result_names density_results_all] = density_results_to_excel(results_dir, excel_file)

if ~strcmp(results_dir(end), filesep)
    results_dir = [results_dir, filesep];
end

results_list = dir([results_dir, '*.mat']);
num_results = length(results_list);

density_results_all = zeros(num_results, 18);
result_names = cell(num_results, 1);

for ii = 1:num_results
    
    load([results_dir, results_list(ii).name]);
    
    result_names{ii} = correct_name(results_list(ii).name(1:end-12), 4);
    if isfield(density_results, 'volume_b'); density_results_all(ii, 1) = density_results.volume_b; end
    if isfield(density_results, 'volume_g'); density_results_all(ii, 2) = density_results.volume_g; end
    if isfield(density_results, 'dense_by_vol_1'); density_results_all(ii, 3) = density_results.dense_by_vol_1; end
    if isfield(density_results, 'area_b'); density_results_all(ii, 4) = density_results.area_b; end
    if isfield(density_results, 'area_g'); density_results_all(ii, 5) = density_results.area_g; end
    if isfield(density_results, 'area_g5'); density_results_all(ii, 6) = density_results.area_g5; end
    if isfield(density_results, 'area_g10'); density_results_all(ii, 7) = density_results.area_g10; end
    if isfield(density_results, 'area_g15'); density_results_all(ii, 8) = density_results.area_g15; end
    if isfield(density_results, 'area_g20'); density_results_all(ii, 9) = density_results.area_g20; end
    if isfield(density_results, 'area_g25'); density_results_all(ii, 10) = density_results.area_g25; end
    if isfield(density_results, 'dense_by_area_1'); density_results_all(ii, 11) = density_results.dense_by_area_1; end
    if isfield(density_results, 'density_area_g5'); density_results_all(ii, 12) = density_results.density_area_g5; end
    if isfield(density_results, 'density_area_g10'); density_results_all(ii, 13) = density_results.density_area_g10; end
    if isfield(density_results, 'density_area_g15'); density_results_all(ii, 14) = density_results.density_area_g15; end
    if isfield(density_results, 'density_area_g20'); density_results_all(ii, 15) = density_results.density_area_g20; end
    if isfield(density_results, 'density_area_g25'); density_results_all(ii, 16) = density_results.density_area_g25; end
    if isfield(density_results, 'min_thickness'); density_results_all(ii, 17) = density_results.min_thickness; end
    if isfield(density_results, 'max_thickness'); density_results_all(ii, 18) = density_results.max_thickness; end
end

xlswrite(excel_file, {'Name'}, 1, 'A1');
xlswrite(excel_file, {'Volume breast'}, 1, 'B1');
xlswrite(excel_file, {'Volume gland'}, 1, 'C1');
xlswrite(excel_file, {'Volume gland (%)'}, 1, 'D1');
xlswrite(excel_file, {'Area breast'}, 1, 'E1');
xlswrite(excel_file, {'Area gland > 0'}, 1, 'F1');
xlswrite(excel_file, {'Area gland > 5'}, 1, 'G1');
xlswrite(excel_file, {'Area gland > 10'}, 1, 'H1');
xlswrite(excel_file, {'Area gland > 15'}, 1, 'I1');
xlswrite(excel_file, {'Area gland > 20'}, 1, 'J1');
xlswrite(excel_file, {'Area gland > 25'}, 1, 'K1');
xlswrite(excel_file, {'Area gland > 0 (%)'}, 1, 'L1');
xlswrite(excel_file, {'Area gland > 5 (%)'}, 1, 'M1');
xlswrite(excel_file, {'Area gland > 10 (%)'}, 1, 'N1');
xlswrite(excel_file, {'Area gland > 15 (%)'}, 1, 'O1');
xlswrite(excel_file, {'Area gland > 20 (%)'}, 1, 'P1');
xlswrite(excel_file, {'Area gland > 25 (%)'}, 1, 'Q1');
xlswrite(excel_file, {'Min thickness'}, 1, 'R1');
xlswrite(excel_file, {'Max thickness'}, 1, 'S1');

xlswrite(excel_file, result_names, 1, 'A2');
xlswrite(excel_file, density_results_all, 1, 'B2');

function corrected_name = correct_name(name, len)
    
    views = {'LCC', 'RCC', 'LML', 'RML'}; 
    pos = []; ii = 1; 
    while isempty(pos); pos = strfind(name,views{ii}); ii = ii+1; end
    corrected_name = name;
    for ii = 1:len - pos + 1
        corrected_name = ['0', corrected_name]; %#ok
    end
    