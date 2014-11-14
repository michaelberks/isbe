density_results_to_excel('J:\stepwedge\results_rosanne\', 'results_rosanne.xls');
density_results_to_excel('J:\stepwedge\results_phreddie\', 'results_phreddie.xls');
density_results_to_excel('J:\stepwedge\results_camilla\', 'results_camilla.xls');

density_results_to_excel('J:\stepwedge\results for phreddie (rox sample)\', 'results_for_phreddie_(rox_sample).xls');
density_results_to_excel('J:\stepwedge\results for rox (rox sample)\', 'results_for_rox_(rox_sample).xls');

for ii = 1:9
    density_results_to_excel(...
        ['J:\stepwedge\results for rox (rox sample)\repeats\repeat ', num2str(ii), '\'],...
        ['repeat ', num2str(ii), '.xls']);
end
%%
%stepwedge_failures('DataPath', 'J:\stepwedge\data', 'ExcelFile', 'failures.xls');
stepwedge_failures('DataPath', 'J:\stepwedge\individual_data\data_phreddie', 'ExcelFile', 'failures.xls');
stepwedge_failures('DataPath', 'J:\stepwedge\individual_data\data_rosanne', 'ExcelFile', 'failures.xls');
stepwedge_failures('DataPath', 'J:\stepwedge\data for phreddie (rox sample)', 'ExcelFile', 'failures.xls');
stepwedge_failures('DataPath', 'J:\stepwedge\data for rox (rox sample)', 'ExcelFile', 'failures.xls');
%%
