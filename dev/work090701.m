read_density_form_batch(3, 'K:\isbe\density\form\readers\reader03\', 50, 'reader03.xls', 1);
read_density_form_batch(5, 'K:\isbe\density\form\readers\reader05\', 50, 'reader05.xls', 0);
read_density_form_batch(6, 'K:\isbe\density\form\readers\reader06\', 50, 'reader06.xls', 0);
read_density_form_batch(11, 'K:\isbe\density\form\readers\reader11\', 50, 'reader11.xls', 0);
%%
rox_list =dir('J:\stepwedge\mike''s stuff\data for rox (rox sample)\*data*');
%phr_list =dir('J:\stepwedge\mike''s stuff\data for phreddie (rox sample)\*data*');

for ii = 1:length(rox_list)
    rox_data = u_load(['J:\stepwedge\mike''s stuff\data for rox (rox sample)\', rox_list(ii).name]);
    phr_data = u_load(['J:\stepwedge\mike''s stuff\data for phreddie (rox sample)\', rox_list(ii).name]);
    
    density_data = rox_data;
    density_data.coarse_edgex = phr_data.coarse_edgex;
    density_data.coarse_edgey = phr_data.coarse_edgey;
    
    save(['J:\stepwedge\mike''s stuff\data_rp\', rox_list(ii).name], 'density_data');
    
    density_data = phr_data;
    density_data.coarse_edgex = rox_data.coarse_edgex;
    density_data.coarse_edgey = rox_data.coarse_edgey;
    
    save(['J:\stepwedge\mike''s stuff\data_pr\', rox_list(ii).name], 'density_data');
    
end
%%
rox_list =dir('J:\stepwedge\mike''s stuff\data for rox (rox sample)\*data*');
%phr_list =dir('J:\stepwedge\mike''s stuff\data for phreddie (rox sample)\*data*');

for ii = 1:length(rox_list)
    rox_data = u_load(['J:\stepwedge\mike''s stuff\data for rox (rox sample)\', rox_list(ii).name]);
    phr_data = u_load(['J:\stepwedge\mike''s stuff\data for phreddie (rox sample)\', rox_list(ii).name]);
    
    density_data = rox_data;
    density_data.wedgevals = phr_data.wedgevals;
    density_data.swx = phr_data.swx;
    density_data.swy = phr_data.swy;
    
    save(['J:\stepwedge\mike''s stuff\data_step_rp\', rox_list(ii).name], 'density_data');
    
    density_data = phr_data;
    density_data.wedgevals = rox_data.wedgevals;
    density_data.swx = rox_data.swx;
    density_data.swy = rox_data.swy;
    
    save(['J:\stepwedge\mike''s stuff\data_step_pr\', rox_list(ii).name], 'density_data');
    
end

%stepwedge_from_saved('SelectionMode', 2);
%stepwedge_from_saved('SelectionMode', 2);
%
density_results_to_excel('J:\stepwedge\mike''s stuff\results_step_pr\', 'results_step_pr.xls');
density_results_to_excel('J:\stepwedge\mike''s stuff\results_step_rp\', 'results_step_rp.xls');
%%
mkdir('J:\stepwedge\mike''s stuff\data_marker_rp\');
mkdir('J:\stepwedge\mike''s stuff\results_marker_rp\'); 
mkdir('J:\stepwedge\mike''s stuff\data_marker_pr\'); 
mkdir('J:\stepwedge\mike''s stuff\results_marker_pr\'); 
for ii = 1:length(rox_list)
    rox_data = u_load(['J:\stepwedge\mike''s stuff\data for rox (rox sample)\', rox_list(ii).name]);
    phr_data = u_load(['J:\stepwedge\mike''s stuff\data for phreddie (rox sample)\', rox_list(ii).name]);
    
    density_data = rox_data;
    density_data.x_b_info = phr_data.x_b_info;
    
    save(['J:\stepwedge\mike''s stuff\data_marker_rp\', rox_list(ii).name], 'density_data');
    
    density_data = phr_data;
    density_data.x_b_info = rox_data.x_b_info;
    
    save(['J:\stepwedge\mike''s stuff\data_marker_pr\', rox_list(ii).name], 'density_data');
    
end

stepwedge_from_saved('SelectionMode', 2);
stepwedge_from_saved('SelectionMode', 2);
density_results_to_excel('J:\stepwedge\mike''s stuff\results_marker_pr\', 'results_marker_pr.xls');
density_results_to_excel('J:\stepwedge\mike''s stuff\results_marker_rp\', 'results_marker_rp.xls');