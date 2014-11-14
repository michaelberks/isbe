data_type = {'abnormals', 'normals'};

for ii = 1:2

    mkdir(['C:\isbe\asymmetry_project\data\k_line_maps\2004_screening\' data_type{ii} '\']);
    mkdir(['C:\isbe\asymmetry_project\data\k_ori_maps\2004_screening\' data_type{ii} '\']);
    mkdir(['C:\isbe\asymmetry_project\data\k_scale_maps\2004_screening\' data_type{ii} '\']);
      
    mammo_dir = ['C:\isbe\asymmetry_project\data\mammograms\\2004_screening_processed\' data_type{ii} '\'];
    mammo_list = dir([mammo_dir '*.mat']);
    mammo_names = get_mammo_info(mammo_list);
    
    for jj = 1:length(mammo_list)
        display(['processing image ' num2str(jj) ' of ' num2str(length(mammo_list))]);
        mammo = imresize(u_load([mammo_dir mammo_list(ii).name]), 0.5, 'bilinear');
        [line_map, orientation_map, scale_map] = gaussian_2nd_derivative_line(mammo, [1 2 4 8]);
        save(['C:\isbe\asymmetry_project\data\k_line_maps\2004_screening\' data_type{ii} '\'...
            mammo_names{ii} '_data.mat'], 'line_map');
        save(['C:\isbe\asymmetry_project\data\k_ori_maps\2004_screening\' data_type{ii} '\'...
            mammo_names{ii} '_data.mat'], 'orientation_map');
        save(['C:\isbe\asymmetry_project\data\k_scale_maps\2004_screening\' data_type{ii} '\'...
            mammo_names{ii} '_data.mat'], 'scale_map');
        
        clear mammo line_map orientation_map scale_map;
    end
end

mammo = imresize(u_load('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\o04_024RCC.mat'), 0.5, 'bilinear');


px_per_mm = 50/9;

scales = [3.2 4.0 5.4 6.2 7.8] * px_per_mm;

blob_map = zeros(size(mammo));
scale_map = zeros(size(mammo));
%%
for ii = 1:length(scales)
    blob_scale = imfilter(double(mammo), fspecial('log', round(6*scales(ii))+1, scales(ii)), 'replicate');
    
    figure; imagesc(blob_scale); axis image; colorbar;
    
    swap_idx = abs(blob_scale) > abs(blob_map);
    blob_map(swap_idx) = blob_scale(swap_idx);
    
    scale_map(swap_idx) = scales(ii);
end
%%

