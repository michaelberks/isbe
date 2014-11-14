mkdir C:\isbe\asymmetry_project\data\synthetic_lines\lines512\g2d_responses

for ii = 1:100
    display(['processing image', num2str(ii)]);
    s = load([asymmetryroot 'data/synthetic_lines/lines512/image' zerostr(ii,3) '.mat']);
    test_image = s.test_image; clear s;
    
    g2d_responses = compute_gaussian_2nd_derivatives(test_image,  [1 2 4 8]);
    
    save([asymmetryroot 'data/synthetic_lines/lines512/g2d_responses/g2d_response' zerostr(ii,3) '.mat'],...
        'g2d_responses');
    
    clear test_image g2d_responses
end
%%
for dist = [16 32 64 128 256] 
%     [P_x x] = mass_maps_hist('2004_screening/abnormals', 'radial_maps_wrf', 'view', zerostr(dist,3), 'sigma', 8); 
%     save(['C:\isbe\asymmetry_project\experiments\radial_maps\cdf_wrf_' zerostr(dist,3) '.mat'], 'P_x', 'x');
%     
%     [P_x x] = mass_maps_hist('2004_screening/abnormals', 'radial_maps_rf', 'view', zerostr(dist,3), 'sigma', 8);
%     save(['C:\isbe\asymmetry_project\experiments\radial_maps\cdf_rf_' zerostr(dist,3) '.mat'], 'P_x', 'x');
    
    [P_x x] = mass_maps_hist('2004_screening_processed/abnormals', 'radial_maps_g2', 'view', zerostr(dist,3), 'sigma', 4);
    save(['C:\isbe\asymmetry_project\experiments\radial_maps\cdf_g2_' zerostr(dist,3) '.mat'], 'P_x', 'x');
    
    [P_x x] = mass_maps_hist('2004_screening_processed/abnormals', 'radial_maps_wg2', 'view', zerostr(dist,3), 'sigma', 4);
    save(['C:\isbe\asymmetry_project\experiments\radial_maps\cdf_wg2_' zerostr(dist,3) '.mat'], 'P_x', 'x');
end
%%
map_types = {'f1_', 'f2_'};
for ii =1:2
    for dist = [60 90 120 150 180]

        [P_x x] = mass_maps_hist('2004_screening/abnormals', 'k_stellate_maps_wrf', 'view', [map_types{ii} zerostr(dist,3)], 'sigma', 8); 
        save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\cdf_wrf_' [map_types{ii} zerostr(dist,3)] '.mat'], 'P_x', 'x');

        [P_x x] = mass_maps_hist('2004_screening/abnormals', 'k_stellate_maps_rf', 'view', [map_types{ii} zerostr(dist,3)], 'sigma', 8);
        save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\cdf_rf_' [map_types{ii} zerostr(dist,3)] '.mat'], 'P_x', 'x');

        [P_x x] = mass_maps_hist('2004_screening_processed/abnormals', 'k_stellate_maps_g2', 'view', [map_types{ii} zerostr(dist,3)], 'sigma', 8);
        save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\cdf_g2_' [map_types{ii} zerostr(dist,3)] '.mat'], 'P_x', 'x');

    end
end
%%
dists = [16 32 64 128 256];
f_x = cell(4,5);
f_P_x = cell(4,5);

for jj = 1:5
    
    dist = dists(jj);
    
    load(['C:\isbe\asymmetry_project\experiments\radial_maps\cdf_rf_' zerostr(dist,3) '.mat']);
    f_x{1, jj} = x;
    f_P_x{1, jj} = cumsum(P_x) / sum(P_x);

    load(['C:\isbe\asymmetry_project\experiments\radial_maps\cdf_wrf_' zerostr(dist,3) '.mat']);
    f_x{2, jj} = x;
    f_P_x{2, jj} = cumsum(P_x) / sum(P_x);

    load(['C:\isbe\asymmetry_project\experiments\radial_maps\cdf_g2_' zerostr(dist,3) '.mat']);
    f_x{3, jj} = x;
    f_P_x{3, jj} = cumsum(P_x) / sum(P_x);

    load(['C:\isbe\asymmetry_project\experiments\radial_maps\cdf_wg2_' zerostr(dist,3) '.mat']);
    f_x{4, jj} = x;
    f_P_x{4, jj} = cumsum(P_x) / sum(P_x);
    clear x P_x
    
end
%%
image_types = {'abnormals', 'normals'};
map_types = {'radial_maps_rf', 'radial_maps_wrf', 'radial_maps_g2', 'radial_maps_wg2'};
mam_types = {'\2004_screening', '\2004_screening', '\2004_screening_processed', '\2004_screening_processed'};
dists = [16 32 64 128 256];
for ll = 1:2
    %
    mam_names = u_load(['C:\isbe\asymmetry_project\data\mam_names\2004_screening_' image_types{ll}]);
    for kk = 3:4
        
        mkdir(['C:\isbe\asymmetry_project\data\' map_types{kk} '_cdf' mam_types{kk} '\' image_types{ll} '\']);
        
        for ii = 1:length(mam_names)

            for jj = 1:5

                dd = dists(jj);
                
                f_j = load_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} mam_types{kk}...
                    '\' image_types{ll} '\' mam_names{ii} '_rad_map_' zerostr(dd,3) '.mat']);

                f_j = imfilter(f_j, fspecial('gaussian', 25, 4), 'replicate');

                f_map = interp1(f_x{kk,jj}, f_P_x{kk,jj}, f_j, 'linear');
                f_map(f_j < min(f_x{kk,jj})) = 0;
                f_map(f_j > max(f_x{kk,jj})) = 1;
            
                save_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} '_cdf' mam_types{kk}...
                    '\' image_types{ll} '\' mam_names{ii} '_rad_map_' zerostr(dd,3) '.mat'], f_map);

            end
        end
    end
end
%%
image_types = {'abnormals', 'normals'};
map_types = {'radial_maps_rf', 'radial_maps_wrf', 'radial_maps_g2', 'radial_maps_wg2'};
mam_types = {'\2004_screening', '\2004_screening', '\2004_screening_processed'};
dists = [16 32 64 128 256];
temp_dists = [30 45 60 75 90];
for ll = 1:2
%     mkdir(['C:\isbe\asymmetry_project\data\radial_maps_rf_cdf\2004_screening\' image_types{ll} '\']);
%     mkdir(['C:\isbe\asymmetry_project\data\radial_maps_wrf_ts\2004_screening\' image_types{ll} '\']);
%     mkdir(['C:\isbe\asymmetry_project\data\radial_maps_g2_ts\2004_screening_processed\' image_types{ll} '\']);
%     mkdir(['C:\isbe\asymmetry_project\data\radial_maps_g2_ts\2004_screening_processed\' image_types{ll} '\']);
    %
    mam_names = u_load(['C:\isbe\asymmetry_project\data\mam_names\2004_screening_' image_types{ll}]);
    for kk = 1:2
        
        mkdir(['C:\isbe\asymmetry_project\data\' map_types{kk} '_cdf' mam_types{kk} '\' image_types{ll} '\']);
        mkdir(['C:\isbe\asymmetry_project\data\' map_types{kk} '_cdf64' mam_types{kk} '\' image_types{ll} '\']);
        
        for ii = 1:length(mam_names)

            f_map = load_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} mam_types{kk}...
                '\' image_types{ll} '\' mam_names{ii} '_rad_map_064.mat']);

            f_map = imfilter(f_map, fspecial('gaussian', 49, 8), 'replicate');

            f_map_cdf = interp1(f_x{kk,3}, f_P_x{kk,3}, f_map, 'linear');
            f_map_cdf(f_map < min(f_x{kk,3})) = 0;
            f_map_cdf(f_map > max(f_x{kk,3})) = 1;
            
            save_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} '_cdf64' mam_types{kk}...
                '\' image_types{ll} '\' mam_names{ii} '_rad_map.mat'], f_map_cdf);
            
            scale_map = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening_processed\' image_types{ll} ...
                '\scales\' mam_names{ii} '_scales.mat']);
            template_map = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening_processed\' image_types{ll}...
                '\' mam_names{ii} '_template.mat']);
            template_map = template_map > 0.3;

            for jj = [1 2 4 5]
                dd = dists(jj);
                swap_idx = (scale_map == temp_dists(jj)) & template_map;
                
                f_j = load_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} mam_types{kk}...
                    '\' image_types{ll} '\' mam_names{ii} '_rad_map_' zerostr(dd,3) '.mat']);

                f_j = imfilter(f_j, fspecial('gaussian', 49, 8), 'replicate');

                f_j_cdf = interp1(f_x{kk,jj}, f_P_x{kk,jj}, f_j, 'linear');
                f_j_cdf(f_j < min(f_x{kk,jj})) = 0;
                f_j_cdf(f_j > max(f_x{kk,jj})) = 1;
            
                f_map_cdf(swap_idx) = f_j_cdf(swap_idx);

            end
%             figure('name', map_types{kk}); 
%             imagesc(f_map_cdf); axis image; colorbar;

            save_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} '_cdf' mam_types{kk}...
                '\' image_types{ll} '\' mam_names{ii} '_rad_map.mat'], f_map_cdf);

        end
    end
end
%%
for ll = 1:2
    mkdir(['C:\isbe\asymmetry_project\data\k_stellate_maps_rf_ms\2004_screening\' image_types{ll} '\']);
    mkdir(['C:\isbe\asymmetry_project\data\k_stellate_maps_wrf_ms\2004_screening\' image_types{ll} '\']);
    mkdir(['C:\isbe\asymmetry_project\data\k_stellate_maps_g2_ms\2004_screening_processed\' image_types{ll} '\']);
    %
    mam_names = u_load(['C:\isbe\asymmetry_project\data\mam_names\2004_screening_' image_types{ll}]);
    for kk = 1:3
        
        for ii = 1:2%length(mam_names)

            f1_map = load_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} mam_types{kk}...
                '\' image_types{ll} '\' mam_names{ii} '_f1_120.mat']);
            f2_map = load_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} mam_types{kk}...
                '\' image_types{ll} '\' mam_names{ii} '_f2_120.mat']);

            f1_map = imfilter(f1_map, fspecial('gaussian', 49, 8), 'replicate');
            f2_map = imfilter(f2_map, fspecial('gaussian', 49, 8), 'replicate');

            f1_map = interp1(f1_x{kk,3}, f1_P_x{kk,3}, f1_map, 'linear');
            f2_map = interp1(f2_x{kk,3}, f2_P_x{kk,3}, f2_map, 'linear');

            for jj = [1 2 4 5]
                dd = dists(jj);
                
                f1_j = load_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} mam_types{kk}...
                    '\' image_types{ll} '\' mam_names{ii} '_f1_' zerostr(2*dd,3) '.mat']);
                f2_j = load_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} mam_types{kk}...
                    '\' image_types{ll} '\' mam_names{ii} '_f2_' zerostr(2*dd,3) '.mat']);

                f1_j = imfilter(f1_j, fspecial('gaussian', 49, 8), 'replicate');
                f2_j = imfilter(f2_j, fspecial('gaussian', 49, 8), 'replicate');

                f1_j = interp1(f1_x{kk,jj}, f1_P_x{kk,jj}, f1_j, 'linear');
                f2_j = interp1(f2_x{kk,jj}, f2_P_x{kk,jj}, f2_j, 'linear');

                swap_idx = f1_j > f1_map;
                f1_map(swap_idx) = f1_j(swap_idx);
                swap_idx = f2_j > f2_map;
                f2_map(swap_idx) = f2_j(swap_idx);

            end
            figure; 
            subplot(1,2,1); imagesc(f1_map); axis image; colorbar;
            subplot(1,2,1); imagesc(f2_map); axis image

%             save_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} '_scaled' mam_types{kk}...
%                 '\' image_types{ll} '\' mam_names{ii} '_f1.mat'], f1_map);
%             save_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} '_scaled' mam_types{kk}...
%                 '\' image_types{ll} '\' mam_names{ii} '_f2.mat'], f2_map);
        end
    end
end
%%
map_types = {'f1', 'f2'};
for ii =1:2

    [P_x x] = mass_maps_hist('2004_screening/normals', 'k_stellate_maps_wrf_scaled', 'view', map_types{ii}, 'sigma', 8); 
    save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\cdf_wrf_' map_types{ii} '_scaled.mat'], 'P_x', 'x');

    [P_x x] = mass_maps_hist('2004_screening/normals', 'k_stellate_maps_rf_scaled', 'view', map_types{ii}, 'sigma', 8);
    save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\cdf_rf_' map_types{ii} '_scaled.mat'], 'P_x', 'x');

    [P_x x] = mass_maps_hist('2004_screening_processed/normals', 'k_stellate_maps_g2_scaled', 'view', map_types{ii}, 'sigma', 8);
    save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\cdf_g2_' map_types{ii} '_scaled.mat'], 'P_x', 'x');

end
%% 
load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\cdf_rf_f1_scaled.mat');
f1_x{1} = x;
f1_P_x{1} = cumsum(P_x) / sum(P_x);

load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\cdf_wrf_f1_scaled.mat');
f1_x{2} = x;
f1_P_x{2} = cumsum(P_x) / sum(P_x);

load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\cdf_g2_f1_scaled.mat');
f1_x{3} = x;
f1_P_x{3} = cumsum(P_x) / sum(P_x);

load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\cdf_rf_f2_scaled.mat');
f2_x{1} = x;
f2_P_x{1} = cumsum(P_x) / sum(P_x);

load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\cdf_wrf_f2_scaled.mat');
f2_x{2} = x;
f2_P_x{2} = cumsum(P_x) / sum(P_x);

load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\cdf_g2_f2_scaled.mat');
f2_x{3} = x;
f2_P_x{3} = cumsum(P_x) / sum(P_x);

%
image_types = {'abnormals', 'normals'};
map_types = {'k_stellate_maps_rf_scaled', 'k_stellate_maps_wrf_scaled', 'k_stellate_maps_g2_scaled'};
mam_types = {'\2004_screening', '\2004_screening', '\2004_screening_processed'};
for ll = 1:2
    mkdir(['C:\isbe\asymmetry_project\data\k_stellate_maps_rf_cdf\2004_screening\' image_types{ll} '\']);
    mkdir(['C:\isbe\asymmetry_project\data\k_stellate_maps_wrf_cdf\2004_screening\' image_types{ll} '\']);
    mkdir(['C:\isbe\asymmetry_project\data\k_stellate_maps_g2_cdf\2004_screening_processed\' image_types{ll} '\']);
    %
    mam_names = u_load(['C:\isbe\asymmetry_project\data\mam_names\2004_screening_' image_types{ll}]);
    for kk = 1:3
        
        for ii = 1:length(mam_names)

            f1_map = load_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} mam_types{kk}...
                '\' image_types{ll} '\' mam_names{ii} '_f1.mat']);
            f2_map = load_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} mam_types{kk}...
                '\' image_types{ll} '\' mam_names{ii} '_f2.mat']);

            f1_map = imfilter(f1_map, fspecial('gaussian', 49, 8), 'replicate');
            f2_map = imfilter(f2_map, fspecial('gaussian', 49, 8), 'replicate');

            f1_map_cdf = interp1(f1_x{kk}, f1_P_x{kk}, f1_map, 'linear');
            f2_map_cdf = interp1(f2_x{kk}, f2_P_x{kk}, f2_map, 'linear');

            f1_map_cdf(f1_map > max(f1_x{kk})) = 1;
            f2_map_cdf(f1_map < min(f1_x{kk})) = 0;
            
%             figure('name', map_types{kk}); 
%             subplot(1,2,1); imagesc(f1_map); axis image;
%             subplot(1,2,2); imagesc(f2_map_cdf); axis image; caxis([0 1]);
            
%             figure('name', map_types{kk}); 
%             subplot(1,2,1); imagesc(f2_map); axis image;
%             subplot(1,2,2); imagesc(f2_map_cdf); axis image; caxis([0 1]);

            save_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk}(1:end-6) 'cdf' mam_types{kk}...
                '\' image_types{ll} '\' mam_names{ii} '_f1.mat'], f1_map_cdf);
            save_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk}(1:end-6) 'cdf' mam_types{kk}...
                '\' image_types{ll} '\' mam_names{ii} '_f2.mat'], f2_map_cdf);
        end
    end
end
%%
data_type = {'abnormals', 'normals'};
    
for ii = 1:2
    mammo_dir = ['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\' data_type{ii} '\'];

    mammo_list = dir([mammo_dir '*ML.mat']);
    mammo_names = get_mammo_info(mammo_list);
    
    for jj = 1:length(mammo_names);
        
        %load segmentation
        segmentation = u_load(['C:\isbe\asymmetry_project\data\segmentations\2004_screening\' data_type{ii} '\'...
            mammo_names{jj} '_segmentation.mat']);
        
        %load K-line map
        mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\' data_type{ii} '\'...
            mammo_names{jj} '_mask.mat']);
        [r c] = size(mask);
        
        %resize the segmentation to the map
        [breast_border breast_air pectoral_edge] = segment_breast_resize([r c], segmentation);

        %Make a mask of the pectoral triangle
        if strcmpi(mammo_names{jj}(4), 'R')
            sx = c;
            sy = 1;
        else
            sx = 1;
            sy = 1;
        end
        mask = mask & ~poly2mask([sx; pectoral_edge(:,1)], [sy; pectoral_edge(:,2)], r, c);

        %figure; imagesc(mask); axis image;
        save(['C:\isbe\asymmetry_project\data\masks_pectoral\2004_screening\' data_type{ii} '\'...
            mammo_names{jj} '_mask.mat'], 'mask');
    end
end
%%
thresh = [(-1:5)/10 (51:90)/100 (901:2:1001)/1000];
map_types = {'f1', 'f2'};
mam_type = {'normals', 'abnormals'};
for ii = 1:2
    for jj = 1:2
        [rf.tp rf.fp rf.fp_pixels] = mass_maps_roc(...
            'k_stellate_maps_rf_cdf', ['2004_screening/' mam_type{jj}], 'map_type', map_types{ii},...
            'mask_dir', 'masks_pectoral', 'thresh', thresh);

        save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rfs_froc_'...
            mam_type{jj} '_' map_types{ii} '.mat'], 'rf');
        
        [wrf.tp wrf.fp wrf.fp_pixels] = mass_maps_roc(...
            'k_stellate_maps_wrf_cdf', ['2004_screening/' mam_type{jj}], 'map_type', map_types{ii},...
            'mask_dir', 'masks_pectoral', 'thresh', thresh);

        save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrfs_froc_'...
            mam_type{jj} '_' map_types{ii} '.mat'], 'wrf');
        
        [g2.tp g2.fp g2.fp_pixels] = mass_maps_roc(...
            'k_stellate_maps_g2_cdf', ['2004_screening_processed/' mam_type{jj}], 'map_type', map_types{ii},...
            'mask_dir', 'masks_pectoral', 'thresh', thresh);

        save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2s_froc_'...
            mam_type{jj} '_' map_types{ii} '.mat'], 'g2');
    end
end
%%
thresh = [(-1:5)/10 (51:90)/100 (901:2:1001)/1000];

mam_type = {'normals', 'abnormals'};

for jj = 2
    for dd = [16 32 64 128 256]
%         [rf.tp rf.fp rf.fp_pixels] = mass_maps_roc(...
%             'radial_maps_rf_cdf', ['2004_screening/' mam_type{jj}], 'thresh', thresh, 'map_type', zerostr(dd,3));
%         save(['C:\isbe\asymmetry_project\experiments\radial_maps\rf_froc_'...
%             mam_type{jj} '_' zerostr(dd,3) '.mat'], 'rf');
% 
%         [wrf.tp wrf.fp wrf.fp_pixels] = mass_maps_roc(...
%             'radial_maps_wrf_cdf', ['2004_screening/' mam_type{jj}], 'thresh', thresh, 'map_type', zerostr(dd,3));
%         save(['C:\isbe\asymmetry_project\experiments\radial_maps\wrf_froc_'...
%             mam_type{jj} '_' zerostr(dd,3) '.mat'], 'wrf');
        
        [g2.tp g2.fp g2.fp_pixels] = mass_maps_roc(...
            'radial_maps_g2_cdf', ['2004_screening_processed/' mam_type{jj}], 'thresh', thresh, 'map_type', zerostr(dd,3), 'sigma', 4, 'resize_factor', 0.25);
        save(['C:\isbe\asymmetry_project\experiments\radial_maps\g2_froc_'...
            mam_type{jj} '_' zerostr(dd,3) '.mat'], 'g2');

        [wg2.tp wg2.fp wg2.fp_pixels] = mass_maps_roc(...
            'radial_maps_wg2_cdf', ['2004_screening_processed/' mam_type{jj}], 'thresh', thresh, 'map_type', zerostr(dd,3), 'sigma', 4, 'resize_factor', 0.25);
        save(['C:\isbe\asymmetry_project\experiments\radial_maps\wg2_froc_'...
            mam_type{jj} '_' zerostr(dd,3) '.mat'], 'wg2');
    
    end   
end
