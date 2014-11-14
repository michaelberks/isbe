%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%********************  KARSSEMEIJER MAPS **********************************
%--------------------------------------------------------------------------
%
% This script contains code relevant to Karssemeijer maps that highlight
% locations in a mammogram likely to be the centre of a radial patterns of
% structures (and thus possibly a spiculated lesion)
%
% We assume that for a set of abnormal and normal mammograms we have
% computed maps of line strength and line orientation. In this script we
% test 2 methods for obtaining these maps - our random forest/DT-CWT method
% and a method using Gaussian 2nd derivative filters (the method originally
% used by Karssemeijer). We also test weighting the line strength maps with
% an abnormality relevance measure
%
% We then use the function compute_k_maps_batch to build karssemeijer maps
% on Hydra, using r_max spacings of 60, 90, 120, 150 and 180. For each
% two feature maps, labelled f1 and f2, are produced.
% Maps are saved in the following dirs:
%   k_stellate_maps_rf - RF/DT-CWT maps, no weighting
%   k_stellate_maps_wrf - RF/DT-CWT maps, relevance weighted
%   k_stellate_maps_g2 - gaussian 2nd derivative maps, no weighting
%
% Having built the maps at each r_max scale we then combine the maps across
% scale into a single map. To do this we use precomputed map of blobs
% in the images (computed using template matching). Where there is evidence
% of a strong blob at some radius R, we take the feature score from the map
% with r_max = 2R, otherwise we use as default the features in the maps
% r_max = 120
%%
%-------------------------------------------------------------------------
%% 1. Combine the maps across scale for each image
map_types = {'k_stellate_maps_rf', 'k_stellate_maps_wrf', 'k_stellate_maps_g2'};
mam_types = {'\2004_screening', '\2004_screening', '\2004_screening_processed'};
ab_types = {'normals', 'abnormals'};
%
for ll = 1:2 %Compute for both normals and normals
    
    %load in the lists of mammograms names
    mam_names = u_load(['C:\isbe\asymmetry_project\data\mam_names\2004_screening_' ab_types{ll}]);
    
    for kk = 1:3 %For each method rf, wrf, g2
        
        %create a directory for the combined maps
        mkdir(['C:\isbe\asymmetry_project\data\' map_types{kk} '_scaled' mam_types{kk}...
            '\' ab_types{ll} '\']);
        
        for ii = 1:length(mam_names) %for each mammogram

            %load in the f1 and f2 maps for r_max = 120 as the deafult base
            %maps
            f1_map = load_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} mam_types{kk}...
                '\' ab_types{ll} '\' mam_names{ii} '_f1_120.mat']);
            f2_map = load_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} mam_types{kk}...
                '\' ab_types{ll} '\' mam_names{ii} '_f2_120.mat']);

            %load in the maps giving the template blob scores and the
            %associated blob radius
            scale_map = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening_processed\' ab_types{ll} '\scales\'...
                mam_names{ii} '_scales.mat']);
            template_map = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening_processed\' ab_types{ll} '\'...
                mam_names{ii} '_template.mat']);
            
            %threshold the blob map
            template_map = template_map > 0.3;

            for jj = [30 45 75 90] %for each r_max value (other than 120)
                
                %Load in the f1 and f2 map for this r_max value
                f1_j = load_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} mam_types{kk}...
                    '\' ab_types{ll} '\' mam_names{ii} '_f1_' zerostr(2*jj,3) '.mat']);
                f2_j = load_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} mam_types{kk}...
                    '\' ab_types{ll} '\' mam_names{ii} '_f2_' zerostr(2*jj,3) '.mat']);

                %Swap any pixels where there is evidence of a blob at this
                %scale
                swap_idx = (scale_map == jj) & template_map;

                f1_map(swap_idx) = f1_j(swap_idx);
                f2_map(swap_idx) = f2_j(swap_idx);

            end
%             figure; 
%             subplot(1,2,1); imagesc(f1_map); axis image;
%             subplot(1,2,2); imagesc(f2_map); axis image;

            %Save the scaled maps
            save_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} '_scaled' mam_types{kk}...
                '\' ab_types{ll} '\' mam_names{ii} '_f1.mat'], f1_map);
            save_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} '_scaled' mam_types{kk}...
                '\' ab_types{ll} '\' mam_names{ii} '_f2.mat'], f2_map);
        end
    end
end

%--------------------------------------------------------------------------
%% 2. For the scaled f1 and f2 maps, compute the distribution of scores
% across the dataset for each method
map_types = {'f1', 'f2'};
for ii =1:2

    [P_x x] = mass_maps_hist('2004_screening/normals', 'k_stellate_maps_wrf_scaled', 'view', map_types{ii}, 'sigma', 8); 
    save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\cdf_wrf_' map_types{ii} '_scaled.mat'], 'P_x', 'x');

    [P_x x] = mass_maps_hist('2004_screening/normals', 'k_stellate_maps_rf_scaled', 'view', map_types{ii}, 'sigma', 8);
    save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\cdf_rf_' map_types{ii} '_scaled.mat'], 'P_x', 'x');

    [P_x x] = mass_maps_hist('2004_screening_processed/normals', 'k_stellate_maps_g2_scaled', 'view', map_types{ii}, 'sigma', 8);
    save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\cdf_g2_' map_types{ii} '_scaled.mat'], 'P_x', 'x');

end
%-------------------------------------------------------------------------
%% 3. Convert the distribution of map scores into CDFs

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

%--------------------------------------------------------------------------
%% 3. Convert map scores into probabilities by interpolating from the CDF
% for each feature map
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
%--------------------------------------------------------------------------
%% 4. Compute FROC scores for each map
thresh = [(-1:5)/10 (51:90)/100 (901:2:1001)/1000];
map_types = {'f1', 'f2'};
mam_type = {'normals', 'abnormals'};
for ii = 1:2
    for jj = 1:2
        
        %RF/DT-CWT, no weighting
        [rf.tp rf.fp rf.fp_pixels] = mass_maps_roc(...
            'k_stellate_maps_rf_cdf', ['2004_screening/' mam_type{jj}], 'map_type', map_types{ii},...
            'mask_dir', 'masks_pectoral', 'thresh', thresh);

        save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rfs_froc_'...
            mam_type{jj} '_' map_types{ii} '.mat'], 'rf');
        
        %RF/DT-CWT, relevance weighted
        [wrf.tp wrf.fp wrf.fp_pixels] = mass_maps_roc(...
            'k_stellate_maps_wrf_cdf', ['2004_screening/' mam_type{jj}], 'map_type', map_types{ii},...
            'mask_dir', 'masks_pectoral', 'thresh', thresh);

        save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrfs_froc_'...
            mam_type{jj} '_' map_types{ii} '.mat'], 'wrf');
        
        %Gaussian, no weighting
        [g2.tp g2.fp g2.fp_pixels] = mass_maps_roc(...
            'k_stellate_maps_g2_cdf', ['2004_screening_processed/' mam_type{jj}], 'map_type', map_types{ii},...
            'mask_dir', 'masks_pectoral', 'thresh', thresh);

        save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2s_froc_'...
            mam_type{jj} '_' map_types{ii} '.mat'], 'g2');
    end
end
%--------------------------------------------------------------------------
%% 5. Display the FROC curves
%
circle_area = pi*(20*50/9)^2;

rf_f1_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rfs_froc_abnormals_f1.mat');
rf_f1_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rfs_froc_normals_f1.mat');

wrf_f1_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrfs_froc_abnormals_f1.mat');
wrf_f1_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrfs_froc_normals_f1.mat');

g2_f1_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2s_froc_abnormals_f1.mat');
g2_f1_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2s_froc_normals_f1.mat');

figure; hold all;
plot(mean(rf_f1_norm.fp), mean(rf_f1_ab.tp > 0), 'x');
plot(mean(wrf_f1_norm.fp), mean(wrf_f1_ab.tp > 0), 'x');
plot(mean(g2_f1_norm.fp), mean(g2_f1_ab.tp > 0), 'x');
legend({'RF maps', 'weighted RF maps', 'gaussian maps'});
%
rf_f2_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rfs_froc_abnormals_f2.mat');
rf_f2_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rfs_froc_normals_f2.mat');

wrf_f2_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrfs_froc_abnormals_f2.mat');
wrf_f2_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrfs_froc_normals_f2.mat');

g2_f2_ab = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2s_froc_abnormals_f2.mat');
g2_f2_norm = u_load('C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2s_froc_normals_f2.mat');

figure; hold all;
plot(mean(rf_f2_norm.fp), mean(rf_f2_ab.tp > 0), 'x');
plot(mean(wrf_f2_norm.fp), mean(wrf_f2_ab.tp > 0), 'x');
plot(mean(g2_f2_norm.fp), mean(g2_f2_ab.tp > 0), 'x');
legend({'RF maps', 'weighted RF maps', 'gaussian maps'});