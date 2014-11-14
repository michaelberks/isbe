%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%**************************  RADIAL MAPS **********************************
%--------------------------------------------------------------------------
%
% This script contains code relevant to radial maps that highlight
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
% We then use the function compute_rad_maps_batch to build the radial maps
% on Hydra, using distance sigmas of 60, 90, 120, 150 and 180.
%
% Maps are saved in the following dirs:
%   radial_maps_rf - RF/DT-CWT maps, no weighting
%   radial_maps_wrf - RF/DT-CWT maps, relevance weighted
%   radial_maps_g2 - gaussian 2nd derivative maps, no weighting
%   radial_maps_wg2 - gaussian 2nd derivative maps, relevance weighted
%%
%-------------------------------------------------------------------------
% 1. Compute the distribution of map scores for each distance 
% map and method
for dist = [16 32 64 128 256] 
    [P_x x] = mass_maps_hist('2004_screening/abnormals', 'radial_maps_wrf', 'view', zerostr(dist,3), 'sigma', 8); 
    save(['C:\isbe\asymmetry_project\experiments\radial_maps\cdf_wrf_' zerostr(dist,3) '.mat'], 'P_x', 'x');
    
    [P_x x] = mass_maps_hist('2004_screening/abnormals', 'radial_maps_rf', 'view', zerostr(dist,3), 'sigma', 8);
    save(['C:\isbe\asymmetry_project\experiments\radial_maps\cdf_rf_' zerostr(dist,3) '.mat'], 'P_x', 'x');
    
    [P_x x] = mass_maps_hist('2004_screening_processed/abnormals', 'radial_maps_g2', 'view', zerostr(dist,3), 'sigma', 8);
    save(['C:\isbe\asymmetry_project\experiments\radial_maps\cdf_g2_' zerostr(dist,3) '.mat'], 'P_x', 'x');
    
    [P_x x] = mass_maps_hist('2004_screening_processed/abnormals', 'radial_maps_wg2', 'view', zerostr(dist,3), 'sigma', 8);
    save(['C:\isbe\asymmetry_project\experiments\radial_maps\cdf_wg2_' zerostr(dist,3) '.mat'], 'P_x', 'x');
end
%%
%-------------------------------------------------------------------------
%% 2. Convert the distribution of map scores into CDFs
dists = [16 32 64 128 256];

%Create a structure to hold the feature score ranges and their associated
%CDFs
f_x = cell(4,5);
f_P_x = cell(4,5);

for jj = 1:5
    %For each distance, load in the distribution of map scores, then take
    %the cumulative sum to compute the CDF.
    dist = dists(jj);
    
    %RF/DT-CWT no weighting
    load(['C:\isbe\asymmetry_project\experiments\radial_maps\cdf_rf_' zerostr(dist,3) '.mat']);
    f_x{1, jj} = x;
    f_P_x{1, jj} = cumsum(P_x) / sum(P_x);

    %RF/DT-CWT relevance weighted
    load(['C:\isbe\asymmetry_project\experiments\radial_maps\cdf_wrf_' zerostr(dist,3) '.mat']);
    f_x{2, jj} = x;
    f_P_x{2, jj} = cumsum(P_x) / sum(P_x);

    %Gaussian no weighting
    load(['C:\isbe\asymmetry_project\experiments\radial_maps\cdf_g2_' zerostr(dist,3) '.mat']);
    f_x{3, jj} = x;
    f_P_x{3, jj} = cumsum(P_x) / sum(P_x);

    %Gaussian relevance weighted
    load(['C:\isbe\asymmetry_project\experiments\radial_maps\cdf_wg2_' zerostr(dist,3) '.mat']);
    f_x{4, jj} = x;
    f_P_x{4, jj} = cumsum(P_x) / sum(P_x);
    clear x P_x
    
end
%--------------------------------------------------------------------------
%% 3. Create a copy of each map with the appropriate CDF probability substituted for the origanl feature score
image_types = {'abnormals', 'normals'};
map_types = {'radial_maps_rf', 'radial_maps_wrf', 'radial_maps_g2', 'radial_maps_wg2'};
mam_types = {'\2004_screening', '\2004_screening', '\2004_screening_processed', '\2004_screening_processed'};

for ll = 1:2
    %
    mam_names = u_load(['C:\isbe\asymmetry_project\data\mam_names\2004_screening_' image_types{ll}]);
    for kk = 1:4
        
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
%--------------------------------------------------------------------------
%% 4. Compute FROC scores for each map
thresh = [(-1:5)/10 (51:90)/100 (901:2:1001)/1000];

mam_type = {'normals', 'abnormals'};

for jj = 2
    for dd = [16 32 64 128 256]
        [rf.tp rf.fp rf.fp_pixels] = mass_maps_roc(...
            'radial_maps_rf_cdf', ['2004_screening/' mam_type{jj}], 'thresh', thresh, 'map_type', zerostr(dd,3), 'sigma', 0);
        save(['C:\isbe\asymmetry_project\experiments\radial_maps\rf_froc_'...
            mam_type{jj} '_' zerostr(dd,3) '.mat'], 'rf');

        [wrf.tp wrf.fp wrf.fp_pixels] = mass_maps_roc(...
            'radial_maps_wrf_cdf', ['2004_screening/' mam_type{jj}], 'thresh', thresh, 'map_type', zerostr(dd,3), 'sigma', 0);
        save(['C:\isbe\asymmetry_project\experiments\radial_maps\wrf_froc_'...
            mam_type{jj} '_' zerostr(dd,3) '.mat'], 'wrf');
        
        [g2.tp g2.fp g2.fp_pixels] = mass_maps_roc(...
            'radial_maps_g2_cdf', ['2004_screening_processed/' mam_type{jj}], 'thresh', thresh, 'map_type', zerostr(dd,3), 'sigma', 0);
        save(['C:\isbe\asymmetry_project\experiments\radial_maps\g2_froc_'...
            mam_type{jj} '_' zerostr(dd,3) '.mat'], 'g2');

        [wg2.tp wg2.fp wg2.fp_pixels] = mass_maps_roc(...
            'radial_maps_wg2_cdf', ['2004_screening_processed/' mam_type{jj}], 'thresh', thresh, 'map_type', zerostr(dd,3), 'sigma', 0);
        save(['C:\isbe\asymmetry_project\experiments\radial_maps\wg2_froc_'...
            mam_type{jj} '_' zerostr(dd,3) '.mat'], 'wg2');
    
    end   
end
%--------------------------------------------------------------------------
%% 5. Display the FROC curves
%
circle_area = pi*(20*50/9)^2;

for dd = [16 32 64 128 256]
    
    rf_ab = u_load(['C:\isbe\asymmetry_project\experiments\radial_maps\rf_froc_abnormals_' zerostr(dd,3) '.mat']);
    rf_norm = u_load(['C:\isbe\asymmetry_project\experiments\radial_maps\rf_froc_normals_' zerostr(dd,3) '.mat']);

    wrf_ab = u_load(['C:\isbe\asymmetry_project\experiments\radial_maps\wrf_froc_abnormals_' zerostr(dd,3) '.mat']);
    wrf_norm = u_load(['C:\isbe\asymmetry_project\experiments\radial_maps\wrf_froc_normals_' zerostr(dd,3) '.mat']);
    
    froc_rf = compute_froc_curve(rf_ab.tp, rf_norm.fp, rf_norm.fp_pixels, circle_area);
    froc_wrf = compute_froc_curve(wrf_ab.tp, wrf_norm.fp, wrf_norm.fp_pixels, circle_area);

%     g2_ab = u_load(['C:\isbe\asymmetry_project\experiments\radial_maps\g2_froc_abnormals_' zerostr(dd,3) '.mat']);
%     g2_norm = u_load(['C:\isbe\asymmetry_project\experiments\radial_maps\g2_froc_normals_' zerostr(dd,3) '.mat']);
%     
%     wg2_ab = u_load(['C:\isbe\asymmetry_project\experiments\radial_maps\wg2_froc_abnormals_' zerostr(dd,3) '.mat']);
%     wg2_norm = u_load(['C:\isbe\asymmetry_project\experiments\radial_maps\wg2_froc_normals_' zerostr(dd,3) '.mat']);

    figure; hold all;
%     plot(mean(rf_norm.fp_pixels) / circle_area, mean(rf_ab.tp > 0), 'x');
%     plot(mean(wrf_norm.fp_pixels) / circle_area, mean(wrf_ab.tp > 0), 'x');
    plot(mean(rf_norm.fp), mean(rf_ab.tp > 0), '-x');
    plot(mean(wrf_norm.fp), mean(wrf_ab.tp > 0), '-x');
    plot(froc_rf(:,1), froc_rf(:,2), '-x');
    plot(froc_wrf(:,1), froc_wrf(:,2), '-x');
%     plot(mean(g2_norm.fp_pixels) / circle_area, mean(g2_ab.tp > 0), 'x');
%     plot(mean(wg2_norm.fp_pixels) / circle_area, mean(wg2_ab.tp > 0), 'x');
%     legend({'RF maps', 'weighted RF maps', 'Gaussian maps', 'weighted Gaussian maps'});
    legend({'RF maps', 'weighted RF maps', 'RF maps', 'weighted RF maps'});
    xlabel('False positives per image')
    ylabel('Sensitivity');
    title(['FROC curve for radial maps, distance \sigma =  ' num2str(dd)]);
end