%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%*********************  KARSSEMEIJER MAPS *********************************
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
%%
%--------------------------------------------------------------------------
%% 1. Combine the maps across scale for each image
map_type = {'g2d_all', 'g2d_1', 'g2d_2', 'g2d_3', 'g2d_rf_all', 'g2d_rf_1', 'g2d_rf_2'};
ab_type = {'normals', 'abnormals'};
offset_r = 0;
offset_c = 0;
template_scales = [30 45 60 75 90];
spacing = 4;
do_template = 0;
do_max = 1;
%
for aa = 1:2 %Compute for both normals and normals
    
    %load in the lists of mammograms names
    mam_names = u_load(['C:\isbe\asymmetry_project\data\mam_names\2004_screening_' ab_type{aa}]);
    
    for mm = 6 %For each method rf, wrf, g2
        
        %create a directory for the combined maps
        mkdir(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
            '_template_scale\2004_screening_processed\' ab_type{aa} '\']);
        mkdir(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
            '_max_scale\2004_screening_processed\' ab_type{aa} '\']);
        
        for ii = 1:length(mam_names) %for each mammogram

            %if ii == 66; continue; end
            %load in the f1 and f2 maps
            f1 = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
                '\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_f1.mat']);
            f2 = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
                '\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_f2.mat']);
            mask = u_load(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
                '\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_mask.mat']);
            mask = mask(offset_r + (1:spacing:end), offset_c + (1:spacing:end));
            
            %Use sigma = 16.2mm as default base for r-max
            if 0%ii < 11
                f1_map = zeros(size(mask));
                f2_map = zeros(size(mask));
                f1_map(mask) = f1(:,3);
                f2_map(mask) = f2(:,3);   
            end
            if do_template
                f1_template_scale = f1(:,3);
                f2_template_scale = f2(:,3);

                %load in the maps giving the template blob scores and the
                %associated blob radius
                scale_map = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening_processed\'...
                    ab_type{aa} '\scales\' mam_names{ii} '_scales.mat']);
                template_map = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening_processed\'...
                    ab_type{aa} '\' mam_names{ii} '_template.mat']);

                if any(size(template_map)~=size(mask))
                    template_map = imresize(template_map, size(mask), 'bilinear');
                    scale_map = imresize(scale_map, size(mask), 'nearest');
                end

                %threshold the blob map
                template_map = template_map > 0.3;

                for jj = [1 2 4 5] %for each r_max value (other than 120)
                    %Swap any pixels where there is evidence of a blob at this
                    %scale
                    swap_mask = (scale_map == template_scales(jj)) & template_map & mask;
                    swap_idx = swap_mask(mask);
                    if ii < 11
                        f1_map(swap_mask) = f1(swap_idx,jj);
                        f2_map(swap_mask) = f2(swap_idx,jj);
                    end

                    f1_template_scale(swap_idx,:) = f1(swap_idx,jj);
                    f2_template_scale(swap_idx,:) = f2(swap_idx,jj);

                end
                if 0%ii < 11
                    figure; 
                    subplot(1,2,1); imagesc(f1_map); axis image;
                    subplot(1,2,2); imagesc(f2_map); axis image;
                end

                %Save the scaled maps
                save_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
                    '_template_scale\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_f1.mat'],...
                    f1_template_scale);
                save_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
                    '_template_scale\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_f2.mat'],...
                    f2_template_scale);

                %Copy over the mask
                copyfile(...
                    ['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
                    '\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_mask.mat'],...
                    ['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
                    '_template_scale\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_mask.mat']);
            end
            if do_max
            
                %Also create maps based on the maximum score over scale
                f1_max_scale = max(f1,[],2);
                f2_max_scale = max(f2,[],2);

                 %Save the scaled maps
                save_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
                    '_max_scale\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_f1.mat'],...
                    f1_max_scale);
                save_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
                    '_max_scale\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_f2.mat'],...
                    f2_max_scale);

                %Copy over the mask
                copyfile(...
                    ['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
                    '\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_mask.mat'],...
                    ['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
                    '_max_scale\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_mask.mat']);
            end
        end
    end
end

%--------------------------------------------------------------------------
%% 2. For the scaled f1 and f2 maps, compute the distribution of scores
% across the dataset for each method
map_type = {'g2d_all', 'g2d_1', 'g2d_2', 'g2d_3', 'g2d_rf_all', 'g2d_rf_1', 'g2d_rf_2'};
scale_type = {'_max_scale', '_template_scale'};
f_type = {'f1', 'f2'};
for mm = 6, for ss=1, for ff=1:2, %#ok

    [p_x x] = mass_maps_hist(...
        '2004_screening_processed/abnormals',...
        ['k_stellate_maps/' map_type{mm} scale_type{ss}],...
        'view', f_type{ff}, 'sigma', 0, 'mask_dir', []);
    
    %Discard x values with no p counts
    empty_x = ~p_x;
    x(empty_x) = [];
    p_x(empty_x) = [];
    P_x = cumsum(p_x) / sum(p_x);
    save(['C:\isbe\asymmetry_project\experiments\k_maps\cdf_'...
        map_type{mm} scale_type{ss} '_' f_type{ff} '.mat'], 'P_x', 'p_x', 'x');

        end
     end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% 3.a. Compute FROC scores for each map
map_type = {'g2d_all', 'g2d_1', 'g2d_2', 'g2d_3', 'g2d_rf_all'};
scale_type = {'_max_scale', '_template_scale'};
f_type = {'f1', 'f2'};
ab_type = {'normals', 'abnormals'};
for mm =5, for ss=1:2, for ff=1:2, for aa = 1 %#ok

    load(['C:\isbe\asymmetry_project\experiments\k_maps\cdf_'...
        map_type{mm} scale_type{ss} '_' f_type{ff} '.mat'], 'P_x', 'x');
    
    thresh = interp1(P_x, x, linspace(0,1,101), 'linear');
    thresh(1) = min(x); thresh(end) = max(x);
     [tp fp fp_pixels] = k_maps_roc(...
        ['k_stellate_maps/' map_type{mm} scale_type{ss}],...
        ['2004_screening_processed/' ab_type{aa}], 'feature_type', f_type{ff},...
        'thresh', thresh, 'plot', 0);%, 'num_mammos', 20
    save(['C:\isbe\asymmetry_project\experiments\k_maps\froc_'...
        map_type{mm} scale_type{ss} '_' f_type{ff} '_' ab_type{aa} '.mat'], 'tp', 'fp', 'fp_pixels');

            end
        end
     end
end
%--------------------------------------------------------------------------
%% 3.b. Compute FROC scores for each map - focus on high end specificity
clear;
map_type = {'g2d_all', 'g2d_1', 'g2d_2', 'g2d_3', 'g2d_rf_all', 'g2d_rf_1', 'g2d_rf_2'};
scale_type = {'_max_scale', '_template_scale'};
f_type = {'f1', 'f2'};
ab_type = {'normals', 'abnormals'};
spacing = 4;

for mm = 6, for ss=1, for ff=1:2, for aa = 1:2 %#ok

    load(['C:\isbe\asymmetry_project\experiments\k_maps\cdf_'...
        map_type{mm} scale_type{ss} '_' f_type{ff} '.mat'], 'P_x', 'x');
    
    thresh = interp1(P_x, x, linspace(0.9,1,50), 'linear');
    thresh(1) = min(x); thresh(end) = max(x);
     [tp fp fp_pixels] = k_maps_roc(...
        ['k_stellate_maps/' map_type{mm} scale_type{ss}],...
        ['2004_screening_processed/' ab_type{aa}], 'feature_type', f_type{ff},...
        'thresh', thresh, 'plot', 0, 'spacing', spacing, 'sigma', 16/spacing);%, 'num_mammos', 20
    save(['C:\isbe\asymmetry_project\experiments\k_maps\froc_'...
        map_type{mm} scale_type{ss} '_' f_type{ff} '_' ab_type{aa} '_h.mat'], 'tp', 'fp', 'fp_pixels');

            end
        end
     end
end
%--------------------------------------------------------------------------
%% 4. Convert map scores into probabilities by interpolating against the CDF of each map type
map_type = {'g2d_all', 'g2d_1', 'g2d_2', 'g2d_3', 'g2d_rf_all'};
scale_type = {'_max_scale', '_template_scale'};
f_type = {'f1', 'f2'};
ab_type = {'normals', 'abnormals'};

for aa = 1:2, for mm = 1:5, for ss=1, for ff=1:2 %#ok
                
    
    mam_names = u_load(['C:\isbe\asymmetry_project\data\mam_names\2004_screening_' ab_type{aa}]);
    
    %Setup dir of existing maps
    map_dir =  ['C:\isbe\asymmetry_project\data\k_stellate_maps\'...
        map_type{mm}  scale_type{ss} '\2004_screening_processed\' ab_type{aa} '\'];
    
    %Make directory for new cdf maps
    cdf_dir =  ['C:\isbe\asymmetry_project\data\k_stellate_maps\'...
        map_type{mm}  scale_type{ss} '_cdf\2004_screening_processed\' ab_type{aa} '\'];
    mkdir(cdf_dir);
    
    %Load cdf info
    load(['C:\isbe\asymmetry_project\experiments\k_maps\cdf_'...
        map_type{mm} scale_type{ss} '_' f_type{ff} '.mat'], 'P_x', 'x');
    for ii = 1:length(mam_names)
        %load map
        f_map = load_uint8([map_dir mam_names{ii} '_' f_type{ff}]);
        f_map_cdf = interp1(x, P_x, f_map, 'linear');
        f_map_cdf(f_map <= min(x)) = 0;
        f_map_cdf(f_map >= max(x)) = 1;
        
        %save the new map
        save_uint8([cdf_dir mam_names{ii} '_' f_type{ff}], f_map_cdf);
        
        %copy across mask
        copyfile([map_dir mam_names{ii} '_mask.mat'],[cdf_dir mam_names{ii} '_mask.mat']);
    end
            end
        end
    end
end     
%--------------------------------------------------------------------------
%% 5. Compute FROC curves for sum and product of CDF maps
map_type = {'g2d_all', 'g2d_1', 'g2d_2', 'g2d_3', 'g2d_rf_all'};
scale_type = {'_max_scale', '_template_scale'};
ab_type = {'normals', 'abnormals'};
spacing = 4;

for mm = 1:5, for ss=1, for aa = 1:2 %#ok
    
    thresh = linspace(0.8,1,100);
    
    %f1*f2
    [tp fp fp_pixels] = k_maps_roc_andor(...
        ['k_stellate_maps/' map_type{mm} scale_type{ss} '_cdf'],...
        ['2004_screening_processed/' ab_type{aa}],...
        'thresh', thresh, 'plot', 0, 'and_or', 1, 'spacing', spacing, 'sigma', 16/spacing);%, 'num_mammos', 20
    save(['C:\isbe\asymmetry_project\experiments\k_maps\froc_'...
        map_type{mm} scale_type{ss} '_' ab_type{aa} '_cdf_hp.mat'], 'tp', 'fp', 'fp_pixels');
    %(f1+f2)/2
    [tp fp fp_pixels] = k_maps_roc_andor(...
        ['k_stellate_maps/' map_type{mm} scale_type{ss} '_cdf'],...
        ['2004_screening_processed/' ab_type{aa}],...
        'thresh', thresh, 'plot', 0, 'and_or', 0, 'spacing', spacing, 'sigma', 16/spacing);%, 'num_mammos', 20
    save(['C:\isbe\asymmetry_project\experiments\k_maps\froc_'...
        map_type{mm} scale_type{ss} '_' ab_type{aa} '_cdf_hs.mat'], 'tp', 'fp', 'fp_pixels');
        end
     end
end
%%
% %Try doing it for mixed maps too (F1 from g2d, F2 from rf)
% for ss=2, for aa = 1:2 %#ok
%     
%     thresh = linspace(0.8,1,100);
%     
%     %f1*f2
%     [tp fp fp_pixels] = k_maps_roc_andor(...
%         ['k_stellate_maps/g2d' scale_type{ss} '_cdf'],...
%         ['2004_screening_processed/' ab_type{aa}],...
%         'f2_map_dir', ['k_stellate_maps/rf_thin' scale_type{ss} '_cdf'],...
%         'thresh', thresh, 'plot', 0, 'and_or', 1);%, 'num_mammos', 20
%     save(['C:\isbe\asymmetry_project\experiments\k_maps\froc_grf'...
%         scale_type{ss} '_' ab_type{aa} '_cdf_hp.mat'], 'tp', 'fp', 'fp_pixels');
%     %(f1+f2)/2
%     [tp fp fp_pixels] = k_maps_roc_andor(...
%         ['k_stellate_maps/g2d' scale_type{ss} '_cdf'],...
%         ['2004_screening_processed/' ab_type{aa}],...
%         'f2_map_dir', ['k_stellate_maps/rf_thin' scale_type{ss} '_cdf'],...
%         'thresh', thresh, 'plot', 0, 'and_or', 0);%, 'num_mammos', 20
%     save(['C:\isbe\asymmetry_project\experiments\k_maps\froc_grf'...
%         scale_type{ss} '_' ab_type{aa} '_cdf_hs.mat'], 'tp', 'fp', 'fp_pixels');
% 
%      end
% end
%--------------------------------------------------------------------------
%% 6. Display the FROC curves
%% Individual maps
circle_area = pi*(20*50/9)^2;
map_type = {'g2d_all', 'g2d_1', 'g2d_2', 'g2d_3', 'g2d_rf_all', 'rf_thin', 'g2d', 'g2d_rf_1'};
scale_type = {'_max_scale', '_template_scale'};
f_type = {'f1', 'f2'};
f_title = {'F1', 'F2'};
figure;
for ff=1:2
    subplot(1,2,ff); hold all; title(f_title{ff});
    legend_text = [];
    for mm = 1:8 
        for ss=1

    fp = load(['C:\isbe\asymmetry_project\experiments\k_maps\froc_'...
        map_type{mm} scale_type{ss} '_' f_type{ff} '_normals_h.mat'], 'fp', 'fp_pixels');
    tp = load(['C:\isbe\asymmetry_project\experiments\k_maps\froc_'...
        map_type{mm} scale_type{ss} '_' f_type{ff} '_abnormals_h.mat'], 'tp');
    
    plot(mean(fp.fp(:,2:end)), mean(tp.tp(:,2:end) > 0), '-x');
    legend_text{end+1} = [map_type{mm} scale_type{ss}];%#ok
    legend_text{end}(legend_text{end}=='_') = ' ';   %#ok     
        end
    end
    axis([0 20 0 1]);
    legend(legend_text, 'location', 'southeast');
    
end 
%--------------------------------------------------------------------------
%% Combined F1 and F2 maps
map_type = {'g2d_all', 'g2d_1', 'g2d_2', 'g2d_3', 'g2d_rf_all'};
scale_type = {'_max_scale', '_template_scale'};
f_type = {'p', 's'};
f_title = {'F1 * F2', 'F1 + F2'};
figure;
for ff=1:2
    subplot(1,2,ff); hold all; title(f_title{ff});
    legend_text = [];
    for mm = 1:5, 
        for ss=1

    fp = load(['C:\isbe\asymmetry_project\experiments\k_maps\froc_'...
        map_type{mm} scale_type{ss} '_normals_cdf_h' f_type{ff} '.mat'], 'fp', 'fp_pixels');
    tp = load(['C:\isbe\asymmetry_project\experiments\k_maps\froc_'...
        map_type{mm} scale_type{ss} '_abnormals_cdf_h' f_type{ff} '.mat'], 'tp');
    
    plot(mean(fp.fp(:,2:end)), mean(tp.tp(:,2:end) > 0), '-x');
    legend_text{end+1} = [map_type{mm} scale_type{ss}];%#ok
    legend_text{end}(legend_text{end}=='_') = ' ';   %#ok     
        end
    end
    axis([0 20 0 1]);
    legend(legend_text, 'location', 'southeast');
end 
%--------------------------------------------------------------------------
%% Combined F1 (g2d) and F2 (rf) maps
scale_type = {'_max_scale', '_template_scale'};
f_type = {'p', 's'};
f_title = {'F1 * F2', 'F1 + F2'};

for ss=1:2,  %#ok
        
    figure; 
    for ff = 1:2
        subplot(1,2,ff); hold all;  title(f_title{ff});
        fp = load(['C:\isbe\asymmetry_project\experiments\k_maps\froc_g2d'...
            scale_type{ss} '_normals_cdf_h' f_type{ff} '.mat'], 'fp', 'fp_pixels');
        tp = load(['C:\isbe\asymmetry_project\experiments\k_maps\froc_g2d'...
            scale_type{ss} '_abnormals_cdf_h' f_type{ff} '.mat'], 'tp');

        plot(mean(fp.fp(:,2:end)), mean(tp.tp(:,2:end) > 0), '-x');
        
        subplot(1,2,ff); hold all;  title(f_title{ff});
        fp = load(['C:\isbe\asymmetry_project\experiments\k_maps\froc_rf_prob'...
            scale_type{ss} '_normals_cdf_h' f_type{ff} '.mat'], 'fp', 'fp_pixels');
        tp = load(['C:\isbe\asymmetry_project\experiments\k_maps\froc_rf_prob'...
            scale_type{ss} '_abnormals_cdf_h' f_type{ff} '.mat'], 'tp');

        plot(mean(fp.fp(:,2:end)), mean(tp.tp(:,2:end) > 0), '-x');

        fp = load(['C:\isbe\asymmetry_project\experiments\k_maps\froc_grf'...
            scale_type{ss} '_normals_cdf_h' f_type{ff} '.mat'], 'fp', 'fp_pixels');
        tp = load(['C:\isbe\asymmetry_project\experiments\k_maps\froc_grf'...
            scale_type{ss} '_abnormals_cdf_h' f_type{ff} '.mat'], 'tp');

        plot(mean(fp.fp(:,2:end)), mean(tp.tp(:,2:end) > 0), '-x');
        axis([0 20 0 1]);
        legend({'g2d', 'rf prob', 'grf'}, 'location', 'southeast');

     end
end
%%
map_type = {'g2d_all', 'g2d_1', 'g2d_2', 'g2d_3', 'g2d_rf_all'};
ab_type = {'normals', 'abnormals'};
%
for aa = 1:2 %Compute for both normals and normals
    
    %load in the lists of mammograms names
    mam_names = u_load(['C:\isbe\asymmetry_project\data\mam_names\2004_screening_' ab_type{aa}]);
    
    for mm = [1 4 6 7] %For each method rf, wrf, g2
        if mm > 5
            spacing = 4;
        else
            spacing = 2;
        end
        
        mkdir(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
            '_max_scale\2004_screening_processed\' ab_type{aa} '\']);
        
        for ii = 1:length(mam_names) %for each mammogram

            if mm == 7 && aa == 2 && ii == 66; continue; end
            
            %load in the f1 and f2 maps
            f1 = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
                '\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_f1.mat']);
            f2 = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
                '\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_f2.mat']);
            
            %Also create maps based on the maximum score over scale
            [dummy f1_max_scale] = max(f1,[],2);
            [dummy f2_max_scale] = max(f2,[],2);
            
             %Save the scaled maps
            save_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
                '_max_scale\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_f1_scale.mat'],...
                f1_max_scale);
            save_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
                '_max_scale\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_f2_scale.mat'],...
                f2_max_scale);
        end
    end
end
%%
