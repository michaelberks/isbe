%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This script includes code to perform a classification experiment on our
% set of masses in the 2004 screening data set, based on the information
% contained in our radial maps
%--------------------------------------------------------------------------
%
% This experiment assumes for each mammogram we have already generated:
%   i) Line maps (using RF 191905 - 1x1 DT conj)
%   ii) Orientation maps (using RF 191934 - 3x3 DT conj)
%   iii) Relevance maps
%   iv) Weighted radial maps, using i-iii) for capture ranges of 8, 16, 32,
%    64, 128 and 256 (+ infinity calculation)
%   v) Mass template maps
%
% See also work101011.m for experiments computing how often local maxima in
% the mass template and radial maps lie within the borders of annotated
% masses and for computing thresholds on the maps
%
%--------------------------------------------------------------------------

% 1) Compute thresholds for the 95th percentile in each radial and template
% map (note this is an estimate of the true value across all data by
% averaging the 95th percentile from each image)
%--------------------------------------------------------------------------

%% Compute thresholds and hits for template maps
[template.thresh template.percentiles] = mass_maps_thresh('template_maps');
[template.hits template.num_maxima template.maxima_ranks template.maxima_scores] =...
    mass_maps_maxima('template_maps', template.thresh);
save('C:\isbe\asymmetry_project\experiments\radial_maps\template_hits.mat', 'template');
%--------------------------------------------------------------------------
%% Compute thresholds and hits for radial maps
[weighted_radial.thresh weighted_radial.percentiles] =...
    mass_maps_thresh('weighted_radial_maps', 'dist_ranges', [16 32 64 128 256], 'sigma', 8);
[weighted_radial.hits weighted_radial.num_maxima weighted_radial.maxima_ranks weighted_radial.maxima_scores] =...
    mass_maps_maxima('weighted_radial_maps', weighted_radial.thresh, 'dist_ranges', [16 32 64 128 256], 'sigma', 8);
save('C:\isbe\asymmetry_project\experiments\radial_maps\weighted_radial_hits.mat', 'weighted_radial');
%
[linop_radial.thresh linop_radial.percentiles] =...
    mass_maps_thresh('linop_radial_maps', 'dist_ranges', [16 32 64 128 256], 'sigma', 8);
[linop_radial.hits linop_radial.num_maxima linop_radial.maxima_ranks linop_radial.maxima_scores] =...
    mass_maps_maxima('linop_radial_maps', linop_radial.thresh, 'dist_ranges', [16 32 64 128 256], 'sigma', 8);
save('C:\isbe\asymmetry_project\experiments\radial_maps\linop_radial_hits.mat', 'linop_radial');

[combi_radial.thresh combi_radial.percentiles] =...
    mass_maps_thresh('combi_radial_maps', 'dist_ranges', [16 32 64 128 256], 'sigma', 8);
[combi_radial.hits combi_radial.num_maxima combi_radial.maxima_ranks combi_radial.maxima_scores] =...
    mass_maps_maxima('combi_radial_maps', combi_radial.thresh, 'dist_ranges', [16 32 64 128 256], 'sigma', 8);
save('C:\isbe\asymmetry_project\experiments\radial_maps\combi_radial_hits.mat', 'combi_radial');
%
[old_weighted_radial.thresh old_weighted_radial.percentiles] =...
    mass_maps_thresh('old_weighted_radial_maps', 'dist_ranges', [16 32 64 128 256], 'sigma', 8);
[old_weighted_radial.hits old_weighted_radial.num_maxima old_weighted_radial.maxima_ranks old_weighted_radial.maxima_scores] =...
    mass_maps_maxima('old_weighted_radial_maps', old_weighted_radial.thresh, 'dist_ranges', [16 32 64 128 256], 'sigma', 8);
save('C:\isbe\asymmetry_project\experiments\radial_maps\old_weighted_radial_hits.mat', 'old_weighted_radial');
%%
[radial.thresh radial.percentiles] =...
    mass_maps_thresh('radial_maps', 'dist_ranges', [16 32 64 128 256], 'sigma', 8);
[radial.hits radial.num_maxima radial.maxima_ranks radial.maxima_scores] =...
    mass_maps_maxima('radial_maps', radial.thresh, 'dist_ranges', [16 32 64 128 256], 'sigma', 8);
save('C:\isbe\asymmetry_project\experiments\radial_maps\radial_hits.mat', 'radial');
%%
[linop_radial2.thresh linop_radial2.percentiles] =...
    mass_maps_thresh('linop_radial_maps2', 'dist_ranges', [16 32 64 128 256], 'sigma', 8);
[linop_radial2.hits linop_radial2.num_maxima linop_radial2.maxima_ranks linop_radial2.maxima_scores] =...
    mass_maps_maxima('linop_radial_maps2', linop_radial2.thresh, 'dist_ranges', [16 32 64 128 256], 'sigma', 8);
save('C:\isbe\asymmetry_project\experiments\radial_maps\linop_radial2_hits.mat', 'linop_radial2');
%%
[mass_class.thresh mass_class.percentiles] =...
    mass_maps_thresh('mass_maps', 'sigma', 16);
[mass_class.hits mass_class.num_maxima mass_class.maxima_ranks mass_class.maxima_scores] =...
    mass_maps_maxima('mass_maps', mass_class.thresh, 'sigma', 16);
save('C:\isbe\asymmetry_project\experiments\radial_maps\mass_class_hits.mat', 'mass_class');
%%
[scaled.thresh scaled.percentiles] = mass_maps_thresh('scaled_wr_maps');
[scaled.hits scaled.num_maxima scaled.maxima_ranks scaled.maxima_scores] =...
    mass_maps_maxima('scaled_wr_maps', scaled.thresh);
save('C:\isbe\asymmetry_project\experiments\radial_maps\scaled_hits.mat', 'scaled');
%%
%**************************************************************************
%**************************************************************************
%%

% 2) Compute masks - based on mass template matches and the values in
% individual radial maps - to select pixels included in the mammograms
% - pixels must be above the 95th percentile in both the 64 AND 128 capture
%   range radial maps, OR be above the 95th percentile in the mass template
%   map
%
% - for masses, the points must also lie within the annotated border of the
%   mass
%--------------------------------------------------------------------------

%% Compute masks for abnormal images
mkdir C:\isbe\asymmetry_project\data\template_masks\2004_screening\abnormals\
spic_thresh = u_load('C:\isbe\asymmetry_project\experiments\radial_maps\weighted_radial_hits.mat');
mass_thresh = u_load('C:\isbe\asymmetry_project\experiments\radial_maps\template_hits.mat');

mammo_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\*.mat');
mammo_names = get_mammo_info(mammo_list);

mass_areas = zeros(146,1);

do_plot = 0;
for ii = 1:length(mammo_names)

    %load the mass and radial maps
    rad_map064 = u_load(['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\' mammo_names{ii} '_rad_map_064.mat']);
    rad_map128 = u_load(['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\' mammo_names{ii} '_rad_map_128.mat']);
    mass_map = u_load(['C:\isbe\asymmetry_project\data\template_maps\2004_screening\abnormals\' mammo_names{ii} '_template.mat']);
    %meta = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta2\' mammo_names{ii} '_meta.mat']);
    
    %Increase the local capture of the radial maps by Gaussian smoothing
    rad_map064 = imfilter(rad_map064, fspecial('gaussian', 40, 8));
    rad_map128 = imfilter(rad_map128, fspecial('gaussian', 40, 8));

    if do_plot
        figure; colormap(jet(256));
        subplot(1,2,1); imagesc(mass_map); axis image;
        subplot(1,2,2); imagesc(mass_map > mass_thresh.thresh); axis image;

        figure; colormap(jet(256));
        subplot(1,2,1); imagesc(rad_map064); axis image;
        subplot(1,2,2); imagesc(rad_map064 > spic_thresh.thresh(3)); axis image;

        figure; colormap(jet(256));
        subplot(1,2,1); imagesc(rad_map128); axis image;
        subplot(1,2,2); imagesc(rad_map128 > spic_thresh.thresh(4)); axis image;
    end

    %Compute the template mask
    template_mask = mass_map > mass_thresh.thresh;
    
    %Compute the radial maps mask
    rad_mask = (rad_map064 > spic_thresh.thresh(3)) & (rad_map128 > spic_thresh.thresh(4));
    
%     %Compute the mask for the mass border (may be multiple masses in the
%     %image)
%     mass_mask = false(size(rad_mask));
%     for kk = 1:length(meta)
%         meta{kk} = meta{kk} / 2;
%         mass_mask = mass_mask | roipoly(mass_mask, meta{kk}(:,1), meta{kk}(:,2)); 
%     end
% 
%     %Compute the final map and save
%     final_mask = mass_mask & (rad_mask | template_mask);
    final_mask = (rad_mask | template_mask);
    
    mass_areas(ii) = sum(final_mask(:));
    save(['C:\isbe\asymmetry_project\data\template_masks\2004_screening\abnormals\' mammo_names{ii} '_mask.mat'], 'final_mask');
end
%--------------------------------------------------------------------------
%% Compute masks for normal images
mkdir C:\isbe\asymmetry_project\data\mass_masks\2004_screening\normals\
spic_thresh = load('C:\isbe\asymmetry_project\experiments\radial_maps\thresh_weighted_radial_maps.mat', 'thresh');
mass_thresh = load('C:\isbe\asymmetry_project\experiments\radial_maps\thresh_mass_maps.mat', 'thresh');

mammo_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\normals\*.mat');
mammo_names = get_mammo_info(mammo_list);

normal_areas = zeros(146,1);

do_plot = 0;
for ii = 1:length(mammo_names)

    %load the template and radial maps
    rad_map064 = u_load(['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\normals\' mammo_names{ii} '_rad_map_064.mat']);
    rad_map128 = u_load(['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\normals\' mammo_names{ii} '_rad_map_128.mat']);
    mass_map = u_load(['C:\isbe\asymmetry_project\data\template_maps\2004_screening\normals\' mammo_names{ii} '_template.mat']);
    
    %Increase the local capture of the radial maps by Gaussian smoothing
    rad_map064 = imfilter(rad_map064, fspecial('gaussian', 40, 8));
    rad_map128 = imfilter(rad_map128, fspecial('gaussian', 40, 8));

    %Compute the template mask
    template_mask = mass_map > mass_thresh.thresh;
    
    %Compute the radial maps mask
    rad_mask = (rad_map064 > spic_thresh.thresh(4)) & (rad_map128 > spic_thresh.thresh(5));

    %Compute the final map and save
    final_mask = rad_mask | template_mask;
    
    normal_areas(ii) = sum(final_mask(:));
    save(['C:\isbe\asymmetry_project\data\mass_masks\2004_screening\normals\' mammo_names{ii} '_mask.mat'], 'final_mask');
end
%%
%--------------------------------------------------------------------------
%% Check data sampling methods are working
[training_data] = sample_mass_feature_data(...
    'num_samples', 1e4,...
    'mask_dir', 'C:\isbe\asymmetry_project\data\mass_masks\2004_screening\abnormals\',...
    'radial_dir', 'C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\',...
    'template_dir', 'C:\isbe\asymmetry_project\data\template_maps\2004_screening\abnormals\',...
    'dist_range', [32 64 128, 256],...
    'sigma_range', [1 2 4 8],...
    'angular_res', 1,...
    'mammo_list', [], ...
    'image_type', '.mat',...
    'save_path', []);
%%
[training_data training_labels] = sample_mass_training_data(... % non-strict mode
    'num_samples', 2e3,...
    'abnormal_data', '2004_screening/abnormals/',...
    'normal_data', '2004_screening/normals/',... 
    'image_dir', [asymmetryroot, 'data/mammograms/'],...
    'mask_dir', [asymmetryroot, 'data/mass_masks/'],...
    'radial_dir', [asymmetryroot, 'data/weighted_radial_maps2/'],...
    'template_dir', [asymmetryroot, 'data/template_maps/'],...
    'fold_id', 1,...
    'num_folds', 10,...
    'view', [],...
    'dist_range', [32 64 128 256],...
    'sigma_range', [1 2 4 8],...
    'angular_res', 12,...
    'save_path', []);
%%
%--------------------------------------------------------------------------
%% Build random forest classifier on the cluster
%% Classify each mammogram using the Random Forests
%%
%--------------------------------------------------------------------------
%% ------------------ Results and analysis --------------------------------
%--------------------------------------------------------------------------
