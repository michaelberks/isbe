function auto_stats = apex_stats_to_xls(auto_stats, varargin)
%ANALYSE_EXTRACTED_APEX_MEASURES *Insert a one line summary here*
%   [] = extract_vessel_centres()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 18-Jun-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
args = u_packargs(varargin,... % the user's input
    0, ... % non-strict mode
    {'image_names'},         ...
    'data_dir',             [nailfoldroot 'data/2_year_study/'],...
    'vessel_centre_dir',    'vessel_centres\',...
    'selected_dir',         [],...
    'metrics_dir',          'apex_metrics',...
    'selected_features', [],...
    'save_dir',         [],...
    'um_per_pix',       1.25,...
    'xls_filename',     'auto_stats.xls');

vessel_centre_dir = [args.data_dir '/' args.vessel_centre_dir '/'];
selected_dir = [args.data_dir '/' args.selected_dir '/'];
metrics_dir = [args.data_dir '/' args.metrics_dir '/'];

if isempty(args.selected_features)
    selected_features = {...
        'num_distal_vessels',...
        'num_nondistal_vessels',...
        'median_apex_width',...
        'mean_mean_width',...
        'mean_weighted_width',...
        'median_weighted_width',...
        'mean_median_width',...
        'mean_std_width',...
        'max_mean_width',...
        'mean_max_width',...
        'mean_min_width',...
        'max_max_width',...
        'std_mean_width',...
        'mean_orientation_entropy',...
        'median_orientation_entropy',...
        'std_orientation_entropy',...
        'mean_connected_orientation_dispersion',...
        'mean_weighted_orientation_dispersion',...
        'dispersion_base_orientation',...
        'dispersion_connected_orientation',...
        'dispersion_weighted_orientation',...
        'total_vessel_prob',...
        'mean_vessel_prob',...
        'total_scores',...
        'mean_scores',...
        'score_density',...
        'vessel_density1',...
        'vessel_density2',...
        'mean_inter_capillary_distance',...
        'std_inter_capillary_distance',...
        'median_inter_capillary_distance'};
else
    selected_features = args.selected_features;
end

if ~exist('auto_stats', 'var') || isempty(auto_stats)
    args.do_auto_stats = 1;
end

if args.do_auto_stats
        
    num_images = length(args.image_names);    

    auto_stats.image_name = args.image_names(:);

    auto_stats.gradeable = false(num_images,1);
    auto_stats.status = zeros(num_images,1);

    auto_stats.num_distal_vessels = zeros(num_images,1);
    auto_stats.num_nondistal_vessels = zeros(num_images,1);

    auto_stats.median_apex_width = nan(num_images, 1);
    auto_stats.mean_mean_width = nan(num_images, 1);
    auto_stats.mean_weighted_width = nan(num_images, 1);
    auto_stats.median_weighted_width = nan(num_images, 1);
    auto_stats.mean_median_width = nan(num_images, 1);
    auto_stats.mean_std_width = nan(num_images, 1);
    auto_stats.max_mean_width = nan(num_images, 1);
    auto_stats.mean_max_width = nan(num_images, 1);
    auto_stats.max_max_width = nan(num_images, 1);
    auto_stats.mean_min_width = nan(num_images, 1);
    auto_stats.std_mean_width = nan(num_images, 1);   

    auto_stats.total_vessel_prob = zeros(num_images,1);
    auto_stats.mean_vessel_prob = zeros(num_images,1);

    auto_stats.total_scores = zeros(num_images,1);
    auto_stats.mean_scores = zeros(num_images,1);

    auto_stats.mean_orientation_entropy = nan(num_images, 1);
    auto_stats.median_orientation_entropy = nan(num_images, 1);
    auto_stats.std_orientation_entropy = nan(num_images, 1);
    
    auto_stats.mean_connected_orientation_dispersion = nan(num_images, 1);
    auto_stats.mean_weighted_orientation_dispersion = nan(num_images, 1);
    
    auto_stats.dispersion_base_orientation = nan(num_images, 1);
    auto_stats.dispersion_connected_orientation = nan(num_images, 1);
    auto_stats.dispersion_weighted_orientation = nan(num_images, 1);

    auto_stats.vessel_density1 = nan(num_images, 1);
    auto_stats.vessel_density2 = nan(num_images, 1);

    auto_stats.mean_inter_capillary_distance = nan(num_images, 1);
    auto_stats.std_inter_capillary_distance = nan(num_images, 1);
    auto_stats.median_inter_capillary_distance = nan(num_images, 1);
        
    auto_stats.score_density = nan(num_images,1);

    for i_im = 1:num_images
        im_name = args.image_names{i_im};
        
        load([metrics_dir im_name '_am.mat']);
        s = load([selected_dir im_name '_sel.mat']);

        auto_stats.gradeable(i_im) = ~isempty(s.selected_distal);
        auto_stats.num_distal_vessels(i_im) = sum(s.selected_distal);
        auto_stats.num_nondistal_vessels(i_im) = sum(s.selected_non_distal); 
        auto_stats.status(i_im) = s.status;

        if auto_stats.num_distal_vessels(i_im)
            auto_stats.median_apex_width(i_im) = naNmedian(apex_measures.distal.width_at_apex) * args.um_per_pix;
            auto_stats.mean_mean_width(i_im) = naNmean(apex_measures.distal.mean_width) * args.um_per_pix;
            auto_stats.mean_weighted_width(i_im) = naNmean(apex_measures.distal.mean_weighted_width) * args.um_per_pix;
            auto_stats.median_weighted_width(i_im) = naNmedian(apex_measures.distal.mean_weighted_width) * args.um_per_pix;
            auto_stats.mean_median_width(i_im) = naNmean(apex_measures.distal.median_width) * args.um_per_pix;
            auto_stats.mean_std_width(i_im) = naNmean(apex_measures.distal.std_width) * args.um_per_pix;
            auto_stats.max_mean_width(i_im) = nanmax(apex_measures.distal.mean_weighted_width) * args.um_per_pix;          
            auto_stats.mean_max_width(i_im) = naNmean(apex_measures.distal.max_width) * args.um_per_pix;
            auto_stats.max_max_width(i_im) = nanmax(apex_measures.distal.max_width) * args.um_per_pix;
            auto_stats.mean_min_width(i_im) = naNmean(apex_measures.distal.min_width) * args.um_per_pix;
            auto_stats.std_mean_width(i_im) = naNstd(apex_measures.distal.mean_weighted_width) * args.um_per_pix;


            ori_entropy = mb_entropy(apex_measures.distal.orientation_hist,2);
            auto_stats.mean_orientation_entropy(i_im) = naNmean(ori_entropy);
            auto_stats.median_orientation_entropy(i_im) = naNmedian(ori_entropy);
            auto_stats.std_orientation_entropy(i_im) = naNstd(ori_entropy);
            
            auto_stats.mean_connected_orientation_dispersion(i_im) = naNmean(abs(apex_measures.distal.connected_orientation));
            auto_stats.mean_weighted_orientation_dispersion(i_im) = naNmean(abs(apex_measures.distal.weighted_orientation));
            
            auto_stats.dispersion_base_orientation(i_im) = abs(naNmean(exp(1i*angle(apex_measures.distal.base_orientation))));
            auto_stats.dispersion_connected_orientation(i_im) = abs(naNmean(exp(1i*angle(apex_measures.distal.connected_orientation))));
            auto_stats.dispersion_weighted_orientation(i_im) = abs(naNmean(exp(1i*angle(apex_measures.distal.weighted_orientation))));

            auto_stats.total_vessel_prob(i_im) = naNsum(apex_measures.distal.total_prob);
            auto_stats.mean_vessel_prob(i_im) = naNmean(apex_measures.distal.total_prob);

            auto_stats.total_scores(i_im) = naNsum(apex_measures.distal.candidate_scores);
            auto_stats.mean_scores(i_im) = naNmean(apex_measures.distal.candidate_scores);
            
            
            %Compute density measures
            load([vessel_centre_dir im_name '_vc.mat'], 'ncols');
            d1 = (ncols - 1) * args.um_per_pix / 1000;
            
            apex_xy = sortrows(apex_measures.distal.apex_xy);
            
            d2 = sqrt(sum((apex_xy(1,:)-apex_xy(end,:)).^2)) * args.um_per_pix / 1000;
            inter_d = sqrt(sum(diff(apex_xy).^2,2)) * args.um_per_pix / 1000;

            auto_stats.vessel_density1(i_im) = auto_stats.num_distal_vessels(i_im) / d1;
            auto_stats.vessel_density2(i_im) = auto_stats.num_distal_vessels(i_im) / d2;

            auto_stats.mean_inter_capillary_distance(i_im) = mean(inter_d);
            auto_stats.std_inter_capillary_distance(i_im) = std(inter_d);
            auto_stats.median_inter_capillary_distance(i_im) = median(inter_d);
            
            auto_stats.score_density(i_im) = auto_stats.total_scores(i_im) / d1;

        end
    end
    
    if ~isempty(args.save_dir)
        create_folder(args.save_dir);
        save([args.save_dir '\auto_stats.mat'], 'auto_stats');
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------     
xls_data = cell(num_images+1, length(selected_features)+1);
xls_data{1,1} = 'Image names';
xls_data(2:end,1) = args.image_names;
for i_f = 1:length(selected_features);

    feature = selected_features{i_f};
    feature_str = feature;
    feature_str(feature_str == '_') = ' ';
    feature_str(1) = feature_str(1) - 32;

    xls_data{1, i_f + 1} = feature_str;

    col_data = auto_stats.(feature);

    if ~iscell(col_data)
        col_data = num2cell(col_data);
    end
    xls_data(2:end,i_f+1) = col_data;
end
xlswrite(args.xls_filename, xls_data, 1, 'A1');            
