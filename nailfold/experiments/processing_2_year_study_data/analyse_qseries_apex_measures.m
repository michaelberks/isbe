function auto_stats = analyse_qseries_apex_measures(auto_stats, varargin)
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
    'do_xls',           1,...
    'do_auto_stats',    0,...
    'do_people_stats',    0,...
    'do_image_plots',   0, ...
    'do_people_plots',  0,...
    'fig_dir',          [],...
    'save_dir',         [],...
    'um_per_pix',       1.6077,...
    'xls_filename_images',     'image_auto_stats.xls',...
    'xls_filename_subjects',     'subject_auto_stats.xls');

metrics_dir = [args.data_dir '/' args.metrics_dir '/'];

if isempty(args.selected_features)
    selected_features = {...
        'num_distal_vessels',...
        'num_nondistal_vessels',...
        'num_giant_vessels',...
        'num_enlarged_vessels',...
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
        'vessel_density',...
        'mean_inter_capillary_distance',...
        'std_inter_capillary_distance',...
        'median_inter_capillary_distance'};
else
    selected_features = args.selected_features;
end

if ~exist('auto_stats', 'var') || isempty(auto_stats)
    args.do_make_stats = 1;
end

if args.do_auto_stats
        
    num_images = length(args.image_names);    

    auto_stats.category = zeros(num_images,1);
    auto_stats.super_category = zeros(num_images,1);
    auto_stats.hand = zeros(num_images,1);
    auto_stats.visit = zeros(num_images,1);
    auto_stats.digit = zeros(num_images,1);
    auto_stats.people_id = zeros(num_images,1);
    auto_stats.people_id_str = cell(num_images,1);

    auto_stats.gradeable = false(num_images,1);
    auto_stats.status = zeros(num_images,1);

    auto_stats.num_distal_vessels = zeros(num_images,1);
    auto_stats.num_nondistal_vessels = zeros(num_images,1);
    auto_stats.num_giant_vessels = zeros(num_images,1);
    auto_stats.num_enlarged_vessels = zeros(num_images,1);

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

    auto_stats.vessel_density = nan(num_images, 1);

    auto_stats.mean_inter_capillary_distance = nan(num_images, 1);
    auto_stats.std_inter_capillary_distance = nan(num_images, 1);
    auto_stats.median_inter_capillary_distance = nan(num_images, 1);
        
    auto_stats.score_density = nan(num_images,1);

    for i_im = 1:num_images
        im_name = args.image_names{i_im};

        switch im_name(1)
            
            case 'c'
                auto_stats.category(i_im) = 1;
                auto_stats.super_category(i_im) = 1;
            case 'd'
                auto_stats.category(i_im) = 2;
                auto_stats.super_category(i_im) = 2;
            case 'l'
                auto_stats.category(i_im) = 3;
                auto_stats.super_category(i_im) = 2;
            case 'n'
                auto_stats.category(i_im) = 6;
                auto_stats.super_category(i_im) = 1;
            case 'p'
                auto_stats.category(i_im) = 4;
                auto_stats.super_category(i_im) = 3;
            case 'u'
                auto_stats.category(i_im) = 5;
                auto_stats.super_category(i_im) = 0;
            case 'D'
                auto_stats.category(i_im) = 2;
                auto_stats.super_category(i_im) = 2;
            case 'L'
                auto_stats.category(i_im) = 3;
                auto_stats.super_category(i_im) = 2;
            case 'f'
                auto_stats.category(i_im) = 7;
                auto_stats.super_category(i_im) = 4;
            case 'F'
                auto_stats.category(i_im) = 7;
                auto_stats.super_category(i_im) = 4;
            otherwise
                continue;               
        end        
        %Try a 3 digit id
        people_id = str2double(im_name(2:5));
        v_idx = 6;
        if isnan(people_id)
            %Try a 3 digit id
            people_id = str2double(im_name(2:4));
            v_idx = 5;
            if isnan(people_id)
                %Try a 2 digit id
                people_id = str2double(im_name(2:3));
                v_idx = 4;
                if isnan(people_id)
                    %Try a 1 digit id
                    people_id = str2double(im_name(2));
                    v_idx = 3;
                end
            end
        end
        auto_stats.people_id_str{i_im} = im_name(1:v_idx-1);
        auto_stats.people_id(i_im) = people_id + 1000*auto_stats.category(i_im);
        auto_stats.visit(i_im) = str2double(im_name(v_idx+1));
        auto_stats.digit(i_im) = str2double(im_name(v_idx+4));
        
        switch im_name(v_idx+2)
            case 'L'
                auto_stats.hand(i_im) = 1;
            case 'R'
                auto_stats.hand(i_im) = 2;
        end
        
        %load([apex_gt_dir  im_name '_gt.mat']);
        load([metrics_dir im_name '_am.mat']);
        
        auto_stats.gradeable(i_im) = ~isempty(apex_measures.distal); 
        if ~auto_stats.gradeable(i_im)
            continue;
        end

        auto_stats.num_distal_vessels(i_im) = sum(~isnan(apex_measures.distal.mean_weighted_width));
        auto_stats.num_nondistal_vessels(i_im) = sum(~isnan(apex_measures.nondistal.mean_weighted_width));

        auto_stats.status(i_im) = auto_stats.num_distal_vessels(i_im) < 3;
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
            
            auto_stats.num_giant_vessels(i_im) = sum((apex_measures.distal.mean_weighted_width* args.um_per_pix) > 50);
            auto_stats.num_enlarged_vessels(i_im) = sum((apex_measures.distal.mean_weighted_width* args.um_per_pix) > 20);

%             if auto_stats.num_nondistal_vessels(i_im)
%                 auto_stats.num_giant_vessels(i_im) = auto_stats.num_giant_vessels(i_im) +...
%                     sum((apex_measures.nondistal.mean_weighted_width* args.um_per_pix) > 50);
%                 auto_stats.num_enlarged_vessels(i_im) = auto_stats.num_enlarged_vessels(i_im) +...
%                     sum((apex_measures.nondistal.mean_weighted_width* args.um_per_pix) > 20);
%             end
            auto_stats.num_enlarged_vessels(i_im) = auto_stats.num_enlarged_vessels(i_im) -...
                auto_stats.num_giant_vessels(i_im);

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
            if auto_stats.num_distal_vessels(i_im) > 2
                apex_xy = sortrows(apex_measures.distal.apex_xy(~isnan(apex_measures.distal.width_at_apex),:));
                d2 = sqrt(sum((apex_xy(1,:)-apex_xy(end,:)).^2)) * args.um_per_pix / 1000;
                auto_stats.vessel_density(i_im) = auto_stats.num_distal_vessels(i_im) / d2;
                auto_stats.score_density(i_im) = auto_stats.total_scores(i_im) / d2;

                inter_d = sqrt(sum(diff(apex_xy).^2,2)) * args.um_per_pix / 1000;
                auto_stats.mean_inter_capillary_distance(i_im) = mean(inter_d);
                auto_stats.std_inter_capillary_distance(i_im) = std(inter_d);
                auto_stats.median_inter_capillary_distance(i_im) = median(inter_d);


            end

        end
    end
    
    if ~isempty(args.save_dir)
        save([args.save_dir '\auto_stats.mat'], 'auto_stats');
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
if args.do_xls
    if ~exist('num_images', 'var')
        num_images = length(args.image_names);
    end
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
    xlswrite(args.xls_filename_images, xls_data, 1, 'A1');
        
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
if args.do_people_stats
    [people_stats.people_ids u_idx] = unique(auto_stats.people_id);
    people_stats.people_ids_str = auto_stats.people_id_str(u_idx);
    num_people = length(people_stats.people_ids);
    num_images = length(auto_stats.category);
    people_stats.category = zeros(num_people,1);
    people_stats.super_category = zeros(num_people,1);

    people_stats.present = false(num_people,5,2,2);

    people_stats.gradeable = false(num_people,5,2,2);
    people_stats.status = zeros(num_people,5,2,2);
    
    for i_f = 1:length(selected_features);

        feature = selected_features{i_f};
        people_stats.(feature) = zeros(num_people,5,2,2);
    end

    for i_im = 1:num_images

        p_id = auto_stats.people_id(i_im);

        if ~p_id
            continue;
        end
        p = people_stats.people_ids==p_id;
        
        v = auto_stats.visit(i_im);
        d = auto_stats.digit(i_im);
        h = auto_stats.hand(i_im);

        people_stats.present(p, d, h, v) = 1;
        people_stats.category(p) = auto_stats.category(i_im);
        people_stats.super_category(p) = auto_stats.super_category(i_im);
        people_stats.gradeable(p, d, h, v) = auto_stats.gradeable(i_im);
        people_stats.status(p, d, h, v) = auto_stats.status(i_im);
        
        for i_f = 1:length(selected_features);
            feature = selected_features{i_f};
            people_stats.(feature)(p, d, h, v) = auto_stats.(feature)(i_im);
        end

    end
    
    if ~isempty(args.save_dir)
        save([args.save_dir '\people_stats.mat'], 'people_stats');
    end
    
    if args.do_xls
        xls_data = cell(num_people+2, 3*length(selected_features)+1);
        xls_data{2,1} = 'Subject IDs';
        
        xls_data(3:end,1) = people_stats.people_ids_str;
        for i_f = 1:length(selected_features);

            feature = selected_features{i_f};
            feature_str = feature;
            feature_str(feature_str == '_') = ' ';
            feature_str(1) = feature_str(1) - 32;

            xls_data{1, 3*i_f - 1} = feature_str;
            xls_data{2, 3*i_f - 1} = '0 mths';
            xls_data{2, 3*i_f} = '12 mths';
            xls_data{2, 3*i_f + 1} = '24 mths';

            col_data = squeeze(sum(sum(people_stats.(feature),2),3));

            if ~iscell(col_data)
                col_data = num2cell(col_data(:,[1,5,6]));
            end
            xls_data(3:end, 3*i_f+(-1:1) ) = col_data;
        end
        xlswrite(args.xls_filename_subjects, xls_data, 1, 'A1');
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
if args.do_image_plots
    
    for i_f = 1:length(selected_features);

            feature = selected_features{i_f};
            feature_str = feature;
            feature_str(feature_str == '_') = ' ';
            feature_str(1) = feature_str(1) - 32;

            dist_ss = ...
                auto_stats.(feature)(auto_stats.super_category == 2 & ~auto_stats.status);
            dist_ss(isnan(dist_ss) | ~(dist_ss < inf)) = [];
            dist_pr = ...
                auto_stats.(feature)(auto_stats.super_category == 3 & ~auto_stats.status);                
            dist_pr(isnan(dist_pr) | ~(dist_pr < inf)) = [];
            dist_hc = ...
                auto_stats.(feature)(auto_stats.super_category == 1 & ~auto_stats.status);
            dist_hc(isnan(dist_hc) | ~(dist_hc < inf)) = [];

            p_ss_hc = ranksum(dist_ss, dist_hc);
            p_ss_pr = ranksum(dist_ss, dist_pr);
            p_hc_pr = ranksum(dist_hc, dist_pr);
            
            if p_ss_hc < 0.0001
                p_ss_hc_str = '< 0.0001';
            else
                p_ss_hc_str = ['= ' num2str(p_ss_hc,2)];
            end
            if p_ss_pr < 0.0001
                p_ss_pr_str = '< 0.0001';
            else
                p_ss_pr_str = ['= ' num2str(p_ss_pr,2)];
            end
            if p_hc_pr < 0.0001
                p_hc_pr_str = '< 0.0001';
            else
                p_hc_pr_str = ['= ' num2str(p_hc_pr,2)];
            end

            display(['Feature: ' feature]);
            display(['SS vs HC: p = ' num2str(p_ss_hc)]);
            display(['SS vs PR: p = ' num2str(p_ss_pr)]);
            display(['PR vs HC: p = ' num2str(p_hc_pr)]);
            display('');

            min_val = min([dist_ss; dist_pr; dist_hc]);
            max_val = max([dist_ss; dist_pr; dist_hc]);
            grid_pts = linspace(min_val, max_val, 100);

            kdist_ss = build_1d_kernel_distribution(dist_ss, grid_pts, 0);
            kdist_pr = build_1d_kernel_distribution(dist_pr, grid_pts, 0);
            kdist_hc = build_1d_kernel_distribution(dist_hc, grid_pts, 0);

            if ~isempty(args.fig_dir)
                figure;
                hold all;
                title({['Feature: \bf' feature_str ', \rm'];...
                    ['SS vs HC: p ' p_ss_hc_str ...
                    ', SS vs PR: p ' p_ss_pr_str ...
                    ', PR vs HC: p ' p_hc_pr_str]},...
                    'fontsize', 18);
                set(gca, 'fontSize', 18);
                plot(kdist_ss.x, kdist_ss.D_f, 'linewidth', 2);
                plot(kdist_pr.x, kdist_pr.D_f, 'linewidth', 2);
                plot(kdist_hc.x, kdist_hc.D_f, 'linewidth', 2);
                legend({'SSc', 'PR', 'HC'});
                saveas(gcf, [args.fig_dir feature '_pdf.fig']);
                exportfig([args.fig_dir feature '_pdf.png']);
            else
                figure; 
                subplot(12,2,1:2);
                title({['Feature: \bf' feature_str ', \rm images as independent samples'];...
                    ['SS vs HC: p = ' num2str(p_ss_hc,3) ...
                    ', SS vs PR: p = ' num2str(p_ss_pr,3) ...
                    ', PR vs HC: p = ' num2str(p_hc_pr,3)]});
                axis off;
                subplot(12,2,3:2:23); hold all;
                title('Kernel estimated PDF');
                plot(kdist_ss.x, kdist_ss.D_f, 'linewidth', 2);
                plot(kdist_pr.x, kdist_pr.D_f, 'linewidth', 2);
                plot(kdist_hc.x, kdist_hc.D_f, 'linewidth', 2);
                legend({'SSc', 'PR', 'HC'});

                subplot(12,2,4:2:24); hold all;
                title('Kernel estimated CDF');

                plot(prctile(dist_ss, 0:100), (0:100)/100, 'linewidth', 2);
                plot(prctile(dist_pr, 0:100), (0:100)/100, 'linewidth', 2);
                plot(prctile(dist_hc, 0:100), (0:100)/100, 'linewidth', 2);
                set(gca, 'ylim', [0 1]);
            end

    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
if args.do_people_plots
       
    ss_idx = people_stats.super_category == 2;
    pr_idx = people_stats.super_category == 3;
    hc_idx = people_stats.super_category == 1;

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
                    'vessel_density',...
      'mean_inter_capillary_distance',...
       'std_inter_capillary_distance',...
    'median_inter_capillary_distance'};
    else
        selected_features = args.selected_features;
    end
    for i_f = 1:length(selected_features);

            feature = selected_features{i_f};
            feature_str = feature;
            feature_str(feature_str == '_') = ' ';
            feature_str(1) = feature_str(1) - 32;

            display('***********************');
            display(['Feature: ' feature]);

            dist_f = people_stats.(feature);
            dist_f = dist_f(:,:);
            valid_ims = sum(people_stats.present(:,:) & ~people_stats.status(:,:)...
                & ~isnan(dist_f) & (dist_f < inf), 2);

            max_f = max(dist_f, [], 2);
            mean_f = naNsum(dist_f, 2) ./  valid_ims;

            p_ss_hc_max = ranksum(max_f(ss_idx), max_f(hc_idx));
            p_ss_pr_max = ranksum(max_f(ss_idx), max_f(pr_idx));
            p_hc_pr_max = ranksum(max_f(hc_idx), max_f(pr_idx));

            p_ss_hc_avg = ranksum(mean_f(ss_idx), mean_f(hc_idx));
            p_ss_pr_avg = ranksum(mean_f(ss_idx), mean_f(pr_idx));
            p_hc_pr_avg = ranksum(mean_f(hc_idx), mean_f(pr_idx));

            display(['SS vs HC: p = ' num2str(p_ss_hc_max)]);
            display(['SS vs PR: p = ' num2str(p_ss_pr_max)]);
            display(['PR vs HC: p = ' num2str(p_hc_pr_max)]);
            display('');

            display(['SS vs HC: p = ' num2str(p_ss_hc_avg)]);
            display(['SS vs PR: p = ' num2str(p_ss_pr_avg)]);
            display(['PR vs HC: p = ' num2str(p_hc_pr_avg)]);
            display('');

            min_val = min(max_f);
            max_val = max(max_f);
            grid_pts = linspace(min_val, max_val, 100);

            kdist_ss = build_1d_kernel_distribution(max_f(ss_idx), grid_pts, 0);
            kdist_pr = build_1d_kernel_distribution(max_f(pr_idx), grid_pts, 0);
            kdist_hc = build_1d_kernel_distribution(max_f(hc_idx), grid_pts, 0);

            figure; 
            subplot(10,2,1:2);
            title({['Feature: \bf' feature_str];...
                ['\bf Max over images: \rm SS vs HC: p = ' num2str(p_ss_hc_max,3) ...
                ', SS vs PR: p = ' num2str(p_ss_pr_max,3) ...
                ', PR vs HC: p = ' num2str(p_hc_pr_max,3)];...
                ['\bf Mean over images: \rm SS vs HC: p = ' num2str(p_ss_hc_avg,3) ...
                ', SS vs PR: p = ' num2str(p_ss_pr_avg,3) ...
                ', PR vs HC: p = ' num2str(p_hc_pr_avg,3)]});
            axis off;

            subplot(10,2,3:2:9); hold all;
            title('Kernel estimated PDF');
            plot(kdist_ss.x, kdist_ss.D_a, 'linewidth', 2);
            plot(kdist_pr.x, kdist_pr.D_a, 'linewidth', 2);
            plot(kdist_hc.x, kdist_hc.D_a, 'linewidth', 2);
            legend({'SSc', 'PR', 'HC'});
            ylabel('Max of all images per subject');

            subplot(10,2,4:2:10); hold all;
            title('Kernel estimated CDF');
            plot(kdist_ss.x, cumsum(kdist_ss.D_a), 'linewidth', 2);
            plot(kdist_pr.x, cumsum(kdist_pr.D_a), 'linewidth', 2);
            plot(kdist_hc.x, cumsum(kdist_hc.D_a), 'linewidth', 2);
            set(gca, 'ylim', [0 1]);       

            %------------------------------------------------------
            min_val = min(mean_f);
            max_val = max(mean_f);
            grid_pts = linspace(min_val, max_val, 100);

            kdist_ss = build_1d_kernel_distribution(mean_f(ss_idx), grid_pts, 0);
            kdist_pr = build_1d_kernel_distribution(mean_f(pr_idx), grid_pts, 0);
            kdist_hc = build_1d_kernel_distribution(mean_f(hc_idx), grid_pts, 0);

            subplot(10,2,13:2:19); hold all;
            plot(kdist_ss.x, kdist_ss.D_a, 'linewidth', 2);
            plot(kdist_pr.x, kdist_pr.D_a, 'linewidth', 2);
            plot(kdist_hc.x, kdist_hc.D_a, 'linewidth', 2);
            ylabel('Mean of all images per subject');

            subplot(10,2,14:2:20); hold all;
            plot(kdist_ss.x, cumsum(kdist_ss.D_a), 'linewidth', 2);
            plot(kdist_pr.x, cumsum(kdist_pr.D_a), 'linewidth', 2);
            plot(kdist_hc.x, cumsum(kdist_hc.D_a), 'linewidth', 2);
            set(gca, 'ylim', [0 1]);

    end

end
