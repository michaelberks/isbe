function auto_stats = analyse_extracted_apex_measures(auto_stats, varargin)
%EXTRACT_VESSEL_CENTRES *Insert a one line summary here*
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
    'do_make_stats',    1,...
    'do_image_plots',   0, ...
    'do_people_plots',  0);

if ~exist('auto_stats', 'var') || isempty(auto_stats)
    args.do_make_stats = 1;
end

if args.do_make_stats
    load C:\isbe\nailfold\data\rsa_study\image_id_data.mat
    image_dir = 'C:\isbe\nailfold\data\rsa_study\images\';
    version_dir = {'mixed_maxima_new_merged\', 'mixed_maxima_new\'};
    cans_dir = {'test_candidates2\', 'selected_candidates\'};
    data_dir = {'test_half', 'final_test'};
    for i_test = 1:2

        test_dir = data_dir{i_test};

        %apex_gt_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_gt\'];
        vessel_centre_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\vessel_centres\full_centres\'];
        apex_measures_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_metrics\'  version_dir{i_test}];
        selected_apexes_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_maps\set12g_half_296655\' version_dir{i_test}  cans_dir{i_test}];

        %im_list = dir([apex_gt_dir '*_gt.mat']);
        im_list = dir([apex_measures_dir '*_am.mat']);

        num_images = length(im_list);    

        auto_stats.(test_dir).category = zeros(num_images,1); 
        auto_stats.(test_dir).hand = zeros(num_images,1);
        auto_stats.(test_dir).visit = zeros(num_images,1);
        auto_stats.(test_dir).digit = zeros(num_images,1);
        auto_stats.(test_dir).people_id = zeros(num_images,1);

        auto_stats.(test_dir).gradeable = false(num_images,1);
        auto_stats.(test_dir).status = zeros(num_images,1);

        auto_stats.(test_dir).num_distal_vessels = zeros(num_images,1);
        auto_stats.(test_dir).num_nondistal_vessels = zeros(num_images,1);

        auto_stats.(test_dir).mean_mean_width = nan(num_images, 1);
        auto_stats.(test_dir).mean_weighted_width = nan(num_images, 1);
        auto_stats.(test_dir).median_weighted_width = nan(num_images, 1);
        auto_stats.(test_dir).mean_median_width = nan(num_images, 1);
        auto_stats.(test_dir).mean_std_width = nan(num_images, 1);
        auto_stats.(test_dir).max_mean_width = nan(num_images, 1);
        auto_stats.(test_dir).mean_max_width = nan(num_images, 1);
        auto_stats.(test_dir).max_max_width = nan(num_images, 1);
        auto_stats.(test_dir).mean_min_width = nan(num_images, 1);
        auto_stats.(test_dir).std_mean_width = nan(num_images, 1);   

        auto_stats.(test_dir).total_vessel_prob = zeros(num_images,1);
        auto_stats.(test_dir).mean_vessel_prob = zeros(num_images,1);

        auto_stats.(test_dir).total_scores = zeros(num_images,1);
        auto_stats.(test_dir).mean_scores = zeros(num_images,1);

        auto_stats.(test_dir).mean_orientation_entropy = nan(num_images, 1);
        auto_stats.(test_dir).median_orientation_entropy = nan(num_images, 1);
        auto_stats.(test_dir).std_orientation_entropy = nan(num_images, 1);

        auto_stats.(test_dir).vessel_density = nan(num_images,1);
        auto_stats.(test_dir).score_density = nan(num_images,1);

        for i_im = 1:num_images
            im_name = im_list(i_im).name(1:6);

            im_idx = find(strcmp(image_id_data.im_names, im_name));
            if isempty(im_idx)
                continue;
            end
            switch image_id_data.category{im_idx}
                case 'S'
                    auto_stats.(test_dir).category(i_im) = 1;
                case 'P'
                    auto_stats.(test_dir).category(i_im) = 2;
                case 'HC'
                    auto_stats.(test_dir).category(i_im) = 3;
                otherwise
                    continue;               
            end
            switch image_id_data.hand{im_idx}
                case 'L'
                    auto_stats.(test_dir).hand(i_im) = 1;
                case 'R'
                    auto_stats.(test_dir).hand(i_im) = 2;
            end
            auto_stats.(test_dir).visit(i_im) = image_id_data.visit(im_idx);
            auto_stats.(test_dir).digit(i_im) = image_id_data.digit(im_idx);
            auto_stats.(test_dir).people_id(i_im) = image_id_data.people_id(im_idx);

            %load([apex_gt_dir  im_name '_gt.mat']);
            load([apex_measures_dir im_name '_am.mat']);
            s = load([selected_apexes_dir im_name '_candidates.mat']);

            auto_stats.(test_dir).gradeable(i_im) = ~isempty(s.kept);
            auto_stats.(test_dir).num_distal_vessels(i_im) = sum(s.kept);
            auto_stats.(test_dir).num_nondistal_vessels(i_im) = sum(s.non_distal); 
            auto_stats.(test_dir).status(i_im) = s.status;

            if auto_stats.(test_dir).num_distal_vessels(i_im)
                auto_stats.(test_dir).mean_mean_width(i_im) = naNmean(apex_measures.mean_width);
                auto_stats.(test_dir).mean_weighted_width(i_im) = naNmean(apex_measures.mean_weighted_width);
                auto_stats.(test_dir).median_weighted_width(i_im) = naNmedian(apex_measures.mean_weighted_width);
                auto_stats.(test_dir).mean_median_width(i_im) = naNmean(apex_measures.median_width);
                auto_stats.(test_dir).mean_std_width(i_im) = naNmean(apex_measures.std_width);
                auto_stats.(test_dir).max_mean_width(i_im) = max(apex_measures.mean_weighted_width);          
                auto_stats.(test_dir).mean_max_width(i_im) = naNmean(apex_measures.max_width);
                auto_stats.(test_dir).max_max_width(i_im) = naNmax(apex_measures.max_width);
                auto_stats.(test_dir).mean_min_width(i_im) = naNmean(apex_measures.min_width);
                auto_stats.(test_dir).std_mean_width(i_im) = naNstd(apex_measures.mean_weighted_width);


                ori_entropy = mb_entropy(apex_measures.orientation_hist,2);
                auto_stats.(test_dir).mean_orientation_entropy(i_im) = naNmean(ori_entropy);
                auto_stats.(test_dir).median_orientation_entropy(i_im) = naNmedian(ori_entropy);
                auto_stats.(test_dir).std_orientation_entropy(i_im) = naNstd(ori_entropy);

                auto_stats.(test_dir).total_vessel_prob(i_im) = naNsum(apex_measures.total_prob);
                auto_stats.(test_dir).mean_vessel_prob(i_im) = naNmean(apex_measures.total_prob);

                auto_stats.(test_dir).total_scores(i_im) = naNsum(apex_measures.candidate_scores);
                auto_stats.(test_dir).mean_scores(i_im) = naNmean(apex_measures.candidate_scores);

                x1 = 1;
                if exist([vessel_centre_dir im_name '_vc.mat'], 'file');
                    load([vessel_centre_dir im_name '_vc.mat'], 'ncols');
                    x2 = ncols;
                else               
                    f = imfinfo([image_dir im_name '.png']);
                    x2 = f.Width;
                end

                d1 = x2 - x1;

                auto_stats.(test_dir).vessel_density(i_im) = auto_stats.(test_dir).num_distal_vessels(i_im) / d1;
                auto_stats.(test_dir).score_density(i_im) = auto_stats.(test_dir).total_scores(i_im) / d1;

            end
        end
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
if args.do_image_plots
    for feature_c = {...
            'num_distal_vessels',...
             'num_nondistal_vessels',...
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
                 'total_vessel_prob',...
                  'mean_vessel_prob',...
                      'total_scores',...
                       'mean_scores',...
                     'score_density',...
                    'vessel_density'};

            feature = feature_c{1};
            feature_str = feature;
            feature_str(feature_str == '_') = ' ';
            feature_str(1) = feature_str(1) - 32;

            dist_ss = [...
                auto_stats.final_test.(feature)(auto_stats.final_test.category == 1 & ~auto_stats.final_test.status);
                auto_stats.test_half.(feature)(auto_stats.test_half.category == 1 & ~auto_stats.test_half.status)];
            dist_ss(isnan(dist_ss)) = [];
            dist_pr = [...
                auto_stats.final_test.(feature)(auto_stats.final_test.category == 2 & ~auto_stats.final_test.status);
                auto_stats.test_half.(feature)(auto_stats.test_half.category == 2 & ~auto_stats.test_half.status)];
            dist_pr(isnan(dist_pr)) = [];
            dist_hc = [...
                auto_stats.final_test.(feature)(auto_stats.final_test.category == 3 & ~auto_stats.final_test.status);
                auto_stats.test_half.(feature)(auto_stats.test_half.category == 3 & ~auto_stats.test_half.status)];
            dist_hc(isnan(dist_hc)) = [];

            p_ss_hc = ranksum(dist_ss, dist_hc);
            p_ss_pr = ranksum(dist_ss, dist_pr);
            p_hc_pr = ranksum(dist_hc, dist_pr);

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
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
if args.do_people_plots
    people_ids = unique([auto_stats.final_test.people_id; auto_stats.test_half.people_id]);
    num_people = max(people_ids);

    people_stats.category = zeros(num_people,1);

    people_stats.present = false(num_people,5,2,2);

    people_stats.gradeable = false(num_people,5,2,2);
    people_stats.status = zeros(num_people,5,2,2);

    people_stats.num_distal_vessels = zeros(num_people,5,2,2);
    people_stats.num_nondistal_vessels = zeros(num_people,5,2,2);

    people_stats.mean_mean_width = nan(num_people,5,2,2);
    people_stats.mean_weighted_width = nan(num_people,5,2,2);
    people_stats.mean_median_width = nan(num_people,5,2,2);
    people_stats.mean_std_width = nan(num_people,5,2,2);
    people_stats.max_mean_width = nan(num_people,5,2,2);
    people_stats.mean_max_width = nan(num_people,5,2,2);
    people_stats.max_max_width = nan(num_people,5,2,2);
    people_stats.mean_min_width = nan(num_people,5,2,2);
    people_stats.std_mean_width = nan(num_people,5,2,2);   

    people_stats.total_vessel_prob = zeros(num_people,5,2,2);
    people_stats.mean_vessel_prob = zeros(num_people,5,2,2);

    people_stats.total_scores = zeros(num_people,5,2,2);
    people_stats.mean_scores = zeros(num_people,5,2,2);

    people_stats.mean_orientation_entropy = nan(num_people,5,2,2);
    people_stats.median_orientation_entropy = nan(num_people,5,2,2);
    people_stats.std_orientation_entropy = nan(num_people,5,2,2);

    people_stats.vessel_density = nan(num_people,5,2,2);
    people_stats.score_density = nan(num_people,5,2,2);

    for test_dirc = {'test_half', 'final_test'}%

        test_dir = test_dirc{1};
        for i_im = 1:length(auto_stats.(test_dir).category)

            p = auto_stats.(test_dir).people_id(i_im);

            if ~p
                continue;
            end

            v = auto_stats.(test_dir).visit(i_im);
            d = auto_stats.(test_dir).digit(i_im);
            h = auto_stats.(test_dir).hand(i_im);

            people_stats.present(p, d, h, v) = 1;
            people_stats.category(p) = auto_stats.(test_dir).category(i_im);

            for feature = {...
                             'gradeable',...
                                'status',...
                    'num_distal_vessels',...
                 'num_nondistal_vessels',...
                       'mean_mean_width',...
                   'mean_weighted_width',...
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
                     'total_vessel_prob',...
                      'mean_vessel_prob',...
                          'total_scores',...
                           'mean_scores',...
                         'score_density',...
                        'vessel_density'};

                people_stats.(feature{1})(p, d, h, v) = auto_stats.(test_dir).(feature{1})(i_im);
            end

        end
    end
    
    ss_idx = people_stats.category == 1;
    pr_idx = people_stats.category == 2;
    hc_idx = people_stats.category == 3;

    for feature_c = {...
            'num_distal_vessels',...
             'num_nondistal_vessels',...
                   'mean_mean_width',...
               'mean_weighted_width',...
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
                 'total_vessel_prob',...
                  'mean_vessel_prob',...
                      'total_scores',...
                       'mean_scores',...
                     'score_density',...
                    'vessel_density'};

            feature = feature_c{1};
            feature_str = feature;
            feature_str(feature_str == '_') = ' ';
            feature_str(1) = feature_str(1) - 32;

            display('***********************');
            display(['Feature: ' feature]);

            dist_f = people_stats.(feature);
            dist_f = dist_f(:,:);
            valid_ims = sum(people_stats.present(:,:) & ~people_stats.status(:,:) & ~isnan(dist_f), 2);

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
