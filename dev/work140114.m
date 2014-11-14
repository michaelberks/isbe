90% specificity threshold for all grades: 0.36676
95% specificity threshold for all grades: 0.47898
99% specificity threshold for all grades: 0.86206

90% specificity threshold for all grades: 0.36712
95% specificity threshold for all grades: 0.4989
99% specificity threshold for all grades: 0.89426
%%
im_list = dir('C:\isbe\nailfold\data\rsa_study\test_half\images\*.mat');

apex_map_dir = 'C:\isbe\nailfold\data\rsa_study\test_half\apex_maps\';
fov_mask_dir = 'C:\isbe\nailfold\data\rsa_study\test_half\fov_masks\';
vessel_centres_dir = 'C:\isbe\nailfold\data\rsa_study\test_half\vessel_centres\full_centres\';
create_folder([apex_map_dir 'mixed_maxima']);

exclsuion_zone = 20;
apex_class_thresh = 0.5;
base_width = 20;

for i_im = 1:length(im_list);
    im_name = im_list(i_im).name(1:6);
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);
    
    load([apex_map_dir 'set12g_half_296655\' im_name '_pred.mat']);
    load([apex_map_dir 'frog\full_centres\' im_name '_pred.mat'], 'apex_class_pred');
    load([vessel_centres_dir im_name '_vc.mat']);
    f_mask = u_load([fov_mask_dir im_name '_f_mask.mat']);
    
    [discard_pts] = discard_edge_preds(vessel_centre, f_mask);
    include_pts = ~discard_pts & (apex_class_pred > apex_class_thresh);
    
    [apex_offset_map] = ...
        transform_apex_offset_preds(apex_class_pred, apex_offset_x_pred, apex_offset_y_pred,...
            vessel_centre, nrows, ncols, base_width, include_pts);
        
    [candidate_xy candidate_scores] = ...
        local_image_maxima(apex_offset_map, exclsuion_zone, f_mask, 0);
    
    save([apex_map_dir 'set12g_half_296655\mixed_maxima\' im_name '_candidates'],...
        'candidate_xy', 'candidate_scores');
end
%%
select_vessels_from_candidates_set('start_i', 1, 'end_i', 601,...
    'do_fill_gaps', 0, 'do_distal_sub', 0, 'do_save', 1, 'do_plot', 0,...
    'data_dir', 'C:\isbe\nailfold\data\rsa_study\test_half\',...
    'input_dir', 'apex_maps\set12g_half_296655\mixed_maxima\',...
    'output_dir', 'apex_maps\set12g_half_296655\mixed_maxima\selected_candidates\',...
    'upper_ydist', -70,...
    'lower_ydist', 45)
%%
make_detection_results_struc(...
    'results_name', 'rf_offset_mixed', ...
    'candidates_dir',    'C:\isbe\nailfold\data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima\',... %mandatory arguments
    'apex_gt_dir', [nailfoldroot 'data/rsa_study/test_half/apex_gt/'],...
    'results_dir', [nailfoldroot 'data/rsa_study/test_half/results/'],...
    'prob_dir', [nailfoldroot 'data/rsa_study/test_half/predictions/detection/rf_classification/296655/'],...
    'selected_gt', 1:601,...
    'selected_candidates', 1:601,...
    'selected_prob', 1:601);
%%
im_by_im_counts_half = analyse_detection_results(...
    'results_name', 'rf_offset_mixed20140129T142615',... rf_offset_half_20131210T123607
    'candidates_dir',    'C:\isbe\nailfold\data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima\',... %mandatory arguments
    'apex_gt_dir', [nailfoldroot 'data/rsa_study/test_half/apex_gt/'],...
    'results_dir', [nailfoldroot 'data/rsa_study/test_half/results/'],...
    'selected_images', [],...
    'use_only_gradeable', 1,...
    'min_num_markers', 2,...
    'max_missing_markers', inf,...
    'overall_summary', 1,...
    'summary_by_grade', 1,...
    'summary_by_shape', 1,...
    'summary_by_size', 1,...
    'compute_rocs', 1,...
    'analysis_by_width', 0,...
    'fixed_counting', 0,...
    'pct_counting', 0,...
    'perfect_counting', 0,...
    'final_selection', 1,...
    'intermediate_selection', 1);
%%
display('************* Results for SSc images ********************');
analyse_detection_results(...
    'results_name', 'rf_offset_mixed20140114T133732',... rf_offset_half_20131210T123607
    'candidates_dir',    'C:\isbe\nailfold\data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima\',... %mandatory arguments
    'apex_gt_dir', [nailfoldroot 'data/rsa_study/test_half/apex_gt/'],...
    'results_dir', [nailfoldroot 'data/rsa_study/test_half/results/'],...
    'selected_images', find(image_id_data.category_idx.ss(image_id_data.dataset_idx.validation)),...
    'use_only_gradeable', 1,...
    'min_num_markers', 2,...
    'max_missing_markers', inf,...
    'overall_summary', 1,...
    'summary_by_grade', 0,...
    'summary_by_shape', 0,...
    'summary_by_size', 0,...
    'compute_rocs', 0,...
    'analysis_by_width', 0,...
    'fixed_counting', 0,...
    'pct_counting', 0,...
    'perfect_counting', 0,...
    'final_selection', 1,...
    'intermediate_selection', 0);
display('');
display('');
display('************* Results for Primary Raynauds images ********************');
analyse_detection_results(...
    'results_name', 'rf_offset_mixed20140114T133732',... rf_offset_half_20131210T123607
    'candidates_dir',    'C:\isbe\nailfold\data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima\',... %mandatory arguments
    'apex_gt_dir', [nailfoldroot 'data/rsa_study/test_half/apex_gt/'],...
    'results_dir', [nailfoldroot 'data/rsa_study/test_half/results/'],...
    'selected_images', find(image_id_data.category_idx.pr(image_id_data.dataset_idx.validation)),...
    'use_only_gradeable', 1,...
    'min_num_markers', 2,...
    'max_missing_markers', inf,...
    'overall_summary', 1,...
    'summary_by_grade', 0,...
    'summary_by_shape', 0,...
    'summary_by_size', 0,...
    'compute_rocs', 0,...
    'analysis_by_width', 0,...
    'fixed_counting', 0,...
    'pct_counting', 0,...
    'perfect_counting', 0,...
    'final_selection', 1,...
    'intermediate_selection', 0);
display('');
display('');
display('************* Results for healthy control images ********************');
analyse_detection_results(...
    'results_name', 'rf_offset_mixed20140114T133732',... rf_offset_half_20131210T123607
    'candidates_dir',    'C:\isbe\nailfold\data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima\',... %mandatory arguments
    'apex_gt_dir', [nailfoldroot 'data/rsa_study/test_half/apex_gt/'],...
    'results_dir', [nailfoldroot 'data/rsa_study/test_half/results/'],...
    'selected_images', find(image_id_data.category_idx.hc(image_id_data.dataset_idx.validation)),...
    'use_only_gradeable', 1,...
    'min_num_markers', 2,...
    'max_missing_markers', inf,...
    'overall_summary', 1,...
    'summary_by_grade', 0,...
    'summary_by_shape', 0,...
    'summary_by_size', 0,...
    'compute_rocs', 0,...
    'analysis_by_width', 0,...
    'fixed_counting', 0,...
    'pct_counting', 0,...
    'perfect_counting', 0,...
    'final_selection', 1,...
    'intermediate_selection', 0);
%%
shape_labels = {'Non-specific', 'Normal', 'Angiogenic','Meandering'};
size_labels = {'Normal', 'Enlarged', 'Giant', 'Irregular', 'Undefined'}; 

image_dir = 'Q:\nailfold\data\rsa_study\images\caps\';
load C:\isbe\nailfold\data\rsa_study\image_id_data.mat

for test_dirc = {'final_test'}%'test_half', 
    
    test_dir = test_dirc{1};
    
    apex_gt_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_gt\'];
    vessel_centre_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\vessel_centres\full_centres\'];

    im_list = dir([apex_gt_dir '*_gt.mat']);

    num_images = length(im_list);
    vessels_stats.(test_dir).category = zeros(num_images,1);
    vessels_stats.(test_dir).gradeable = false(num_images,1);
    vessels_stats.(test_dir).grade = cell(num_images, 1);

    vessels_stats.(test_dir).num_distal_vessels = zeros(num_images,1);
    vessels_stats.(test_dir).num_nondistal_vessels = zeros(num_images,1);
    vessels_stats.(test_dir).num_undefined_vessels = zeros(num_images,1);

    vessels_stats.(test_dir).mean_width = nan(num_images,1);
    vessels_stats.(test_dir).median_width = nan(num_images,1);
    vessels_stats.(test_dir).max_width = nan(num_images,1);
    vessels_stats.(test_dir).std_width = nan(num_images,1);

    vessels_stats.(test_dir).shape_counts = nan(num_images, 4);
    vessels_stats.(test_dir).size_counts = nan(num_images, 5);

    vessels_stats.(test_dir).vessel_density = nan(num_images,1);

    for i_im = 1:num_images
        im_name = im_list(i_im).name(1:6);
        im_idx = find(strcmp(image_id_data.im_names, im_name));
        switch image_id_data.category{im_idx}
            case 'S'
                vessels_stats.(test_dir).category(i_im) = 1;
            case 'P'
                vessels_stats.(test_dir).category(i_im) = 2;
            case 'HC'
                vessels_stats.(test_dir).category(i_im) = 3;
        end
        
        load([apex_gt_dir  im_name '_gt.mat']);

        vessels_stats.(test_dir).gradeable(i_im) = gradeable;
        vessels_stats.(test_dir).grade{i_im} = majority_grade;

        vessels_stats.(test_dir).num_distal_vessels(i_im) = sum(is_distal);
        vessels_stats.(test_dir).num_nondistal_vessels(i_im) = sum(is_non_distal);
        vessels_stats.(test_dir).num_undefined_vessels(i_im) = sum(is_undefined);

        if vessels_stats.(test_dir).num_distal_vessels(i_im)
            vessels_stats.(test_dir).mean_width(i_im) = mean(apex_widths(is_distal));
            vessels_stats.(test_dir).median_width(i_im) = median(apex_widths(is_distal));
            vessels_stats.(test_dir).max_width(i_im) = max(apex_widths(is_distal));
            vessels_stats.(test_dir).std_width(i_im) = std(apex_widths(is_distal));

            for i_sh = 1:length(shape_labels)
                vessels_stats.(test_dir).shape_counts(i_im, i_sh) = ...
                    sum(strcmpi(apex_shape, shape_labels{i_sh}) & is_distal);
            end
            for i_sz = 1:length(size_labels)
                vessels_stats.(test_dir).size_counts(i_im, i_sz) = ...
                    sum(strcmpi(apex_size, size_labels{i_sz}) & is_distal);
            end
            
            x1 = 1;
            if exist([vessel_centre_dir im_name '_vc.mat'], 'file');
                load([vessel_centre_dir im_name '_vc.mat'], 'ncols');
                x2 = ncols;
            else               
                f = imfinfo([image_dir im_name '.bmp']);
                x2 = f.Width;
            end
            
%             p = polyfit(apex_xy(is_distal,1), apex_xy(is_distal,2),1);
%             y1 = polyval(p, x1);
%             y2 = polyval(p, x2);
%             d1 = sqrt((x1-x2)^2 + (y1-y2)^2);

            d1 = x2 - x1;

            vessels_stats.(test_dir).vessel_density(i_im) = vessels_stats.(test_dir).num_distal_vessels(i_im) / d1;
            
        end
    end
    [vessels_stats.(test_dir).grade_idx vessels_stats.(test_dir).image_grade_labels] = grp2idx(vessels_stats.(test_dir).grade);
    vessels_stats.(test_dir).shape_pct = bsxfun(@rdivide, vessels_stats.test_half.shape_counts, vessels_stats.test_half.num_distal_vessels);
    vessels_stats.(test_dir).size_pct = bsxfun(@rdivide, vessels_stats.test_half.size_counts, vessels_stats.test_half.num_distal_vessels);
end
%%               
for feature_c = {...
        'num_distal_vessels',...
         'num_nondistal_vessels',...
         'num_undefined_vessels',...
                    'mean_width',...
                  'median_width',...
                     'max_width',...
                     'std_width',...
                  'shape_counts',...
                   'size_counts',...
                   'shape_pct',...
                   'size_pct',...
                'vessel_density',...
                'grade_idx'};
            
        feature = feature_c{1};
        
        %for test_dir = {'test_half', 'final_test'}%

            for i_dim = 1:size(vessels_stats.final_test.(feature), 2)
                dist_ss = [...
                    vessels_stats.final_test.(feature)(vessels_stats.final_test.category == 1, i_dim);
                    vessels_stats.test_half.(feature)(vessels_stats.test_half.category == 1, i_dim)];
                dist_ss(isnan(dist_ss)) = [];
                dist_pr = [...
                    vessels_stats.final_test.(feature)(vessels_stats.final_test.category == 2, i_dim);
                    vessels_stats.test_half.(feature)(vessels_stats.test_half.category == 2, i_dim)];
                dist_pr(isnan(dist_pr)) = [];
                dist_hc = [...
                    vessels_stats.final_test.(feature)(vessels_stats.final_test.category == 3, i_dim);
                    vessels_stats.test_half.(feature)(vessels_stats.test_half.category == 3, i_dim)];
                dist_hc(isnan(dist_hc)) = [];

                p_ss_hc = ranksum(dist_ss, dist_hc);
                p_ss_pr = ranksum(dist_ss, dist_pr);
                p_hc_pr = ranksum(dist_hc, dist_pr);

                display(['Feature: ' feature '(' num2str(i_dim) ')']);
                display(['SS vs HC: p = ' num2str(p_ss_hc)]);
                display(['SS vs PR: p = ' num2str(p_ss_pr)]);
                display(['PR vs HC: p = ' num2str(p_hc_pr)]);
                display('');

                if strcmp(feature, 'grade_idx')
                    bins = 1:7;
                    counts_ss = hist(dist_ss, bins);
                    counts_ss = counts_ss / sum(counts_ss);

                    counts_pr = hist(dist_pr, bins);
                    counts_pr = counts_pr / sum(counts_pr);

                    counts_hc = hist(dist_hc, bins);
                    counts_hc = counts_hc / sum(counts_hc);
                    
                    figure; bar(bins, [counts_ss; counts_pr; counts_hc]');
                    set(gca, 'xticklabel', image_grade_labels);
                else
                    min_val = min([dist_ss; dist_pr; dist_hc]);
                    max_val = max([dist_ss; dist_pr; dist_hc]);
                    grid_pts = linspace(min_val, max_val, 100);

                    kdist_ss = build_1d_kernel_distribution(dist_ss, grid_pts, 0);
                    kdist_pr = build_1d_kernel_distribution(dist_pr, grid_pts, 0);
                    kdist_hc = build_1d_kernel_distribution(dist_hc, grid_pts, 0);

                    figure; hold all;
                    title(['Feature: ' feature]);
                    plot(kdist_ss.x, kdist_ss.D_a, 'linewidth', 2);
                    plot(kdist_pr.x, kdist_pr.D_a, 'linewidth', 2);
                    plot(kdist_hc.x, kdist_hc.D_a, 'linewidth', 2);
                end
 
            end
        %end
end
%%
%--------------------------------------------------------------------------
test_dir = 'test_half';
ss_idx = vessels_stats.(test_dir).category == 1;
hc_idx = vessels_stats.(test_dir).category == 3;

X_train = [...
    vessels_stats.(test_dir).grade_idx(ss_idx,:) ...
    vessels_stats.(test_dir).num_distal_vessels(ss_idx,:)...
    vessels_stats.(test_dir).num_nondistal_vessels(ss_idx,:)...
    vessels_stats.(test_dir).num_undefined_vessels(ss_idx,:)...
    vessels_stats.(test_dir).mean_width(ss_idx,:)...
    vessels_stats.(test_dir).max_width(ss_idx,:)...
    vessels_stats.(test_dir).std_width(ss_idx,:)...
    vessels_stats.(test_dir).shape_counts(ss_idx,:)...
    vessels_stats.(test_dir).size_counts(ss_idx,:)...
    vessels_stats.(test_dir).vessel_density(ss_idx,:); 
    ...
    vessels_stats.(test_dir).grade_idx(hc_idx,:) ...
    vessels_stats.(test_dir).num_distal_vessels(hc_idx,:)...
    vessels_stats.(test_dir).num_nondistal_vessels(hc_idx,:)...
    vessels_stats.(test_dir).num_undefined_vessels(hc_idx,:)...
    vessels_stats.(test_dir).mean_width(hc_idx,:)...
    vessels_stats.(test_dir).max_width(hc_idx,:)...
    vessels_stats.(test_dir).std_width(hc_idx,:)...
    vessels_stats.(test_dir).shape_counts(hc_idx,:)...
    vessels_stats.(test_dir).size_counts(hc_idx,:)...
    vessels_stats.(test_dir).vessel_density(hc_idx,:) ];

y_train = [true(sum(ss_idx), 1); false(sum(hc_idx),1)];

model_dir = 'C:\isbe\nailfold\models\image_class\';
mkdir model_dir;
rf_class_dir = [model_dir datestr(now, 30)];

warning('off', 'ASYM:unexpectedArgument');
rf_class_args.prediction_type = 'rf_classification';
rf_class_args.n_trees = 100;
rf_class_args.d = [];
rf_class_args.w_prior = 0;
rf_class_args.impure_thresh = 1.0000e-004;
rf_class_args.split_min = 100;
rf_class_args.end_cut_min = 25;
rf_class_args.do_ubound = 0;
rf_class_args.quiet = 1;
rf_class_args.overwrite = 0;
rf_class_args.minimise_size = 0;
rf_class_args.split_criterion = 'gdi';
rf_class_args.var_criterion = 'mabs';

rf_class_args.sampling_args.sampling_method = 'sample_saved_training_data';
rf_class_args.decomposition_args = [];

%Save the RF classification args (before we add in the data)
create_folder(rf_class_dir);
save([rf_class_dir 'rf_args.mat'], 'rf_class_args');
rf_class_args.tree_dir = [rf_class_dir 'trees/']; 

%Add the data/labels and train the forest
rf_class_args.sampling_args.y = y_train;
rf_class_args.sampling_args.X = X_train;
apex_class_rf = random_forest_class_train(rf_class_args);

%
test_dir = 'final_test';
ss_idx = vessels_stats.(test_dir).category == 1;
hc_idx = vessels_stats.(test_dir).category == 3;

X_test = [...
    vessels_stats.(test_dir).grade_idx(ss_idx,:) ...
    vessels_stats.(test_dir).num_distal_vessels(ss_idx,:)...
    vessels_stats.(test_dir).num_nondistal_vessels(ss_idx,:)...
    vessels_stats.(test_dir).num_undefined_vessels(ss_idx,:)...
    vessels_stats.(test_dir).mean_width(ss_idx,:)...
    vessels_stats.(test_dir).max_width(ss_idx,:)...
    vessels_stats.(test_dir).std_width(ss_idx,:)...
    vessels_stats.(test_dir).shape_counts(ss_idx,:)...
    vessels_stats.(test_dir).size_counts(ss_idx,:)...
    vessels_stats.(test_dir).vessel_density(ss_idx,:); 
    ...
    vessels_stats.(test_dir).grade_idx(hc_idx,:) ...
    vessels_stats.(test_dir).num_distal_vessels(hc_idx,:)...
    vessels_stats.(test_dir).num_nondistal_vessels(hc_idx,:)...
    vessels_stats.(test_dir).num_undefined_vessels(hc_idx,:)...
    vessels_stats.(test_dir).mean_width(hc_idx,:)...
    vessels_stats.(test_dir).max_width(hc_idx,:)...
    vessels_stats.(test_dir).std_width(hc_idx,:)...
    vessels_stats.(test_dir).shape_counts(hc_idx,:)...
    vessels_stats.(test_dir).size_counts(hc_idx,:)...
    vessels_stats.(test_dir).vessel_density(hc_idx,:) ];

y_test = [true(sum(ss_idx), 1); false(sum(hc_idx),1)];

[~,votes] = random_forest_class_predict(apex_class_rf, X_test);
image_class_pred = votes(:,2) / rf_class_args.n_trees; clear votes;

valid = ~any(isnan(X_test),2);

[roc_pts auc] = calculate_roc_curve(image_class_pred(valid), y_test(valid));

figure; axis equal; hold on;
plot(roc_pts(:,1), roc_pts(:,2), '-');
plot(roc_pts(:,1), roc_pts(:,2), 'rx');
axis([0 1 0 1]);
%% ------------------------------------------------------------------------
%%
load C:\isbe\nailfold\data\rsa_study\image_id_data.mat
image_dir = 'C:\isbe\nailfold\data\rsa_study\images\';
version_dir = {'mixed_maxima_new_merged\test_candidates2\', 'mixed_maxima_new\selected_candidates\'};
data_dir = {'test_half', 'final_test'};
for i_test = 1:2
    
    test_dir = data_dir{i_test};
    
    %apex_gt_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_gt\'];
    vessel_centre_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\vessel_centres\full_centres\'];
    apex_measures_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_metrics\'];
    selected_apexes_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_maps\set12g_half_296655\' version_dir{i_test}];
    
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
            
%             p = polyfit(apex_xy(is_distal,1), apex_xy(is_distal,2),1);
%             y1 = polyval(p, x1);
%             y2 = polyval(p, x2);
%             d1 = sqrt((x1-x2)^2 + (y1-y2)^2);

            d1 = x2 - x1;

            auto_stats.(test_dir).vessel_density(i_im) = auto_stats.(test_dir).num_distal_vessels(i_im) / d1;
            auto_stats.(test_dir).score_density(i_im) = auto_stats.(test_dir).total_scores(i_im) / d1;
            
        end
    end
end
%%
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
        %plot(kdist_ss.x, cumsum(kdist_ss.D_a), 'linewidth', 2);
        %plot(kdist_pr.x, cumsum(kdist_pr.D_a), 'linewidth', 2);
        %plot(kdist_hc.x, cumsum(kdist_hc.D_a), 'linewidth', 2);
        
        plot(prctile(dist_ss, 0:100), (0:100)/100, 'linewidth', 2);
        plot(prctile(dist_pr, 0:100), (0:100)/100, 'linewidth', 2);
        plot(prctile(dist_hc, 0:100), (0:100)/100, 'linewidth', 2);
        set(gca, 'ylim', [0 1]);
        
end
%%
for feature_c = {...
        'num_distal_vessels',...
         'num_nondistal_vessels',...
                    'mean_width',...
           'mean_weighted_width',...
             'mean_median_width',...
                'mean_width_std',...
                     'max_width',...
                     'std_width',...
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
    figure; 
    a1 = subplot(1,2,1); hold all; 
    a2 = subplot(1,2,2); hold all; 
    title(['Feature: ' feature]);
    
    for i_gr = 1:7
        feature_dist = auto_stats.test_half.(feature)(vessels_stats.test_half.grade_idx == i_gr);
        feature_dist(isnan(feature_dist)) = [];
        
        min_val = min(feature_dist);
        max_val = max(feature_dist);
        grid_pts = linspace(min_val, max_val, 100);

        kdist_ss = build_1d_kernel_distribution(feature_dist, grid_pts, 0);
    
        plot(a1, kdist_ss.x, kdist_ss.D_a, 'linewidth', 2);
        plot(a2, kdist_ss.x, cumsum(kdist_ss.D_a), 'linewidth', 2);
        
    end
    legend(a1, vessels_stats.test_half.image_grade_labels(1:7));
    set(a2, 'ylim', [0 1]);   
end

%%
%%
%--------------------------------------------------------------------------
test_dir = 'test_half';
ss_idx = auto_stats.(test_dir).category == 1 & ~auto_stats.(test_dir).status;
hc_idx = auto_stats.(test_dir).category == 3 & ~auto_stats.(test_dir).status;
            
X_train = [...
    auto_stats.(test_dir).num_distal_vessels(ss_idx,:) ...
    auto_stats.(test_dir).num_nondistal_vessels(ss_idx,:)...
    auto_stats.(test_dir).mean_width(ss_idx,:)...
    auto_stats.(test_dir).mean_weighted_width(ss_idx,:)...
    auto_stats.(test_dir).mean_median_width(ss_idx,:)...
    auto_stats.(test_dir).mean_width_std(ss_idx,:)...
    auto_stats.(test_dir).max_width(ss_idx,:)...
    auto_stats.(test_dir).std_width(ss_idx,:)...
    auto_stats.(test_dir).mean_orientation_entropy(ss_idx,:)...
    auto_stats.(test_dir).median_orientation_entropy(ss_idx,:)...
    auto_stats.(test_dir).std_orientation_entropy(ss_idx,:)...
    auto_stats.(test_dir).vessel_density(ss_idx,:)... 
    auto_stats.(test_dir).total_vessel_prob(ss_idx,:)...
    auto_stats.(test_dir).mean_vessel_prob(ss_idx,:)...
    auto_stats.(test_dir).total_scores(ss_idx,:)...
    auto_stats.(test_dir).mean_scores(ss_idx,:)...
    auto_stats.(test_dir).score_density(ss_idx,:);
    ...
    auto_stats.(test_dir).num_distal_vessels(hc_idx,:) ...
    auto_stats.(test_dir).num_nondistal_vessels(hc_idx,:)...
    auto_stats.(test_dir).mean_width(hc_idx,:)...
    auto_stats.(test_dir).mean_weighted_width(hc_idx,:)...
    auto_stats.(test_dir).mean_median_width(hc_idx,:)...
    auto_stats.(test_dir).mean_width_std(hc_idx,:)...
    auto_stats.(test_dir).max_width(hc_idx,:)...
    auto_stats.(test_dir).std_width(hc_idx,:)...
    auto_stats.(test_dir).mean_orientation_entropy(hc_idx,:)...
    auto_stats.(test_dir).median_orientation_entropy(hc_idx,:)...
    auto_stats.(test_dir).std_orientation_entropy(hc_idx,:)...
    auto_stats.(test_dir).vessel_density(hc_idx,:)... 
    auto_stats.(test_dir).total_vessel_prob(hc_idx,:)...
    auto_stats.(test_dir).mean_vessel_prob(hc_idx,:)...
    auto_stats.(test_dir).total_scores(hc_idx,:)...
    auto_stats.(test_dir).mean_scores(hc_idx,:)...
    auto_stats.(test_dir).score_density(hc_idx,:) ];

for i_dim = 1:size(X_train,2);
    valid_pts = ~isnan(X_train(:,i_dim));
    valid_vals = X_train(valid_pts,i_dim);
    r_idx = ceil(sum(valid_pts)*rand(sum(~valid_pts),1));
    X_train(~valid_pts,i_dim) = valid_vals(r_idx,1);
end

y_train = [true(sum(ss_idx), 1); false(sum(hc_idx),1)];

model_dir = 'C:\isbe\nailfold\models\image_class\';
mkdir model_dir;
rf_class_dir = [model_dir datestr(now, 30)];

warning('off', 'ASYM:unexpectedArgument');
rf_class_args.prediction_type = 'rf_classification';
rf_class_args.n_trees = 100;
rf_class_args.d = [];
rf_class_args.w_prior = 0;
rf_class_args.impure_thresh = 1.0000e-004;
rf_class_args.split_min = 100;
rf_class_args.end_cut_min = 25;
rf_class_args.do_ubound = 0;
rf_class_args.quiet = 1;
rf_class_args.overwrite = 0;
rf_class_args.minimise_size = 0;
rf_class_args.split_criterion = 'gdi';
rf_class_args.var_criterion = 'mabs';

rf_class_args.sampling_args.sampling_method = 'sample_saved_training_data';
rf_class_args.decomposition_args = [];

%Save the RF classification args (before we add in the data)
create_folder(rf_class_dir);
save([rf_class_dir 'rf_args.mat'], 'rf_class_args');
rf_class_args.tree_dir = [rf_class_dir 'trees/']; 

%Add the data/labels and train the forest
rf_class_args.sampling_args.y = y_train;
rf_class_args.sampling_args.X = X_train;
apex_class_rf = random_forest_class_train(rf_class_args);

%
test_dir = 'final_test';
ss_idx = auto_stats.(test_dir).category == 1 & ~auto_stats.(test_dir).status;
hc_idx = auto_stats.(test_dir).category == 3 & ~auto_stats.(test_dir).status;

X_test = [...
    auto_stats.(test_dir).num_distal_vessels(ss_idx,:) ...
    auto_stats.(test_dir).num_nondistal_vessels(ss_idx,:)...
    auto_stats.(test_dir).mean_width(ss_idx,:)...
    auto_stats.(test_dir).mean_weighted_width(ss_idx,:)...
    auto_stats.(test_dir).mean_median_width(ss_idx,:)...
    auto_stats.(test_dir).mean_width_std(ss_idx,:)...
    auto_stats.(test_dir).max_width(ss_idx,:)...
    auto_stats.(test_dir).std_width(ss_idx,:)...
    auto_stats.(test_dir).mean_orientation_entropy(ss_idx,:)...
    auto_stats.(test_dir).median_orientation_entropy(ss_idx,:)...
    auto_stats.(test_dir).std_orientation_entropy(ss_idx,:)...
    auto_stats.(test_dir).vessel_density(ss_idx,:)... 
    auto_stats.(test_dir).total_vessel_prob(ss_idx,:)...
    auto_stats.(test_dir).mean_vessel_prob(ss_idx,:)...
    auto_stats.(test_dir).total_scores(ss_idx,:)...
    auto_stats.(test_dir).mean_scores(ss_idx,:)...
    auto_stats.(test_dir).score_density(ss_idx,:);
    ...
    auto_stats.(test_dir).num_distal_vessels(hc_idx,:) ...
    auto_stats.(test_dir).num_nondistal_vessels(hc_idx,:)...
    auto_stats.(test_dir).mean_width(hc_idx,:)...
    auto_stats.(test_dir).mean_weighted_width(hc_idx,:)...
    auto_stats.(test_dir).mean_median_width(hc_idx,:)...
    auto_stats.(test_dir).mean_width_std(hc_idx,:)...
    auto_stats.(test_dir).max_width(hc_idx,:)...
    auto_stats.(test_dir).std_width(hc_idx,:)...
    auto_stats.(test_dir).mean_orientation_entropy(hc_idx,:)...
    auto_stats.(test_dir).median_orientation_entropy(hc_idx,:)...
    auto_stats.(test_dir).std_orientation_entropy(hc_idx,:)...
    auto_stats.(test_dir).vessel_density(hc_idx,:)... 
    auto_stats.(test_dir).total_vessel_prob(hc_idx,:)...
    auto_stats.(test_dir).mean_vessel_prob(hc_idx,:)...
    auto_stats.(test_dir).total_scores(hc_idx,:)...
    auto_stats.(test_dir).mean_scores(hc_idx,:)...
    auto_stats.(test_dir).score_density(hc_idx,:) ];

for i_dim = 1:size(X_test,2);
    valid_pts = ~isnan(X_test(:,i_dim));
    valid_vals = X_test(valid_pts,i_dim);
    r_idx = ceil(sum(valid_pts)*rand(sum(~valid_pts),1));
    X_test(~valid_pts,i_dim) = valid_vals(r_idx,1);
end

y_test = [true(sum(ss_idx), 1); false(sum(hc_idx),1)];

[~,votes] = random_forest_class_predict(apex_class_rf, X_test);
image_class_pred = votes(:,2) / rf_class_args.n_trees; clear votes;
valid = ~any(isnan(X_test),2);

[roc_pts auc] = calculate_roc_curve(image_class_pred(valid), y_test(valid));

figure; axis equal; hold on;
plot(roc_pts(:,1), roc_pts(:,2), '-');
plot(roc_pts(:,1), roc_pts(:,2), 'rx');
axis([0 1 0 1]);
%%
test_dir = 'final_test';
ss_idx = auto_stats.(test_dir).category == 1 & ~auto_stats.(test_dir).status;
hc_idx = auto_stats.(test_dir).category == 3 & ~auto_stats.(test_dir).status;

figure; hold on;
plot(auto_stats.(test_dir).total_scores(hc_idx,:), auto_stats.(test_dir).median_orientation_entropy(hc_idx,:), 'gx');
plot(auto_stats.(test_dir).total_scores(ss_idx,:), auto_stats.(test_dir).median_orientation_entropy(ss_idx,:), 'r+');
%%
for i_im = find(auto_stats.final_test.median_orientation_entropy > 5 & hc_idx)'
    im_idx = find(strcmp(image_id_data.im_names, im_list(i_im).name(1:6)));
    view_subject(image_id_data.people_id(im_idx),'visit', image_id_data.visit(im_idx));
    set(gcf, 'name', im_list(i_im).name(1:6));
end
%%
i_ims = find(auto_stats.final_test.total_scores < 5 & hc_idx)';
for ii = 1:10
    i_im = i_ims(ii);
    im_idx = find(strcmp(image_id_data.im_names, im_list(i_im).name(1:6)));
    view_subject(image_id_data.people_id(im_idx),'visit', image_id_data.visit(im_idx));
    set(gcf, 'name', im_list(i_im).name(1:6));
end
%%
i_ims = find(auto_stats.final_test.total_scores > 5 &...
    auto_stats.final_test.score_density > 4e-3 &...
    ismember(vessels_stats.final_test.grade_idx, [1 7]))';
for ii = 1:10
    i_im = i_ims(ii);
    im_idx = find(strcmp(image_id_data.im_names, im_list(i_im).name(1:6)));
    view_subject(image_id_data.people_id(im_idx),'visit', image_id_data.visit(im_idx));
    set(gcf, 'name', im_list(i_im).name(1:6));
end
%%
figure; hold all;
plot(auto_stats.final_test.mean_orientation_entropy(hc_idx), auto_stats.final_test.mean_weighted_width(hc_idx), 'rx');
plot(auto_stats.final_test.mean_orientation_entropy(ss_idx), auto_stats.final_test.mean_weighted_width(ss_idx), 'go');
%%
figure; hold all;
plot(auto_stats.final_test.mean_orientation_entropy(hc_idx), auto_stats.final_test.num_distal_vessels(hc_idx), 'rx');
plot(auto_stats.final_test.mean_orientation_entropy(ss_idx), auto_stats.final_test.num_distal_vessels(ss_idx), 'go');
%%
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
%%
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