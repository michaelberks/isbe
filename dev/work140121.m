prob_bins = linspace(0, 1, 50);
width_bins = linspace(0, 100, 50);

test_dir = 'test_half';
im_list = dir(['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_gt\*.mat']);
num_images = length(im_list);


prob_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\predictions\detection\rf_classification\296655\'];
width_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\predictions\width\rf_regression\297037\'];
measures_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_measures\'];

prob_dist = zeros(num_images, 50);
prob_sums = zeros(num_images, 1);
prob_means = zeros(num_images, 1);
vessel_sums = zeros(num_images, 1);

width_dist = zeros(num_images, 50);

for i_im = 1:num_images
    width_im = u_load([width_dir im_list(i_im).name(1:6) '_pred.mat']);
    prob_im = u_load([prob_dir im_list(i_im).name(1:6) '_pred.mat']);
    mask = width_im > 0;
    
    prob_dist(i_im,:) = hist(prob_im(mask), prob_bins);
    %width_dist(i_im,:) = hist(width_im(mask), width_bins);
    width_dist(i_im,:) = compute_weighted_histogram(width_im(mask), prob_im(mask), width_bins);
    prob_sums(i_im) = sum(prob_im(mask));
    prob_means(i_im) = mean(prob_im(mask));
    
    load([apex_measures_dir im_list(i_im).name(1:6) '_am.mat'], 'apex_measures');
    vessel_sums(i_im) = sum(apex_measures.total_prob);
end

prob_distn = bsxfun(@rdivide, prob_dist, sum(prob_dist,2));
width_distn = bsxfun(@rdivide, width_dist, sum(width_dist,2));
%%
ss_idx = vessels_stats.(test_dir).category == 1;
pr_idx = vessels_stats.(test_dir).category == 2;
hc_idx = vessels_stats.(test_dir).category == 3;

figure; hold all;
plot(prob_bins(2:end), mean(prob_distn(:,2:end)));
plot([prob_bins(2:end); prob_bins(2:end)],...
    [mean(prob_distn(:,2:end)) + std(prob_distn(:,2:end)); mean(prob_distn(:,2:end)) - std(prob_distn(:,2:end))], 'r');

figure; hold all;
plot(prob_bins(2:end), mean(prob_distn(ss_idx,2:end)));
plot([prob_bins(2:end); prob_bins(2:end)],...
    [mean(prob_distn(ss_idx,2:end)) + std(prob_distn(ss_idx,2:end)); mean(prob_distn(ss_idx,2:end)) - std(prob_distn(ss_idx,2:end))], 'r');

figure; hold all;
plot(prob_bins(2:end), mean(prob_distn(pr_idx,2:end)));
plot([prob_bins(2:end); prob_bins(2:end)],...
    [mean(prob_distn(pr_idx,2:end)) + std(prob_distn(pr_idx,2:end)); mean(prob_distn(pr_idx,2:end)) - std(prob_distn(pr_idx,2:end))], 'r');

figure; hold all;
plot(prob_bins(2:end), mean(prob_distn(hc_idx,2:end)));
plot([prob_bins(2:end); prob_bins(2:end)],...
    [mean(prob_distn(hc_idx,2:end)) + std(prob_distn(hc_idx,2:end)); mean(prob_distn(hc_idx,2:end)) - std(prob_distn(hc_idx,2:end))], 'r');
%%
figure; hold all;
plot(width_bins(:), mean(width_distn(:,:)));
plot([width_bins; width_bins],...
    [mean(width_distn(:,:)) + std(width_distn(:,:)); mean(width_distn(:,:)) - std(width_distn(:,:))], 'r');

figure; hold all;
plot(width_bins(:), mean(width_distn(ss_idx,:)));
% plot([width_bins; width_bins],...
%     [mean(width_distn(ss_idx,:)) + std(width_distn(ss_idx,:)); mean(width_distn(ss_idx,:)) - std(width_distn(ss_idx,:))], 'r');

% figure; hold all;
plot(width_bins(:), mean(width_distn(pr_idx,:)));
% plot([width_bins; width_bins],...
%     [mean(width_distn(pr_idx,:)) + std(width_distn(pr_idx,:)); mean(width_distn(pr_idx,:)) - std(width_distn(pr_idx,:))], 'r');

% figure; hold all;
plot(width_bins(:), mean(width_distn(hc_idx,:)));
% plot([width_bins; width_bins],...
%     [mean(width_distn(hc_idx,:)) + std(width_distn(hc_idx,:)); mean(width_distn(hc_idx,:)) - std(width_distn(hc_idx,:))], 'r');

%%
figure; hold all; a1 = gca;
figure; hold all; a2 = gca;
for i_gr = 1:7

    gr_idx = vessels_stats.final_test.grade_idx == i_gr;
    plot(a1, prob_bins(:), mean(prob_distn(gr_idx,:)));
    plot(a2, width_bins(:), mean(width_distn(gr_idx,:)));
end
legend(image_grade_labels(1:7));
%%
prob_sum_bins = linspace(min(prob_sums), max(prob_sums), 50);
vessel_sum_bins = linspace(min(vessel_sums), max(vessel_sums), 50);

figure; hold all; a1 = gca;
figure; hold all; a2 = gca;
figure; hold all; a3 = gca;
for i_gr = 1:7
    
    gr_idx = vessels_stats.(test_dir).grade_idx == i_gr;
    kdist_sum = build_1d_kernel_distribution(prob_sums(gr_idx), prob_sum_bins, 0);
    kdist_mean = build_1d_kernel_distribution(prob_means(gr_idx), prob_bins, 0);
    kdist_vsum = build_1d_kernel_distribution(vessel_sums(gr_idx), vessel_sum_bins, 0);
    
    plot(a1, prob_sum_bins, kdist_sum.D_a);
    plot(a2, prob_bins, kdist_mean.D_a);
    plot(a3, vessel_sum_bins, kdist_vsum.D_a);
end
legend(a1, vessels_stats.(test_dir).image_grade_labels(1:7));
legend(a3, vessels_stats.(test_dir).image_grade_labels(1:7));
legend(a3, vessels_stats.(test_dir).image_grade_labels(1:7));
    
%%
prob_sum_bins = linspace(min(prob_sums), max(prob_sums), 50);
figure; hold all; a1 = gca;
figure; hold all; a2 = gca;

for i_gr = 1:3
    
    gr_idx = vessels_stats.final_test.category == i_gr & vessels_stats.final_test.gradeable;
    kdist_sum = build_1d_kernel_distribution(prob_sums(gr_idx), prob_sum_bins, 0);
    kdist_mean = build_1d_kernel_distribution(prob_means(gr_idx), prob_bins, 0);
    
    plot(a1, prob_sum_bins, kdist_sum.D_a);
    plot(a2, prob_bins, kdist_mean.D_a);
    
    gr_idx = vessels_stats.final_test.category == i_gr & ~vessels_stats.final_test.gradeable;
    kdist_sum = build_1d_kernel_distribution(prob_sums(gr_idx), prob_sum_bins, 0);
    kdist_mean = build_1d_kernel_distribution(prob_means(gr_idx), prob_bins, 0);
    
    plot(a1, prob_sum_bins, kdist_sum.D_a, '--');
    plot(a2, prob_bins, kdist_mean.D_a, '--');
end
%%
figure; hold all;
u_idx = vessels_stats.test_half.grade_idx == 2;
plot(auto_stats.test_half.mean_vessel_prob(u_idx), auto_stats.test_half.num_distal_vessels(u_idx), 'rx');
plot(auto_stats.test_half.mean_vessel_prob(~u_idx), auto_stats.test_half.num_distal_vessels(~u_idx), 'g+');
%%
naughty = find(vessels_stats.final_test.gradeable & prob_sums < 5e4 & ~auto_stats.test_half.status);
for ii = 1:20
    im = u_load(['C:\isbe\nailfold\data\rsa_study\' test_dir '\images\' im_list(naughty(ii)).name(1:6) '.mat']);
    mask = u_load(['C:\isbe\nailfold\data\rsa_study\' test_dir '\fov_masks\' im_list(naughty(ii)).name(1:6) '_f_mask.mat']);
    
    figure;
    imgray(im);
    g_min = min(im(mask)) + 5;
    g_max = max(im(mask)) - 5;
    caxis([g_min g_max]);
end
%%
naughty = find(vessels_stats.test_half.gradeable & auto_stats.test_half.total_vessel_prob < 1e4 & ~auto_stats.test_half.status);
for ii = 1:23
    im = u_load(['C:\isbe\nailfold\data\rsa_study\' test_dir '\images\' im_list(naughty(ii)).name(1:6) '.mat']);
    mask = u_load(['C:\isbe\nailfold\data\rsa_study\' test_dir '\fov_masks\' im_list(naughty(ii)).name(1:6) '_f_mask.mat']);
    
    figure;
    imgray(im);
    g_min = min(im(mask)) + 5;
    g_max = max(im(mask)) - 5;
    caxis([g_min g_max]);
    load(['C:\isbe\nailfold\data\rsa_study\test\apex_clusters\' im_list(naughty(ii)).name(1:6) '_apex_clusters.mat']);
    plot_apex_clusters(vessels, gcf, 2)
    xlabel(num2str(auto_stats.(test_dir).mean_weighted_width(naughty(ii))));
end
%%
naughty = find(vessels_stats.(test_dir).gradeable & auto_stats.test_half.num_distal_vessels < 5 & ~auto_stats.test_half.status);
for ii = 1:20
    im = u_load(['C:\isbe\nailfold\data\rsa_study\' test_dir '\images\' im_list(naughty(ii)).name(1:6) '.mat']);
    mask = u_load(['C:\isbe\nailfold\data\rsa_study\' test_dir '\fov_masks\' im_list(naughty(ii)).name(1:6) '_f_mask.mat']);
    
    figure;
    imgray(im);
    g_min = min(im(mask)) + 5;
    g_max = max(im(mask)) - 5;
    caxis([g_min g_max]);
    load(['C:\isbe\nailfold\data\rsa_study\test\apex_clusters\' im_list(naughty(ii)).name(1:6) '_apex_clusters.mat']);
    plot_apex_clusters(vessels, gcf, 2)
    xlabel(num2str(auto_stats.(test_dir).mean_weighted_width(naughty(ii))));
end
%%
naughty = find(vessels_stats.(test_dir).gradeable & auto_stats.test_half.vessel_density < 0.004 & ~auto_stats.test_half.status);
for ii = 1:20
    im = u_load(['C:\isbe\nailfold\data\rsa_study\' test_dir '\images\' im_list(naughty(ii)).name(1:6) '.mat']);
    mask = u_load(['C:\isbe\nailfold\data\rsa_study\' test_dir '\fov_masks\' im_list(naughty(ii)).name(1:6) '_f_mask.mat']);
    
    figure;
    imgray(im);
    g_min = min(im(mask)) + 5;
    g_max = max(im(mask)) - 5;
    caxis([g_min g_max]);
    load(['C:\isbe\nailfold\data\rsa_study\test\apex_clusters\' im_list(naughty(ii)).name(1:6) '_apex_clusters.mat']);
    plot_apex_clusters(vessels, gcf, 2)
    xlabel(num2str(auto_stats.(test_dir).mean_weighted_width(naughty(ii))));
end
%%
naughty = find(vessels_stats.(test_dir).gradeable & auto_stats.test_half.score_density < 0.0035 & ~auto_stats.test_half.status);
for ii = 1:26
    im = u_load(['C:\isbe\nailfold\data\rsa_study\' test_dir '\images\' im_list(naughty(ii)).name(1:6) '.mat']);
    mask = u_load(['C:\isbe\nailfold\data\rsa_study\' test_dir '\fov_masks\' im_list(naughty(ii)).name(1:6) '_f_mask.mat']);
    
    figure;
    imgray(im);
    g_min = min(im(mask)) + 5;
    g_max = max(im(mask)) - 5;
    caxis([g_min g_max]);
    load(['C:\isbe\nailfold\data\rsa_study\test\apex_clusters\' im_list(naughty(ii)).name(1:6) '_apex_clusters.mat']);
    plot_apex_clusters(vessels, gcf, 2)
    xlabel([vessels_stats.(test_dir).grade{naughty(ii)} ': ' num2str(auto_stats.(test_dir).num_distal_vessels(naughty(ii)))]);
end
%%
naughty = find(~vessels_stats.(test_dir).gradeable & auto_stats.test_half.vessel_density > 0.004 & ~auto_stats.test_half.status);
for ii = 1:20
    im = u_load(['C:\isbe\nailfold\data\rsa_study\' test_dir '\images\' im_list(naughty(ii)).name(1:6) '.mat']);
    mask = u_load(['C:\isbe\nailfold\data\rsa_study\' test_dir '\fov_masks\' im_list(naughty(ii)).name(1:6) '_f_mask.mat']);
    
    figure;
    imgray(im);
    g_min = min(im(mask)) + 5;
    g_max = max(im(mask)) - 5;
    caxis([g_min g_max]);
    load(['C:\isbe\nailfold\data\rsa_study\test\apex_clusters\' im_list(naughty(ii)).name(1:6) '_apex_clusters.mat']);
    plot_apex_clusters(vessels, gcf, 2)
    xlabel(num2str(auto_stats.(test_dir).mean_weighted_width(naughty(ii))));
    
    load(['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_maps\set12g_half_296655\mixed_maxima\'...
        im_list(naughty(ii)).name(1:6) '_candidates.mat']);
    plot(candidate_xy(kept,1), candidate_xy(kept,2), 'r*');
end
%%
apex_gt_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_gt\'];
apex_measures_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_metrics\'];
candidates_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_maps\set12g_half_296655\mixed_maxima\'];
results_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\results\'];

im_list = dir([apex_gt_dir '*.mat']);
num_images = length(im_list);

load([results_dir 'rf_offset_mixed20140114T133732.mat']);

tp_prob_sums = zeros(0,1);
fp_prob_sums = zeros(0,1);
tp_can_scores = zeros(0,1);
fp_can_scores = zeros(0,1);

for i_im = 1:num_images   
        
    im_name = im_list(i_im).name(1:6);
    
    %Load in the GT structure
    load([apex_gt_dir im_name '_gt.mat']);

    load([candidates_dir im_name '_candidates.mat'], 'candidate_xy', 'kept', 'candidate_scores');%, 'non_distal', 'intermediate_selections'    
    load([apex_measures_dir im_name '_am.mat'], 'apex_measures');
    
    candidate_scores = candidate_scores(kept);
    tp = detections{i_im,3}(kept) > 0;
    fp = ~tp;
    
    tp_prob_sums = [tp_prob_sums; apex_measures.total_prob(tp)]; %#ok
    fp_prob_sums = [fp_prob_sums; apex_measures.total_prob(fp)]; %#ok
    tp_can_scores = [tp_can_scores; candidate_scores(tp)]; %#ok
    fp_can_scores = [fp_can_scores; candidate_scores(fp)]; %#ok    
    
end

figure; hold on;
plot(tp_can_scores, tp_prob_sums, 'gx');
plot(fp_can_scores, fp_prob_sums, 'r+');
%%
prob_sum_bins = linspace(0, 5000, 50);
kdist_fp = build_1d_kernel_distribution(fp_prob_sums, prob_sum_bins, 0);
kdist_tp = build_1d_kernel_distribution(tp_prob_sums, prob_sum_bins, 0);

figure; hold all;
plot(prob_sum_bins, kdist_tp.D_a, 'g', 'linewidth', 2);
plot(prob_sum_bins, kdist_fp.D_a, 'r', 'linewidth', 2);
%%
total_pts = 0;
for i_im = 1:num_images   
        
    total_pts = total_pts + length(detections{i_im,3});
end
%%
feature_X = zeros(total_pts, 10);
feature_y = false(total_pts, 1);
curr_pt = 0;
for i_im = 1:num_images   
        
    im_name = im_list(i_im).name(1:6);
    
    %Load in the GT structure
    load([candidates_dir im_name '_candidates.mat'], 'candidate_scores', 'candidate_displacements', 'vessel_centre_sum');%, 'non_distal', 'intermediate_selections'    
    load([apex_measures_dir 'all/' im_name '_am.mat'], 'apex_measures');
    
    idx = curr_pt + (1:length(candidate_scores));
    curr_pt = idx(end);
    
    if isempty(candidate_displacements)
        candidate_displacements = nan(size(candidate_scores));
    end
    
    feature_X(idx,1) = apex_measures.mean_width;
    feature_X(idx,2) = apex_measures.median_width;
    feature_X(idx,3) = apex_measures.std_width;
    feature_X(idx,4) = apex_measures.mean_weighted_width;
    feature_X(idx,5) = apex_measures.mean_weighted_prob;
    feature_X(idx,6) = apex_measures.total_prob;
    feature_X(idx,7) = mb_entropy(apex_measures.orientation_hist,2);
    feature_X(idx,8) = angle(apex_measures.base_orientation)/2;
    feature_X(idx,9) = candidate_scores;
    feature_X(idx,10) = candidate_displacements;
    %feature_X(idx,11) = vessel_centre_sum;
    
    feature_y(idx,1) = detections{i_im,3} > 0;
end
%%
bad_image_pts = isnan(feature_X(:,10));
feature_X(bad_image_pts,:) = [];
feature_y(bad_image_pts,:) = [];

total_pts = size(feature_X,1);

for i_dim = 1:9
    valid_pts = ~isnan(feature_X(:,i_dim));
    valid_vals = feature_X(valid_pts,i_dim);
    r_idx = ceil(sum(valid_pts)*rand(sum(~valid_pts),1));
    feature_X(~valid_pts,i_dim) = valid_vals(r_idx,1);
end
%%
train_idx = 1:floor(total_pts/2);
test_idx = floor(total_pts/2)+1:total_pts;
model_dir = 'C:\isbe\nailfold\models\apex\final\';
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
rf_class_args.sampling_args.y = feature_y(train_idx,:);
rf_class_args.sampling_args.X = feature_X(train_idx,:);
apex_class_rf = random_forest_class_train(rf_class_args);

[~,votes] = random_forest_class_predict(apex_class_rf, feature_X(test_idx,:));
apex_class_pred = votes(:,2) / rf_class_args.n_trees; clear votes;   

[roc_pts auc] = calculate_roc_curve(apex_class_pred, feature_y(test_idx,:));

figure; axis equal; hold on;
plot(roc_pts(:,1), roc_pts(:,2), '-');
plot(roc_pts(:,1), roc_pts(:,2), 'rx');
axis([0 1 0 1]);

