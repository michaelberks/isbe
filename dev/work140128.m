test_dir = 'test_half';
apex_measures_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_metrics\'];
selected_apexes_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_maps\set12g_half_296655\mixed_maxima\'];
image_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\images\'];
f_mask_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\fov_masks\'];
    
im_list = dir([apex_measures_dir '*_am.mat']);


cat_names = {'SSc', 'PR', 'HC'};

features = {...
        'num_distal_vessels',...
         'median_weighted_width',...
    'median_orientation_entropy',...
                 'total_scores'};
prcts = [10 50 90];

for i_f = [2]
            
    feature = features{i_f};
    feature_str = feature;
    feature_str(feature_str == '_') = ' ';
    feature_str(1) = feature_str(1) - 32;
        
    for i_cat = 1:3
             
        valid_ims = auto_stats.(test_dir).category == i_cat & ~auto_stats.(test_dir).status & auto_stats.(test_dir).total_scores > 5;
        cat_list = im_list(valid_ims);
        dist_f = auto_stats.(test_dir).(feature)(valid_ims);
        
        pct_f = prctile(dist_f, prcts);
        
        for i_p = 1:3
            
            [~,im_idx] = min(abs(dist_f - pct_f(i_p)));
            
            im_name = cat_list(im_idx).name(1:6);
            
            im = u_load([image_dir im_name '.mat']);
            mask = u_load([f_mask_dir im_name '_f_mask.mat']);
            load([selected_apexes_dir im_name '_candidates.mat']);
            load([apex_measures_dir im_name '_am.mat']);
            
            figure;
            imgray(im);
            g_min = min(im(mask)) + 5;
            g_max = max(im(mask)) - 5;
            caxis([g_min g_max]);
            
            can_xy = candidate_xy(kept,:);
            
            for i_can = 1:sum(kept)

                text(can_xy(i_can,1), can_xy(i_can,2), num2str(apex_measures.candidate_scores(i_can),3), 'color', 'r');
                 text(can_xy(i_can,1), can_xy(i_can,2)+10, num2str(mb_entropy(apex_measures.orientation_hist(i_can,:),2),3), 'color', 'g');
                 text(can_xy(i_can,1), can_xy(i_can,2)+20, num2str(apex_measures.mean_weighted_width(i_can),3), 'color', 'b');
            end
            title([im_name ' from ' cat_names{i_cat} ' group, ' feature_str ' = ' num2str(dist_f(im_idx))  ' (' num2str(prcts(i_p)) '%-ile)']);
        end
    end
end

%%
test_dir = 'final_test';
hog_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_hogs\set12g_half_296655\mixed_maxima\']; 
mkdir(hog_dir);
im_list = dir(['C:\isbe\nailfold\data\rsa_study\' test_dir '\images\*.mat']);
for i_im = 1:length(im_list)
    
    im_name = im_list(i_im).name(1:6);
    if ~exist(['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_gt\' im_name '_gt.mat'], 'file')
        [candidates_hogs candidates_scores] =...
            compute_apex_candidate_hogs(...
            'selected_ims',         i_im,...
            'data_dir',             [nailfoldroot 'data/rsa_study/' test_dir '/'],...
            'feature_im_dir',       'images/',...
            'ori_dir',              'rf_regression/296621/',...
            'width_dir',            'rf_regression/297037/',...
            'candidates_dir',       'apex_maps\set12g_half_296655\mixed_maxima',...
            'feature_sigma',        0,...
            'ori_sigma',            0,...
            'width_sigma',          2,...
            'num_cells',            8,...
            'cell_sz',              8,... %Size of HoG cells in blocks
            'block_sz',             [2 2],...%Size of blocks in cells
            'num_ori_bins',         12,... %Number of bins in orientation histograms
            'norm_method',          'l1-sqrt',... %Method for local normalisation
            'block_spacing',        8,... %Separation between blocks in pixels - controls level of overlap
            'gradient_operator',    [-1 0 1],...
            'spatial_sigma',        0, ...
            'angle_wrap',           0,...
            'base_width',           20, ...
            'dist_thresh',          24^2);
        save([hog_dir im_name '.mat'], 'candidates*');
    end
end
%%
clear; pack;
im_list = dir('C:\isbe\nailfold\data\rsa_study\test_half\images\*.mat');
for i_set = 1:2    
    selected_ims = (i_set-1)*100 + (1:100);
    
    rf_class_args.sampling_args.y = [];
    rf_class_args.sampling_args.X = [];

    %candidates_hogs = [];
    for i_im = 1:100
        im_name = im_list(selected_ims(i_im)).name(1:6);
        s = load(['C:\isbe\nailfold\data\rsa_study\test_half\apex_hogs\set12g_half_296655\mixed_maxima\' im_name '.mat'],...
        'candidates_hogs', 'candidates_class', 'candidates_scores');
    
        rf_class_args.sampling_args.y = [rf_class_args.sampling_args.y; ...
            s.candidates_class]; clear candidates_class ;
        rf_class_args.sampling_args.X = [rf_class_args.sampling_args.X; ...
            s.candidates_hogs s.candidates_scores]; clear candidates_hogs candidates_scores;
    end

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
    model_dir = 'C:\isbe\nailfold\models\apex\final_pass\';
    rf_class_dir = [model_dir 'im_set_no_wrap' num2str(i_set) '\'];
    create_folder(rf_class_dir);
    save([rf_class_dir 'rf_args.mat'], 'rf_class_args');
    rf_class_args.tree_dir = [rf_class_dir 'trees/']; 

    apex_class_rf = random_forest_class_train(rf_class_args);

    save([rf_class_dir '\rf.mat'], 'apex_class_rf');
    clear rf_class_args;

    counts = zeros(1, 1:apex_class_rf.D);
    for i_tr = 1:length(apex_class_rf.trees)
        tree = u_load([apex_class_rf.tree_dir apex_class_rf.trees{i_tr}]);
        counts = counts + hist(tree.var(tree.var > 0), 1:apex_class_rf.D);
    end
    figure; bar(1:apex_class_rf.D, counts);
end
%%
im_list = dir('C:\isbe\nailfold\data\rsa_study\test_half\images\*.mat');
old_apexes_dir = 'C:\isbe\nailfold\data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima\';
new_apexes_dir = 'C:\isbe\nailfold\data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima_new2\';
mkdir(new_apexes_dir);    

for i_im = 1:601
    if i_im == 1
        load C:\isbe\nailfold\models\apex\final_pass\im_set_no_wrap2\rf.mat
    elseif i_im == 101
        load C:\isbe\nailfold\models\apex\final_pass\im_set_no_wrap1\rf.mat
    end
    
    im_name = im_list(i_im).name(1:6);
    load(['C:\isbe\nailfold\data\rsa_study\test_half\apex_hogs\set12g_half_296655\mixed_maxima\' im_name '.mat'],...
        'candidates_hogs', 'candidates_scores');
    load([old_apexes_dir im_name '_candidates.mat'], 'candidate_xy');

    [~,votes] = random_forest_class_predict(apex_class_rf, [candidates_hogs candidates_scores]);
    candidate_scores = votes(:,2) / 100; clear votes;
    
    save([new_apexes_dir im_name '_candidates.mat'], 'candidate_xy', 'candidate_scores');

end
%%
test_dir = 'final_test';

im_list = dir(['C:\isbe\nailfold\data\rsa_study\' test_dir '\images\*.mat']);

old_apexes_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_maps\set12g_half_296655\mixed_maxima\'];
new_apexes_dir = ['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_maps\set12g_half_296655\mixed_maxima_new\'];
mkdir(new_apexes_dir);    

load C:\isbe\nailfold\models\apex\final_pass\im_set_no_wrap2\rf.mat

for i_im = 1:length(im_list)
    
    display(['Processing image ' num2str(i_im)]);
    im_name = im_list(i_im).name(1:6);
    
    if ~exist(['C:\isbe\nailfold\data\rsa_study\' test_dir '\apex_gt\' im_name '_gt.mat'], 'file')
        
        load([hog_dir im_name '.mat'],...
            'candidates_hogs', 'candidates_scores');
        load([old_apexes_dir im_name '_candidates.mat'], 'candidate_xy');

        [~,votes] = random_forest_class_predict(apex_class_rf, [candidates_hogs candidates_scores]);
        candidate_scores = votes(:,2) / 100; clear votes;

        save([new_apexes_dir im_name '_candidates.mat'], 'candidate_xy', 'candidate_scores');
    end
end
%%
im_list = dir('C:\isbe\nailfold\data\rsa_study\test_half\images\*.mat');
for i_set = 1:6
    
    selected_ims = (i_set-1)*100 + (1:100);
    
    apex_class_pred = [];
    candidates_class = false(0,1);
    candidates_scores = [];
    %candidates_displacements = [];
    
    for i_im = 1:100
        im_name = im_list(selected_ims(i_im)).name(1:6);
        s = load(['C:\isbe\nailfold\data\rsa_study\test_half\apex_hogs\set12g_half_296655\mixed_maxima\' im_name '.mat'],...
            'candidates_class', 'candidates_scores', 'candidates_displacements');
        t = load(['C:\isbe\nailfold\data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima_new\' im_name '_candidates.mat'],...
            'candidate_scores');
        apex_class_pred = [apex_class_pred; t.candidate_scores]; %#ok
        candidates_class = [candidates_class; s.candidates_class]; %#ok
        candidates_scores = [candidates_scores; s.candidates_scores]; %#ok
        %candidates_displacements = [candidates_displacements; s.candidates_displacements]; %#ok
        
    end

    [roc_pts auc_new] = calculate_roc_curve(apex_class_pred, candidates_class);

    figure; axis equal; hold all;
    plot(roc_pts(:,1), roc_pts(:,2), '-');
    %plot(roc_pts(:,1), roc_pts(:,2), 'rx');
    axis([0 1 0 1]);

    [roc_pts auc_old] = calculate_roc_curve(candidates_scores, candidates_class);
    plot(roc_pts(:,1), roc_pts(:,2), 'g-');
    %plot(roc_pts(:,1), roc_pts(:,2), 'yx');
    axis([0 1 0 1]);
    legend({['New A_z = ' num2str(auc_new,3)], ['Old A_z = ' num2str(auc_old,3)]}, 'location', 'southeast');
end
%%
select_vessels_from_candidates_set('start_i', 1, 'end_i', [],...
    'do_fill_gaps', 0, 'do_distal_sub', 0, 'do_save', 1, 'do_plot', 0,...
    'data_dir', 'C:\isbe\nailfold\data\rsa_study\test_half\',...
    'input_dir', 'apex_maps\set12g_half_296655\mixed_maxima_new2\',...
    'output_dir', 'apex_maps\set12g_half_296655\mixed_maxima_new2\selected_candidates2\',...
    'upper_ydist', -90,...
    'lower_ydist', 30,...
    'initial_thresh', 0.7,...
    'class_map', class_map,...
    'weak_vessel_thresh', 0.5,...
    'strong_vessel_thresh', 0.9,...
    'bad_image_thresh', 0.7,...
    'angle_discard_thresh', pi/2.5);
%%
make_detection_results_struc(...
    'results_name', 'detection_results', ...
    'candidates_dir',    [nailfoldroot 'data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima_new2\selected_candidates\'],... %mandatory arguments
    'apex_gt_dir', [nailfoldroot 'data/rsa_study/test_half/apex_gt/'],...
    'results_dir', [nailfoldroot 'data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima_new2\selected_candidates\results\'],...
    'prob_dir', [nailfoldroot 'data/rsa_study/test_half/predictions/detection/rf_classification/296655/'],...
    'selected_gt', 1:601,...
    'selected_candidates', 1:601,...
    'selected_prob', 1:601);
%%
im_by_im_counts_half = analyse_detection_results(...
    'results_name', 'detection_results',... rf_offset_half_20131210T123607
    'candidates_dir', [nailfoldroot 'data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima_new\selected_candidates3\'],... %mandatory arguments
    'apex_gt_dir', [nailfoldroot 'data/rsa_study/test_half/apex_gt/'],...
    'results_dir', [nailfoldroot 'data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima_new\selected_candidates\results/'],...
    'selected_images', [],...
    'use_only_gradeable', 1,...
    'min_num_markers', 1,...
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
    'intermediate_selection', 1);
%%
im_list = dir('C:\isbe\nailfold\data\rsa_study\test_half\images\*.mat');
apex_dir = ...
    [nailfoldroot 'data\rsa_study\test_half\apex_maps\set12g_half_296655\'];
apex_gt_dir = ...
    [nailfoldroot 'data\rsa_study\test_half\apex_gt\'];

version_dir = {'mixed_maxima\selected_candidates', 'mixed_maxima_new_merged\test_candidates'};

results = cell(2,1);
results{1} = load([apex_dir 'mixed_maxima\selected_candidates\results\detection_results.mat']);
results{2} = load([apex_dir 'mixed_maxima_new_merged\selected_candidates\results\detection_results.mat']);

can_class = zeros(0,2);
can_scores = zeros(0,2);
can_displacements = zeros(0,2);
curr_row = 0;
for i_im = 1:601   

    im_name = im_list(i_im).name(1:6);
    load([apex_gt_dir im_name '_gt.mat'], 'is_distal', 'is_non_distal');
        
    for i_ver = 1:2
        can = load([apex_dir  version_dir{i_ver} '\' im_name '_candidates.mat'],...
            'candidate_scores', 'candidate_displacements');
        
        num_cans = length(can.candidate_scores);
        
        can_class_i = ones(num_cans,1);
        tp = results{i_ver}.detections{i_im,3}>0;
        apex_idx = results{i_ver}.detections{i_im,3}(tp);
        can_class_i(tp) = can_class_i(tp) + is_distal(apex_idx) + 2*is_non_distal(apex_idx);                
        can_class(curr_row+(1:num_cans),i_ver) = can_class_i;

        can_scores(curr_row+(1:num_cans),i_ver) = can.candidate_scores;
        if isempty(can.candidate_displacements)
            can.candidate_displacements = nan(size(can.candidate_scores));
        end
        can_displacements(curr_row+(1:num_cans),i_ver) = can.candidate_displacements;
    end
    curr_row = curr_row + num_cans;
end

score_dist = cell(2,3);
disp_dist = cell(2,3);
for i_ver = 1:2    
    grid_pts = linspace(0,max(can_scores(:,i_ver)),100);    
    
    for i_cl = 1:3
        score_dist{i_ver,i_cl} = build_1d_kernel_distribution(can_scores(can_class(:,i_ver)==i_cl,i_ver), grid_pts, 0);
        disp_dist{i_ver,i_cl} = build_1d_kernel_distribution(...
            can_displacements(can_class(:,i_ver)==i_cl & ~isnan(can_displacements(:,i_ver)),i_ver), linspace(-300,200,100), 0);
    end
end
%%
fig_dir = 'M:\nailfold\weekly_presentations\figures\new_candidate_scoring_method\';
colors = 'rgb';
figure; 
subplot(2,2,1);
hold all;
for i_cl = 1:3
    plot(score_dist{1,i_cl}.x, score_dist{1,i_cl}.D_f, colors(i_cl), 'linewidth', 2);
end
legend({'FP', 'Distal', 'Non-distal'});
xlabel('X');
ylabel('p(X|C_i)');
title('Marginal class PDFs for sum of RF location votes');

subplot(2,2,2);
hold all;
for i_cl = 1:3
    plot(score_dist{2,i_cl}.x, score_dist{2,i_cl}.D_f, colors(i_cl), 'linewidth', 2);
end
xlabel('X');
ylabel('p(X|C_i)');
title('Marginal class PDFs for new apex classification RF scores');

subplot(2,2,3);
hold all;
for i_cl = 1:3
    plot(disp_dist{1,i_cl}.x, disp_dist{1,i_cl}.D_f, colors(i_cl), 'linewidth', 2);
end
set(gca, 'ylim', [0 0.08]);

xlabel('X');
ylabel('p(X|C_i)');
title('Marginal class PDFs for displacements to estimated distal line: old method');

subplot(2,2,4);
hold all;
for i_cl = 1:3
    plot(disp_dist{2,i_cl}.x, disp_dist{2,i_cl}.D_f, colors(i_cl), 'linewidth', 2);
end
set(gca, 'ylim', [0 0.08]);

xlabel('X');
ylabel('p(X|C_i)');
title('Marginal class PDFs for displacements to estimated distal line: new method');
saveas(gcf, [fig_dir 'marginal_pdfs.fig']);

%%
% disp_diffs = nan(601,1);
% for i_im = 1:601   
% 
%     im_name = im_list(i_im).name(1:6);
%     can = cell(2,1);
%     for i_ver = 1:2
%         can{i_ver} = load([apex_dir  version_dir{i_ver} '\selected_candidates\' im_name '_candidates.mat'],...
%             'candidate_xy', 'candidate_displacements');
%     end
%     
%     idx2 = find(...
%         can{1}.candidate_xy(1,1) == can{2}.candidate_xy(:,1) & ...
%         can{1}.candidate_xy(1,2) == can{2}.candidate_xy(:,2));
%     
%     if ~isempty(can{1}.candidate_displacements) && ~isempty(can{2}.candidate_displacements)
%         disp_diffs(i_im) = can{1}.candidate_displacements(1) - can{2}.candidate_displacements(idx2);
%     end
% end
    
    
%%
combined_dist = cell(2,3);
for i_ver = 1:2
           
    [xx yy] = meshgrid(...
        linspace(0, max(can_scores(:,i_ver)), 100),...
        linspace(-200, 200, 100));
    grid_xy = [xx(:) yy(:)]; clear xx yy;    
    
    figure;
    for i_cl = 1:3
        cl_idx = can_class(:,i_ver)==i_cl & ~isnan(can_displacements(:,i_ver));        
        combined_dist{i_ver,i_cl} = build_2d_kernel_distribution(...
            [can_scores(cl_idx,i_ver), can_displacements(cl_idx,i_ver)], grid_xy, []);
        
        subplot(1,3,i_cl); imgray(reshape(combined_dist{i_ver,i_cl}.D_f,100,100) );
    end    
end
%%
colors = 'rgb';
label = {'old', 'new'};
for i_ver = 1:2
    
    figure; hold on;
    for i_cl = 1:3
        mesh(reshape(combined_dist{i_ver,i_cl}.x,100,100),  reshape(combined_dist{i_ver,i_cl}.y,100,100),...
            reshape(combined_dist{i_ver,i_cl}.D_f,100,100), 'edgecolor', colors(i_cl));
    end
    xlabel('X (candidate scores)');
    ylabel('Y (candidate displacements)');
    zlabel('p(X&Y|C_i)');
    title(['Joint class pdfs over candidate displacement and score:' label{i_ver} ' method']);
    saveas(gcf, [fig_dir 'joint_pdfs_' label{i_ver} '.fig']);
end
%%
conditional_class_probs = cell(2,1);
conditional_total_prob = cell(2,1);
posterior_class_probs = cell(2,1);
posterior_class = cell(2,1);
max_posterior_probs = cell(2,1);
prior_class_probs = cell(2,1);

colors = 'rgb';
for i_ver = 1:2
    conditional_class_probs{i_ver} = zeros(100,100,3);
    conditional_total_prob{i_ver} = zeros(100,100) + eps;
    prior_class_probs{i_ver} = zeros(3,1);
    
    for i_cl = 1:3
        prior_class_probs{i_ver}(i_cl) =...
            sum(can_class(:,i_ver)==i_cl) / size(can_class,1);
        
        conditional_class_probs{i_ver}(:,:,i_cl) = ...
            reshape(combined_dist{i_ver,i_cl}.D_f,100,100) * prior_class_probs{i_ver}(i_cl);
                    
        conditional_total_prob{i_ver} = conditional_total_prob{i_ver} +...
            conditional_class_probs{i_ver}(:,:,i_cl);
    end
    
    posterior_class_probs{i_ver} = zeros(100,100,3);
    
    figure; hold on;
    for i_cl = 1:3
        posterior_class_probs{i_ver}(:,:,i_cl) = ...
            conditional_class_probs{i_ver}(:,:,i_cl) ./ conditional_total_prob{i_ver};
        
        mesh(reshape(combined_dist{i_ver,i_cl}.x,100,100),  reshape(combined_dist{i_ver,i_cl}.y,100,100),...
            posterior_class_probs{i_ver}(:,:,i_cl), 'edgecolor', colors(i_cl));
    end
        
    xlabel('X (candidate scores)');
    ylabel('Y (candidate displacements)');
    zlabel('p(Ci|X&Y)');
    title(['Posterior class probabilities over candidate displacement and score:' label{i_ver} ' method']);
    saveas(gcf, [fig_dir 'posterior_pdfs_' label{i_ver} '.fig']);
    
    [max_posterior_probs{i_ver} posterior_class{i_ver}] = ...
        max(posterior_class_probs{i_ver}, [], 3);   
    
    figure; mesh(reshape(combined_dist{i_ver,i_cl}.x,100,100),  reshape(combined_dist{i_ver,i_cl}.y,100,100), ...
        conditional_total_prob{i_ver});
    xlabel('X (candidate scores)');
    ylabel('Y (candidate displacements)');
    zlabel('p(X&Y)');
    title(['Joint PDFs over candidate displacement and score:' label{i_ver} ' method']);
    saveas(gcf, [fig_dir 'data_pdfs_' label{i_ver} '.fig']);
    
    figure; imgray(posterior_class{i_ver}); colormap(eye(3));
    xlabel('X (candidate scores)');
    ylabel('Y (candidate displacements)');
    title(['Class map of MAP:' label{i_ver} ' method']);
end
%%
clear class_map;
class_map.x = reshape(combined_dist{2,1}.x,100,100);
class_map.y = reshape(combined_dist{2,1}.y,100,100);
class_map.post_class = posterior_class{2} - 1;
class_map.post_probs = max_posterior_probs{2};
%%
probs = cell(2,2);
for i_ver = 1:2
    valid_idx = ~isnan(can_displacements(:,i_ver));
    
    probs{i_ver,1} = interp2(...
        reshape(combined_dist{i_ver,1}.x,100,100),...
        reshape(combined_dist{i_ver,1}.y,100,100),...
        reshape(combined_dist{i_ver,1}.D_f,100,100), can_scores(valid_idx,i_ver), can_displacements(valid_idx,i_ver));
    probs{i_ver,2} = interp2(...
        reshape(combined_dist{i_ver,2}.x,100,100),...
        reshape(combined_dist{i_ver,2}.y,100,100),...
        reshape(combined_dist{i_ver,2}.D_f,100,100), can_scores(valid_idx,i_ver), can_displacements(valid_idx,i_ver));
    
    kept = probs{i_ver,1} > 3*probs{i_ver,2};
    display(['# TPs = ' num2str(sum(kept & can_class(valid_idx,i_ver)))]);
    display(['# FPs = ' num2str(sum(kept & ~can_class(valid_idx,i_ver)))]);
end

%%
select_vessels_from_candidates_set('start_i', 12, 'end_i', 12,...
    'do_fill_gaps', 0, 'do_distal_sub', 0, 'do_save', 1, 'do_plot', 0,...
    'data_dir', 'C:\isbe\nailfold\data\rsa_study\test_half\',...
    'candidates_dir', 'apex_maps\set12g_half_296655\mixed_maxima\',...
    'upper_ydist', -70,...
    'lower_ydist', 45)

%%
old_apexes_dir = 'C:\isbe\nailfold\data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima\';
new_apexes_dir = 'C:\isbe\nailfold\data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima_new\';

candidates_class = false(0,1);
candidates_scores = [];
candidates_displacements = [];
    
for i_im = 1:601
    im_name = im_list(i_im).name(1:6);
    
    s = load([old_apexes_dir im_name '_candidates.mat'], 'candidate_scores', 'candidate_displacements');    
    t = load(['C:\isbe\nailfold\data\rsa_study\test_half\apex_hogs\set12g_half_296655\mixed_maxima\' im_name '.mat'],...
        'candidates_class');

    candidates_class = [candidates_class; t.candidates_class]; %#ok
    candidates_scores = [candidates_scores; s.candidate_scores]; %#ok
    if isempty(s.candidate_displacements)
        s.candidate_displacements = nan(size(s.candidate_scores));
    end
    candidates_displacements = [candidates_displacements; s.candidate_displacements]; %#ok
end

tp = candidates_class;
fp =~tp;
figure; 
%subplot(1,3,1);
hold on;
plot(candidates_displacements(tp), candidates_scores(tp), 'g.');
plot(candidates_displacements(fp), candidates_scores(fp), 'r.');
%set(gca, 'xlim', [-200 200]);

kdist_tp_d = build_1d_kernel_distribution(candidates_displacements(tp & ~isnan(candidates_displacements)), linspace(-200,300,100), 0);
kdist_fp_d = build_1d_kernel_distribution(candidates_displacements(fp & ~isnan(candidates_displacements)), linspace(-200,300,100), 0);

figure; hold all;
plot(kdist_fp_d.x, kdist_fp_d.D_f, 'linewidth', 2);
plot(kdist_tp_d.x, kdist_tp_d.D_f, 'linewidth', 2);
%%
im_list = dir('C:\isbe\nailfold\data\rsa_study\test_half\images\*.mat');

apexes_dir = {'C:\isbe\nailfold\data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima_new\selected_candidates3\';
	'C:\isbe\nailfold\data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima_new_merged\test_candidates2\'};
image_dir = 'C:\isbe\nailfold\data\rsa_study\test_half\images\';
f_mask_dir = 'C:\isbe\nailfold\data\rsa_study\test_half\fov_masks\';

ap_markers = 'xo';
ap_colors = 'rg';
for i_im = (20*3)+(1:20)
    
    im_name = im_list(i_im).name(1:6);
            
    im = u_load([image_dir im_name '.mat']);
    mask = u_load([f_mask_dir im_name '_f_mask.mat']);

    figure;
    imgray(im);
    g_min = min(im(mask)) + 5;
    g_max = max(im(mask)) - 5;
    caxis([g_min g_max]);
    
    for i_ap = 1:2
        load([apexes_dir{i_ap} im_name '_candidates.mat']);
    
        can_xy = candidate_xy(kept,:);
        reject_xy = candidate_xy(~kept,:);   
        can_scores = candidate_scores(kept);
        
        if ~isempty(candidate_displacements)
            can_disp = candidate_displacements(kept);
        else
            can_disp = inf(size(kept));
        end
    
        plot(can_xy(:,1), can_xy(:,2), ['y' ap_markers(i_ap)]);
        plot(reject_xy(:,1), reject_xy(:,2), 'c.', 'markersize', 2);
        
        for i_can = 1:sum(kept) 
            
            text(can_xy(i_can,1)+10, can_xy(i_can,2) + (i_ap-1)*10,...
                [num2str(can_scores(i_can),3) ', ' num2str(can_disp(i_can),3)], 'color', ap_colors(i_ap));
            
        end
    end
end
%%
im_list = dir('C:\isbe\nailfold\data\rsa_study\test_half\images\*.mat');

prob_dir = 'C:\isbe\nailfold\data\rsa_study\test_half\predictions\detection\rf_classification\296655\';
apex_map_dir = 'C:\isbe\nailfold\data\rsa_study\test_half\apex_maps\set12g_half_296655\';
create_folder([apex_map_dir 'mixed_maxima_new_merged']);

dist_thresh = 60;
connect_thresh = 0.5;
n_connect_pts = 20;

g = gaussian_filters_1d(1);
g = g / sum(g);

for i_im = 1:length(im_list);

    im_name = im_list(i_im).name(1:6);
    display(['Processing image :' num2str(i_im)]);

    load([apex_map_dir 'mixed_maxima_new\' im_name '_candidates'],...
        'candidate_xy', 'candidate_scores');

    vessel_prob = u_load([prob_dir im_name '_pred.mat']);
    vessel_prob = conv2(g', g, vessel_prob, 'same');
    [candidate_xy, candidate_scores] = post_merge_candidates(candidate_xy, candidate_scores, ...
        vessel_prob, dist_thresh, connect_thresh, n_connect_pts, 0);

    save([apex_map_dir 'mixed_maxima_new_merged\' im_name '_candidates'],...
        'candidate_xy', 'candidate_scores');
end
%%
select_vessels_from_candidates_set('start_i', 1, 'end_i', [],...
    'do_fill_gaps', 0, 'do_distal_sub', 0, 'do_save', 1, 'do_plot', 0,...
    'data_dir', 'C:\isbe\nailfold\data\rsa_study\test_half\',...
    'input_dir', 'apex_maps\set12g_half_296655\mixed_maxima_new_merged\',...
    'output_dir', 'apex_maps\set12g_half_296655\mixed_maxima_new_merged\selected_candidates\',...
    'upper_ydist', -90,...
    'lower_ydist', 30,...
    'initial_thresh', 0.7,...
    'class_map', class_map,...
    'weak_vessel_thresh', 0.5,...
    'strong_vessel_thresh', 0.9,...
    'bad_image_thresh', 0.7,...
    'angle_discard_thresh', pi/2.5);
%%
make_detection_results_struc(...
    'results_name', 'detection_results', ...
    'candidates_dir',    [nailfoldroot 'data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima_new_merged\selected_candidates\'],... %mandatory arguments
    'apex_gt_dir', [nailfoldroot 'data/rsa_study/test_half/apex_gt/'],...
    'results_dir', [nailfoldroot 'data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima_new_merged\selected_candidates\results\'],...
    'prob_dir', [nailfoldroot 'data/rsa_study/test_half/predictions/detection/rf_classification/296655/'],...
    'selected_gt', 1:601,...
    'selected_candidates', 1:601,...
    'selected_prob', 1:601);
%%
analyse_detection_results(...
    'results_name', 'detection_results',... rf_offset_half_20131210T123607
    'candidates_dir', [nailfoldroot 'data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima_new_merged\test_candidates2\'],... %mandatory arguments
    'apex_gt_dir', [nailfoldroot 'data/rsa_study/test_half/apex_gt/'],...
    'results_dir', [nailfoldroot 'data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima_new_merged\selected_candidates\results/'],...
    'selected_images', [],...
    'use_only_gradeable', 1,...
    'min_num_markers', 1,...
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
analyse_detection_results(...
    'results_name', 'detection_results',... rf_offset_half_20131210T123607
    'candidates_dir', [nailfoldroot 'data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima\selected_candidates\'],... %mandatory arguments
    'apex_gt_dir', [nailfoldroot 'data/rsa_study/test_half/apex_gt/'],...
    'results_dir', [nailfoldroot 'data\rsa_study\test_half\apex_maps\set12g_half_296655\mixed_maxima\selected_candidates\results/'],...
    'selected_images', [],...
    'use_only_gradeable', 1,...
    'min_num_markers', 1,...
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
extract_apex_measures_set(...
        'data_dir', 'C:\isbe\nailfold\data\rsa_study\test_half\',...
        'num_jobs', 601, 'task_id', 140,...
        'prob_dir',             'rf_classification/296655/',...
        'ori_dir',              'rf_regression/296621/',...
        'width_dir',            'rf_regression/297037/',...
        'candidates_dir',       'apex_maps\set12g_half_296655\mixed_maxima_new_merged\test_candidates2',...
        'metrics_dir',          'apex_metrics\mixed_maxima_new_merged',...
        'fov_mask_dir',         'fov_masks/',...
        'prob_sigma',           1,...
        'ori_sigma',            0,...
        'width_sigma',          1,...
        'plot', 0);
%%
[candidates_hogs  candidates_class candidates_scores] =...
        compute_apex_candidate_hogs(...
        'results_name',         'rf_offset_mixed20140129T142615',...
        'selected_ims',         26,...
        'data_dir',             [nailfoldroot 'data/rsa_study/test_half/'],...
        'feature_im_dir',       'images/',...
        'prob_dir',             'rf_classification/296655/',...
        'ori_dir',              'rf_regression/296621/',...
        'width_dir',            'rf_regression/297037/',...
        'candidates_dir',       'apex_maps\set12g_half_296655\mixed_maxima',...
        'results_dir',          'apex_gt',...
        'feature_sigma',        0,...
        'ori_sigma',            0,...
        'width_sigma',          2,...
        'num_cells',            8,...
        'cell_sz',              8,... %Size of HoG cells in blocks
        'block_sz',             [2 2],...%Size of blocks in cells
        'num_ori_bins',         9,... %Number of bins in orientation histograms
        'norm_method',          'l1-sqrt',... %Method for local normalisation
        'block_spacing',        8,... %Separation between blocks in pixels - controls level of overlap
        'gradient_operator',    [-1 0 1],...
        'spatial_sigma',        0, ...
        'angle_wrap',           1,...
        'base_width',           20, ...
        'dist_thresh',          24^2);
%%
load C:\isbe\nailfold\models\apex\final_pass\class_map.mat
select_vessels_from_candidates_set('start_i', 1, 'end_i', [],...
    'do_fill_gaps', 0, 'do_distal_sub', 0, 'do_save', 1, 'do_plot', 0,...
    'data_dir', 'C:\isbe\nailfold\data\rsa_study\test_half\',...
    'input_dir', 'apex_maps\set12g_half_296655\mixed_maxima_new_merged\',...
    'output_dir', 'apex_maps\set12g_half_296655\mixed_maxima_new_merged\test_candidates2\',...
    'upper_ydist', -70,...
    'lower_ydist', 45,...
    'initial_thresh', 0.5,...
    'class_map', class_map,...
    'weak_vessel_thresh', 0.5,...
    'strong_vessel_thresh', 0.9,...
    'bad_image_thresh', 0.7,...
    'angle_discard_thresh', pi/2.5);
%%
load C:\isbe\nailfold\models\apex\final_pass\class_map.mat


select_vessels_from_candidates_set('start_i', 1, 'end_i', [],...
    'do_fill_gaps', 0, 'do_distal_sub', 0, 'do_save', 1, 'do_plot', 0,...
    'data_dir', 'C:\isbe\nailfold\data\rsa_study\final_test\',...
    'input_dir', 'apex_maps\set12g_half_296655\mixed_maxima_new\',...
    'output_dir', 'apex_maps\set12g_half_296655\mixed_maxima_new\selected_candidates\',...
    'upper_ydist', -70,...
    'lower_ydist', 45,...
    'initial_thresh', 0.5,...
    'class_map', class_map,...
    'weak_vessel_thresh', 0.5,...
    'strong_vessel_thresh', 0.9,...
    'bad_image_thresh', 0.7,...
    'angle_discard_thresh', pi/2.5);
%%
apex_gt_dir = [nailfoldroot 'data/rsa_study/final_test/apex_gt/'];
candidate_dir = [nailfoldroot 'data\rsa_study\final_test\apex_maps\set12g_half_296655\mixed_maxima_new\selected_candidates\'];

apex_gt_list = dir([apex_gt_dir '*.mat']);
candidate_list = dir([candidate_dir '*.mat']);

apex_gt_cell = cell(length(apex_gt_list),1);
for ii = 1:length(apex_gt_list)
    apex_gt_cell{ii} = apex_gt_list(ii).name(1:6);
end
candidate_cell = cell(length(candidate_list),1);
for ii = 1:length(candidate_list)
    candidate_cell{ii} = candidate_list(ii).name(1:6);
end
selected_candidates = find(ismember(candidate_cell, apex_gt_cell));

mkdir([nailfoldroot 'data/rsa_study/final_test/results/']); 
make_detection_results_struc(...
    'results_name', 'detection_results', ...
    'candidates_dir',    [nailfoldroot 'data\rsa_study\final_test\apex_maps\set12g_half_296655\mixed_maxima_new\selected_candidates\'],... %mandatory arguments
    'apex_gt_dir', [nailfoldroot 'data/rsa_study/final_test/apex_gt/'],...
    'results_dir', [nailfoldroot 'data\rsa_study\final_test\apex_maps\set12g_half_296655\mixed_maxima_new\selected_candidates\results\'],...
    'prob_dir', [nailfoldroot 'data/rsa_study/final_test/predictions/detection/rf_classification/296655/'],...
    'selected_gt', [],...
    'selected_candidates', selected_candidates,...
    'selected_prob', selected_candidates);
%%
extract_apex_measures_set(...
        'data_dir', 'C:\isbe\nailfold\data\rsa_study\final_test\',...
        'num_jobs', 1, 'task_id', 1,...
        'prob_dir',             'rf_classification/296655/',...
        'ori_dir',              'rf_regression/296621/',...
        'width_dir',            'rf_regression/297037/',...
        'candidates_dir',       'apex_maps\set12g_half_296655\mixed_maxima_new\selected_candidates',...
        'metrics_dir',          'apex_metrics\mixed_maxima_new',...
        'fov_mask_dir',         'fov_masks/',...
        'prob_sigma',           1,...
        'ori_sigma',            0,...
        'width_sigma',          1,...
        'plot', 0);
%%
analyse_detection_results(...
    'results_name', 'detection_results',... rf_offset_half_20131210T123607
    'candidates_dir', [nailfoldroot 'data\rsa_study\final_test\apex_maps\set12g_half_296655\mixed_maxima_new\selected_candidates\'],... %mandatory arguments
    'apex_gt_dir', [nailfoldroot 'data/rsa_study/final_test/apex_gt/'],...
    'results_dir', [nailfoldroot 'data\rsa_study\final_test\apex_maps\set12g_half_296655\mixed_maxima_new\selected_candidates\results/'],...
    'selected_images', [],...
    'use_only_gradeable', 0,...
    'min_num_markers', 1,...
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
    