%people_stats.(feature{1})(p, d, h, v)
im_dir = 'C:\isbe\nailfold\data\2_year_study\images\';
for feature = {...
                'num_distal_vessels',...
             'num_nondistal_vessels',...
               'mean_weighted_width',...
                    'max_mean_width',...
        'median_orientation_entropy',...
                    'vessel_density'};
      
    feature_str = feature{1};
    feature_str(feature_str == '_') = ' ';
    feature_str(1) = feature_str(1) - 32;
            
    figure;
    a = zeros(3,1);
    for ii = 1:6; a(ii) = subplot(2,3,ii); hold all; end
    title(a(2), ['Feature: \bf' feature_str ', \rm']);
    
    for i_p = 1:length(people_stats.super_category)

        s_cat = people_stats.super_category(i_p);
        if ~s_cat; continue; end

        baseline_im = squeeze(people_stats.present(i_p,:,:,1) & ~people_stats.status(i_p,:,:,1));

        [digit hand] = find(baseline_im);

        if isempty(digit)
            continue;
        elseif length(digit) > 1
            digit = digit(1);
            hand = hand(1);
        end

        followup_visits = find(squeeze(people_stats.present(i_p,digit,hand,5:end) & ~people_stats.status(i_p,digit,hand,5:end)))+1;
        feature_vals = squeeze(people_stats.(feature{1})(i_p, digit, hand, [1; followup_visits+3]));

        plot(a(s_cat), [1; followup_visits], feature_vals, '-x');
        plot(a(s_cat+3), [1; followup_visits], feature_vals - feature_vals(1), '-x');
        
        cat = people_stats.category(i_p);
        
        switch cat
            
            case 1
                cat_str = 'c';
            case 2
                cat_str = 'd';
            case 3
                cat_str = 'l';
            case 6
                cat_str = 'n';
            case 4
                cat_str = 'p';
                
            otherwise
                continue;               
        end        
        people_str = [cat_str zerostr(rem(people_stats.people_ids(i_p),1000),2)];
        visit_str = '156';
        switch hand
            case 1
                hand_str = 'L';
            case 2
                 hand_str = 'R';
        end 
        digit_str = ['D' num2str(digit) 'X3LrgMosaic.mat'];
        
        if i_p <= 20
            figure;

            for i_v = 1:3
                if ismember(i_v, [1; followup_visits])
                    im_name = [people_str 'V' visit_str(i_v) hand_str digit_str];
                    nailfold = u_load([im_dir im_name]);
                    subplot(3,1,i_v); imgray(nailfold);
                end               
            end
        end
            
    end
    
    a_min1 = inf;
    a_max1 = -inf;
    a_min2 = inf;
    a_max2 = -inf;
    
    for ii = 1:3
        ylim1 = get(a(ii), 'ylim');
        ylim2 = get(a(ii+3), 'ylim');
        
        a_min1 = min(a_min1, ylim1(1));
        a_max1 = max(a_max1, ylim1(2));
        a_min2 = min(a_min2, ylim2(1));
        a_max2 = max(a_max2, ylim2(2));
    end
    for ii = 1:3
        set(a(ii), 'ylim', [a_min1 a_max1]);
        set(a(ii+3), 'ylim', [a_min2 a_max2]);
    end
end
%%
for i_p = 10

    s_cat = people_stats.super_category(i_p);
    if ~s_cat; continue; end

    baseline_im = squeeze(people_stats.present(i_p,:,:,1) & ~people_stats.status(i_p,:,:,1));

    [digit hand] = find(baseline_im);

    if isempty(digit)
        continue;
    elseif length(digit) > 1
        digit = digit(1);
        hand = hand(1);
    end

    followup_visits = find(squeeze(people_stats.present(i_p,digit,hand,5:end) & ~people_stats.status(i_p,digit,hand,5:end)))+1;       
    cat = people_stats.category(i_p);

    switch cat

        case 1
            cat_str = 'c';
        case 2
            cat_str = 'd';
        case 3
            cat_str = 'l';
        case 6
            cat_str = 'n';
        case 4
            cat_str = 'p';

        otherwise
            continue;               
    end        
    people_str = [cat_str zerostr(rem(people_stats.people_ids(i_p),1000),2)];
    visit_str = '156';
    switch hand
        case 1
            hand_str = 'L';
        case 2
             hand_str = 'R';
    end 
    digit_str = ['D' num2str(digit) 'X3LrgMosaic'];

    visit_names = cell(3,1);
    for i_v = 1:3
        if ismember(i_v, [1; followup_visits])
            visit_names{i_v} = [people_str 'V' visit_str(i_v) hand_str digit_str];
        end               
    end
    
    display_automated_markup(... % non-strict mode
        'image_names',         visit_names,...
        'data_dir',             [nailfoldroot 'data/2_year_study/'],...
        'image_dir',            'images',...
        'vessel_centre_dir',    'vessel_centres',...
        'metrics_dir',          'apex_maps\set12g_half_296655\miccai_maxima\apex_metrics',...
        'candidates_dir',       'apex_maps\set12g_half_296655\miccai_maxima\rescores',...
        'selected_dir',         'apex_maps\set12g_half_296655\miccai_maxima\selected_apexes',...
        'selected_features', [],...
        'do_xls',           1,...
        'do_make_stats',    0,...
        'do_image_plots',   0, ...
        'do_people_plots',  0,...
        'fig_dir',          [],...
        'um_per_pix',       1.25,...
        'xls_filename',     'auto_stats.xls',...
        'aam_thresh',       -2e4,...
        'plot_rejected',    1,...
        'plot_r', 3,...
        'plot_c', 1);    
                
end
%%
compute_apex_candidate_hogs( ... % the user's input
    'image_names',          {'d49V1LD4X3LrgMosaic'},...
    'data_dir',             [nailfoldroot 'data/2_year_study/'],...
    'feature_im_dir',       'images',...
    'ori_dir',              'rf_regression/296621/',...
    'width_dir',            'rf_regression/297037/',...
    'candidates_dir',       'apex_maps\set12g_half_296655\island_maxima',...
    'hog_dir',              'apex_maps\set12g_half_296655\island_maxima\hogs\',...
    'feature_sigma',        0,...
    'prob_sigma',           2,...
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
    'angle_wrap',           1,...
    'base_width',           20, ...
    'dist_thresh',          24^2,...
    'overwrite',            1);
%%
rf = u_load('C:\isbe\nailfold\models\apex\rescoring\miccai_all\rf.mat');
predict_apex_candidate_rescore(...
    'image_names',          {'d49V1LD4X3LrgMosaic'},...
    'apex_class_rf',        rf,...
    'data_dir',             [nailfoldroot 'data/2_year_study/'],...
    'candidates_dir',       'apex_maps\set12g_half_296655\island_maxima',...
    'rescore_dir',          'apex_maps\set12g_half_296655\island_maxima\rescores',...
    'hog_dir',              'apex_maps\set12g_half_296655\island_maxima\hogs');
    
%%
load('C:\isbe\nailfold\data\rsa_study\data_lists\image_id_data.mat');
load('C:\isbe\nailfold\data\rsa_study\data_lists\miccai_lists.mat', 'miccai_selection');
im_names = image_id_data.im_names(miccai_selection.validation);
%%
[co_occurrence_n] = miccai_results_cooc_fun(...
    'image_names',          im_names,...
    'data_dir',             [nailfoldroot 'data/rsa_study/master_set/'],...
    'prob_dir',             'rf_classification/296655',...
    'cluster_dir',          'apex_clusters_merged',...
    'candidates_dir',       'apex_maps\set12g_half_296655\standard_maxima\rescores',...(new)
    'selected_dir',         'apex_maps\set12g_half_296655\standard_maxima\selected_apexes',...(new)
    'metrics_dir',          'apex_maps\set12g_half_296655\standard_maxima\apex_metrics(rf)',...(new)
    'width_feature',        'width_at_apex',...median_width
    'width_fudge', - 1,...
    'do_distal',            1,...
    'do_nondistal',         1,...
    'plot', 0);

%%
counts_s_d = hist(detected_widths_s, 1:30);
counts_s_a = hist([detected_widths_s; missed_widths_s], 1:30);
figure; 
subplot(1,2,1); bar(1:30, counts_s_a);
subplot(1,2,2); bar(1:30, counts_s_d ./ counts_s_a);
%
counts_n_d = hist(detected_widths_n, 1:30);
counts_n_a = hist([detected_widths_n; missed_widths_n], 1:30);
figure; 
subplot(1,2,1); bar(1:30, counts_n_a);
subplot(1,2,2); bar(1:30, counts_n_d ./ counts_n_a);
%
counts_i_d = hist(detected_widths_i, 1:30);
counts_i_a = hist([detected_widths_i; missed_widths_i], 1:30);
figure; 
subplot(1,2,1); bar(1:30, counts_i_a);
subplot(1,2,2); bar(1:30, counts_i_d ./ counts_i_a);
%%
im_list = dir('C:\isbe\nailfold\data\rsa_study\set12g_half\images\*.mat');
im_names = cell(450,1);
for i_im = 1:450;
    im_names{i_im} = im_list(i_im+450).name(1:end-4);
end
%%

fg_mask_dir = 'C:\isbe\nailfold\data\rsa_study\set12g_half\vessel_masks\';
label_dir = 'C:\isbe\nailfold\data\rsa_study\set12g_half\width_maps\';

prediction_dir = 'C:\isbe\nailfold\data\rsa_study\set12g_half\predictions\width\rf_regression\297037\';
[width_errors predicted_widths gt_widths] =...
    compute_image_width_errors(im_names, prediction_dir, fg_mask_dir, ... % non-strict mode
    'label_type',	[], ...
    'gt_widths',      [],...
    'centre_idx',   [],...
    'label_dir',	label_dir, ...
    'fov_mask_dir',	[],...
    'do_log', 0);

% prediction_dir = 'C:\isbe\nailfold\data\rsa_study\set12g_half\predictions\width\rf_regression\422303\';
% [width_errors predicted_widths gt_widths] =...
%     compute_image_width_errors(im_names, prediction_dir, fg_mask_dir, ... % non-strict mode
%     'label_type',	[], ...
%     'gt_widths',      [],...
%     'centre_idx',   [],...
%     'label_dir',	label_dir, ...
%     'fov_mask_dir',	[],...
%     'do_log', 1);


figure; plot(predicted_widths, width_errors, 'x');
hold all;
[~,m,b] = regression(predicted_widths, width_errors, 'one');
x = [0 80];
y = m*x + b;
plot(x, y, 'y--');

calibrated_widths = predicted_widths - (predicted_widths*m + b);
calibrated_errors = calibrated_widths - gt_widths;

figure; plot(gt_widths, width_errors, 'x');
figure; plot(gt_widths, calibrated_errors, 'x');
%%
prediction_dir = 'C:\isbe\nailfold\data\rsa_study\set12g_half\predictions\width\rf_regression\422303\';
label_dir = 'C:\isbe\nailfold\data\rsa_study\set12g_half\log_width_maps\';
[log_width_errors predicted_log_widths gt_log_widths] =...
    compute_image_width_errors(im_names, prediction_dir, fg_mask_dir, ... % non-strict mode
    'label_type',	[], ...
    'gt_widths',      [],...
    'centre_idx',   [],...
    'label_dir',	label_dir, ...
    'fov_mask_dir',	[],...
    'do_log', 0);


figure; plot(predicted_log_widths, log_width_errors, 'x');
hold all;
[~,m,b] = regression(predicted_log_widths, log_width_errors, 'one');
x = [0 10];
y = m*x + b;
plot(x, y, 'y--');

calibrated_log_widths = predicted_log_widths - (predicted_log_widths*m + b);
calibrated_log_errors = calibrated_log_widths - gt_log_widths;

figure; plot(gt_log_widths, log_width_errors, 'x');
figure; plot(gt_log_widths, calibrated_log_errors, 'x');
%%
calibrated_exp_widths = exp(calibrated_log_widths);
calibrated_exp_errors = calibrated_exp_widths - gt_widths;

figure; plot(gt_widths, width_errors, 'x');
figure; plot(gt_widths, calibrated_exp_errors, 'x');
%%
width_dir = 'C:\isbe\nailfold\data\rsa_study\set12g_half\width_maps\';
log_dir = 'C:\isbe\nailfold\data\rsa_study\set12g_half\log_width_maps\';
create_folder(log_dir);

w_list = dir([width_dir '*.mat']);
for i_w = 1:length(w_list)
    width_map = u_load([width_dir w_list(i_w).name]);
    log_width_map = log(width_map);
    log_width_map(~width_map) = 0;
    save([log_dir w_list(i_w).name], 'log_width_map');
end
%%
figure; 
subplot(1,2,1); hold all;
plot(sort(width_errors), linspace(0,1,length(width_errors)));
plot(sort(calibrated_errors), linspace(0,1,length(calibrated_errors)));
plot(sort(calibrated_exp_errors), linspace(0,1,length(calibrated_exp_errors)));
subplot(1,2,2); hold all;
plot(sort(abs(width_errors)), linspace(0,1,length(width_errors)));
plot(sort(abs(calibrated_errors)), linspace(0,1,length(calibrated_errors)));
plot(sort(abs(calibrated_exp_errors)), linspace(0,1,length(calibrated_exp_errors)));
%%
DECOMP_TYPE="{'g2dia','h2dia'}" NUM_ANGLES=6 WIN_SIZE=3 SIGMA_RANGE="[2 5]" OUTPUT_TYPE="width" PREDICTION_TYPE="rf_regression" WIDTH_DIR="log_width_maps" NUM_TREES=10 IMAGE_TYPE="real" DATA_ROOT="scratch/nailfold/" IMAGE_ROOT="data/rsa_study/set12g_half" NUM_SAMPLES=50000 MODEL_ROOT="models/vessel" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DATA_ROOT="scratch/nailfold/models/" MODEL_ROOT="vessel/width/rf_regression" MODEL_PATH="422303" MAKE_SAMPLED_MAPS=1 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
DATA_ROOT="scratch/nailfold/" MODEL_ROOT="models/vessel/width/rf_regression" MODEL_PATH="422303" NUM_JOBS=900 IMAGE_ROOT="data/rsa_study/set12g_half" USE_SAMPLED_MAPS=1 OVERWRITE=0 qsub -l twoday -hold_jid 422772 -V -t 451-900 matlab_code/trunk/hydra/cuc/predict_image_set.sh
%%
e_idx = false(450,1);
n_idx = false(450,1);
g_idx = false(450,1);
for i_im = 1:450
    e_idx(i_im) = im_names{i_im}(1)=='e';
    n_idx(i_im) = im_names{i_im}(1)=='n';
    g_idx(i_im) = im_names{i_im}(1)=='g';
end

display(sum(e_idx));
display(sum(n_idx));
display(sum(g_idx));
display(sum(e_idx | g_idx | n_idx));

discard_half = rand(450,1) < 0.25;
keep_idx = (e_idx&discard_half) | (n_idx&discard_half) | g_idx;

%%
im_names2 = im_names(keep_idx);
[~,~, gt_widths2] =...
    compute_image_width_errors(im_names2, prediction_dir, fg_mask_dir, ... % non-strict mode
    'label_type',	[], ...
    'gt_widths',      [],...
    'centre_idx',   [],...
    'label_dir',	label_dir, ...
    'fov_mask_dir',	[],...
    'do_log', 0);
%%
figure; 

im_names_e = im_names(e_idx);
[~,~, gt_widths_e] =...
    compute_image_width_errors(im_names_e, prediction_dir, fg_mask_dir, ... % non-strict mode
    'label_type',	[], ...
    'gt_widths',      [],...
    'centre_idx',   [],...
    'label_dir',	label_dir, ...
    'fov_mask_dir',	[],...
    'do_log', 0);

subplot(2,3,1); hist(gt_widths_e, 1:120);
subplot(2,3,4); hist(log(gt_widths_e), linspace(0,5,50));

im_names_n = im_names(n_idx);
[~,~, gt_widths_n] =...
    compute_image_width_errors(im_names_n, prediction_dir, fg_mask_dir, ... % non-strict mode
    'label_type',	[], ...
    'gt_widths',      [],...
    'centre_idx',   [],...
    'label_dir',	label_dir, ...
    'fov_mask_dir',	[],...
    'do_log', 0);
subplot(2,3,2); hist(gt_widths_n, 1:120);
subplot(2,3,5); hist(log(gt_widths_n), linspace(0,5,50));

im_names_g = im_names(g_idx);
[~,~, gt_widths_g] =...
    compute_image_width_errors(im_names_g, prediction_dir, fg_mask_dir, ... % non-strict mode
    'label_type',	[], ...
    'gt_widths',      [],...
    'centre_idx',   [],...
    'label_dir',	label_dir, ...
    'fov_mask_dir',	[],...
    'do_log', 0);
subplot(2,3,3); hist(gt_widths_g, 1:120);
subplot(2,3,6); hist(log(gt_widths_g), linspace(0,5,50));