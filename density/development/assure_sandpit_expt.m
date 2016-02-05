image_root = 'c:\isbe\density\assure\375_125\';
model_root = 'c:\isbe\density\assure\models\';
image_dir = 'c:\isbe\density\assure\375_125\images\';
full_mask_dir = 'c:\isbe\density\assure\375_125\full_masks\';
compressed_mask_dir = 'c:\isbe\density\assure\375_125\compressed_masks\';
label_mask_dir = 'c:\isbe\density\assure\375_125\label_masks\';

for i_fold = 1:5
    shift_images = 500 - i_fold*100;
    build_predictor(...
        'predictor_name',       'rf', ...
        'task_id',				1, ...
        'job_id',				['test_forest_' zerostr(i_fold,1)], ...
        ... % Folders
        'image_root',           image_root,...
        'model_root',           model_root, ...
        ... % Output parameters
        'output_type',          'class_label', ...
        ... % Sampling parameters
        'num_samples',			1000, ...
        'max_n_images',         400, ...
        'shift_images',         shift_images, ...
        'bg_ratio',				0, ...
        'replace_sample',       false, ...
        'sampling_method',      'generate_training_data', ...
        'image_type', 			'real', ...
            ... % Image sampling parameters
            'image_dir',         'images',...
            'fg_mask_dir',       'compressed_masks',...
            'fov_mask_dir',      [],...
            'class_label_dir',   'label_masks',...
            'make_resampling_maps',  0, ...
        ... % Image feature/decomposition parameters
        'num_levels', 			6, ...
        'rgb_channel',          'rgb', ...
        'normalise', 			0, ...
        'win_size',				1, ...
        'decomp_type', 			'dt', ...
            ... % DTCWT parameters
            'feature_shape', 		'rect', ...
            'feature_type',			'conj', ...
        ... % Predictor parameters
        'prediction_type',		'rf_classification', ...
            ... % Tree/Forest parameters
            'n_trees',				2, ...
            'split_min',			10, ...
        ... % Miscellaneous parameters
        'overwrite',			true,...
        'quiet',                1);
end
%%
predict_image_set(...
    'model_id',         'class_label/rf_classification/test_forest_1',...
    'image_dir',        'images',...
    'prediction_dir',   'predictions',...
    'model_name',       'rf01',...
    'model_root',		model_root, ...
    'image_root',       image_root,...
    'task_id',			1, ...
    'num_jobs',			5, ...
    'mask_dir',			'compressed_masks', ...
    'max_size',			1024, ...
    'use_sampled_maps', 0,...
    'use_probs',		1, ...
    'overwrite',		false);
%%
DECOMP_TYPE="{'dt'}" NUM_LEVELS=6 WIN_SIZE=1 OUTPUT_TYPE="class_label" PREDICTION_TYPE="rf_classification" FG_MASK_DIR="compressed_masks" FOV_MASK_DIR="" CLASS_DIR="label_masks" NUM_TREES=10 IMAGE_TYPE="real" DATA_ROOT="scratch/assure/" IMAGE_ROOT="375_125" NUM_SAMPLES=100000 MODEL_ROOT="models" MAX_N_IMAGES=400 SHIFT_IMAGES=400 qsub -V -l twoday -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="{'dt'}" NUM_LEVELS=6 WIN_SIZE=1 OUTPUT_TYPE="class_label" PREDICTION_TYPE="rf_classification" FG_MASK_DIR="compressed_masks" FOV_MASK_DIR="" CLASS_DIR="label_masks" NUM_TREES=10 IMAGE_TYPE="real" DATA_ROOT="scratch/assure/" IMAGE_ROOT="375_125" NUM_SAMPLES=100000 MODEL_ROOT="models" MAX_N_IMAGES=400 SHIFT_IMAGES=300 qsub -V -l twoday -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="{'dt'}" NUM_LEVELS=6 WIN_SIZE=1 OUTPUT_TYPE="class_label" PREDICTION_TYPE="rf_classification" FG_MASK_DIR="compressed_masks" FOV_MASK_DIR="" CLASS_DIR="label_masks" NUM_TREES=10 IMAGE_TYPE="real" DATA_ROOT="scratch/assure/" IMAGE_ROOT="375_125" NUM_SAMPLES=100000 MODEL_ROOT="models" MAX_N_IMAGES=400 SHIFT_IMAGES=200 qsub -V -l twoday -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="{'dt'}" NUM_LEVELS=6 WIN_SIZE=1 OUTPUT_TYPE="class_label" PREDICTION_TYPE="rf_classification" FG_MASK_DIR="compressed_masks" FOV_MASK_DIR="" CLASS_DIR="label_masks" NUM_TREES=10 IMAGE_TYPE="real" DATA_ROOT="scratch/assure/" IMAGE_ROOT="375_125" NUM_SAMPLES=100000 MODEL_ROOT="models" MAX_N_IMAGES=400 SHIFT_IMAGES=100 qsub -V -l twoday -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="{'dt'}" NUM_LEVELS=6 WIN_SIZE=1 OUTPUT_TYPE="class_label" PREDICTION_TYPE="rf_classification" FG_MASK_DIR="compressed_masks" FOV_MASK_DIR="" CLASS_DIR="label_masks" NUM_TREES=10 IMAGE_TYPE="real" DATA_ROOT="scratch/assure/" IMAGE_ROOT="375_125" NUM_SAMPLES=100000 MODEL_ROOT="models" MAX_N_IMAGES=400 SHIFT_IMAGES=0 qsub -V -l twoday -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
%%
DATA_ROOT="scratch/assure/models/" MODEL_ROOT="class_label/rf_classification" MODEL_PATH="908684" MAKE_SAMPLED_MAPS=0 qsub -V -l short -hold_jid 908684 matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
DATA_ROOT="scratch/assure/models/" MODEL_ROOT="class_label/rf_classification" MODEL_PATH="908685" MAKE_SAMPLED_MAPS=0 qsub -V -l short -hold_jid 908685 matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
DATA_ROOT="scratch/assure/models/" MODEL_ROOT="class_label/rf_classification" MODEL_PATH="908686" MAKE_SAMPLED_MAPS=0 qsub -V -l short -hold_jid 908686 matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
DATA_ROOT="scratch/assure/models/" MODEL_ROOT="class_label/rf_classification" MODEL_PATH="908687" MAKE_SAMPLED_MAPS=0 qsub -V -l short -hold_jid 908687 matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
DATA_ROOT="scratch/assure/models/" MODEL_ROOT="class_label/rf_classification" MODEL_PATH="908688" MAKE_SAMPLED_MAPS=0 qsub -V -l short -hold_jid 908688 matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
%%
DATA_ROOT="scratch/assure/" MODEL_ROOT="models/class_label/rf_classification" MODEL_PATH="908684" NUM_JOBS=50 IMAGE_ROOT="375_125" USE_SAMPLED_MAPS=0 MASK_DIR="compressed_masks" qsub -V -t 2-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
DATA_ROOT="scratch/assure/" MODEL_ROOT="models/class_label/rf_classification" MODEL_PATH="908685" NUM_JOBS=50 IMAGE_ROOT="375_125" USE_SAMPLED_MAPS=0 MASK_DIR="compressed_masks" qsub -V -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
DATA_ROOT="scratch/assure/" MODEL_ROOT="models/class_label/rf_classification" MODEL_PATH="908686" NUM_JOBS=50 IMAGE_ROOT="375_125" USE_SAMPLED_MAPS=0 MASK_DIR="compressed_masks" qsub -V -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
DATA_ROOT="scratch/assure/" MODEL_ROOT="models/class_label/rf_classification" MODEL_PATH="908687" NUM_JOBS=50 IMAGE_ROOT="375_125" USE_SAMPLED_MAPS=0 MASK_DIR="compressed_masks" qsub -V -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
DATA_ROOT="scratch/assure/" MODEL_ROOT="models/class_label/rf_classification" MODEL_PATH="908688" NUM_JOBS=50 IMAGE_ROOT="375_125" USE_SAMPLED_MAPS=0 MASK_DIR="compressed_masks" qsub -V -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
%%
data_list_dir = 'C:\isbe\density\assure\375_125\data_lists\';
load([data_list_dir 'case_fold_names.mat'], 'case_fold_names', 'control_fold_names', 'case_fold_names');
[~, case_fold_idx] = sort(case_fold_names);

prediction_dir = 'C:\isbe\density\assure\375_125\predictions\class_label\rf_classification\';
fold_dirs = {
    [prediction_dir '908684\']
    [prediction_dir '908685\']
    [prediction_dir '908686\']
    [prediction_dir '908687\']
    [prediction_dir '908688\']};

%%
prediction_scores = zeros(length(case_fold_names),1);
for i_fold = 1:5
    for i_case = 1:100
        i_idx = 100*(i_fold-1) + i_case;
        i_case = case_fold_idx(i_idx);
        compressed_mask = u_load([compressed_mask_dir case_fold_names{i_case} '_c_mask.mat']);
        prediction_image = u_load([fold_dirs{i_fold} case_fold_names{i_case} '_pred.mat']);
        prediction_scores(i_case) = mean(prediction_image(compressed_mask));
        display(prediction_scores(i_case));
    end
    %figure; 
    %subplot(1,2,1); imgray(compressed_mask);
    %subplot(1,2,2); imgray(prediction_image);
end
%%
ims1 = dir([prediction_dir '908684\fold_1*.mat']);
ims2 = dir([prediction_dir '908685\fold_1*.mat']);
%%
scores = zeros(100,2);
for i_case = 1:100
    p1 = u_load([prediction_dir '908684\' ims1(i_case).name]);
    p2 = u_load([prediction_dir '908685\' ims2(i_case).name]);
    i_case = case_fold_idx(i_case);
    compressed_mask = u_load([compressed_mask_dir case_fold_names{i_case} '_c_mask.mat']);
    
    scores(i_case,1) = mean(p1(compressed_mask));
    scores(i_case,2) = mean(p2(compressed_mask));
%     figure; 
%     subplot(1,2,1); imgray(compressed_mask);
%     subplot(1,2,2); imgray(p1);
end
%%
score_labels = [false(375,1); true(125,1)];
[roc_pts, auc, ~, ~, auc_se] = calculate_roc_curve(prediction_scores,score_labels(randperm(500)));

figure; plot(roc_pts(:,1), roc_pts(:,2), '-x');
axis([0 1 0 1]); axis equal;
title(['New model: AUC = ' num2str(auc, 2) ' \pm ' num2str(auc_se,2)]);
for i_case = 1:10
    [roc_pts, auc, ~, ~, auc_se] = calculate_roc_curve(prediction_scores,score_labels(randperm(500)));

    figure; plot(roc_pts(:,1), roc_pts(:,2), '-x');
    axis([0 1 0 1]); axis equal;
    title(['New model random permutations: AUC = ' num2str(auc, 2) ' \pm ' num2str(auc_se,2)]);
end
%%
previous_scores = zeros(500,1);
%score_dir = 'A:\D_2.2_Datasets\SCORES\SCORES_VDM\';
%scores_dir = 'A:\D_2.2_Datasets\SCORES\SCORES_COMBINED\125_375_VDM\';
scores_dir = 'A:\D_2.2_Datasets\SCORES\125_AM\SCORES_VDM\';
scores_list = dir([scores_dir '*.mat']);
for i_case = 1:500
    score = load([scores_dir scores_list(i_case).name]);
    previous_scores(i_case) = score.score_class_1;
end
%%
score_labels = [true(125,1); false(375,1)];
[roc_pts, auc, ~, ~, auc_se] = calculate_roc_curve(previous_scores,score_labels);

figure; plot(roc_pts(:,1), roc_pts(:,2), '-x');
axis([0 1 0 1]); axis equal;
title(['Old model: AUC = ' num2str(auc, 2) ' \pm ' num2str(auc_se,2)]);
    
for i_case = 1:10
    [roc_pts, auc, ~, ~, auc_se] = calculate_roc_curve(previous_scores,score_labels(randperm(500)));

    figure; plot(roc_pts(:,1), roc_pts(:,2), '-x');
    axis([0 1 0 1]); axis equal;
    title(['Old model random permutations AUC = ' num2str(auc, 2) ' \pm ' num2str(auc_se,2)]);
end
%%
figure; hold all;
grid_pts = linspace(30, 55, 100);

kdist = build_1d_kernel_distribution(previous_scores, grid_pts, 0);
plot(kdist.x, kdist.D_f, 'b', 'linewidth', 2);
kdist = build_1d_kernel_distribution(100*prediction_scores, grid_pts, 0);
plot(kdist.x, kdist.D_f, 'r', 'linewidth', 2);
%%
figure; hold all;
grid_pts = linspace(30, 55, 100);

kdist = build_1d_kernel_distribution(previous_scores(126:500), grid_pts, 0);
plot(kdist.x, kdist.D_f, 'b', 'linewidth', 2);
kdist = build_1d_kernel_distribution(previous_scores(1:125), grid_pts, 0);
plot(kdist.x, kdist.D_f, 'r', 'linewidth', 2);
title('Previous model')
%%
figure; hold all;
grid_pts = linspace(30, 55, 100);
kdist = build_1d_kernel_distribution(100*prediction_scores(1:375), grid_pts, 0);
plot(kdist.x, kdist.D_f, 'b', 'linewidth', 2);
kdist = build_1d_kernel_distribution(100*prediction_scores(376:500), grid_pts, 0);
plot(kdist.x, kdist.D_f, 'r', 'linewidth', 2);
title('New model');
%%

%score_dir = 'A:\D_2.2_Datasets\SCORES\SCORES_VDM\';
%scores_dir = 'A:\D_2.2_Datasets\SCORES\SCORES_COMBINED\125_375_VDM\';
scores_dir = 'A:\D_8.3_Datasets\SCORES\SCORES_VDM\';
scores_list = dir([scores_dir '*.mat']);
num_cases = length(scores_list);
previous_scores = zeros(num_cases,1);
case_labels = false(num_cases,1);
for i_case = 1:num_cases
    score = load([scores_dir scores_list(i_case).name]);
    previous_scores(i_case) = score.score_class_1;
    
    if ~isempty(strfind(scores_list(i_case).name, 'Cancer'))
        case_labels(i_case) = 1;
    end
end
%%
[roc_pts, auc, ~, ~, auc_se] = calculate_roc_curve(previous_scores, case_labels);

figure; plot(roc_pts(:,1), roc_pts(:,2), '-x');
axis([0 1 0 1]); axis equal;
title(['Old model: AUC = ' num2str(auc, 2) ' \pm ' num2str(auc_se,2)]);
%%
decomposition_args.decomp_type = {'dt'};      %Use the dual-tree
decomposition_args.win_size = 1;            %Window size about pixel to sample features
decomposition_args.levels = 1:6;            %Levels to include
decomposition_args.feature_shape = 'rect';  %Keep this as rect, although there are other things you can play with (see sample_dt_data)
decomposition_args.feature_type = 'conj';   %Again, 'conj' works best, but can use 'mag' for magnitude only, or 'phase' (see convert_complex_representation)
decomposition_args.do_max = 0;              %Keep all 6 sub-bands (0) or only keep the maximum (1)
decomposition_args.rotate = 0;              %Try and make rotation invariant (best left switched off)
decomposition_args.use_nag = 0;             %Use the NAG toolbox to interpolate (best left swicthed off)
decomposition_args.normalise = 0;

rf_mat = u_load('A:\D_2.6_Datasets\Nijmegen_21_1_2015_20_RF_MB\rf.mat');

%image_dir = 'A:\D_8.3_Datasets\ASSURE_PROCAS_300CancerVolpara\';
 
%mam_vdm = imread([image_dir 'PROCAS_CANCERS_00001\PROCAS_CANCERS-00001-LMLO115232-20121113-RAW_hint_densityMap_H_36.67mm.pgm']);
%mam_mask = imread([image_dir 'PROCAS_CANCERS_00001\PROCAS_CANCERS-00001-LMLO115232-20121113-RAW_hint_segmentation.pgm']);
mam_vdm = u_load('C:\isbe\density\assure\procas318\images\cancer0320.mat');
mam_vdm = (double(mam_vdm)./ 65535 - 0.5) .* 2.0;
mam_mask = u_load('C:\isbe\density\assure\procas318\compressed_masks\cancer0320_c_mask.mat');
%%
display('Class prediction');
tic;
[prediction_image_c] = predict_image(... % non-strict mode
    'image_in', mam_vdm, ... % the mandatory arguments
    'decomposition_args', decomposition_args, ...
    'predictor', rf_mat,...
    'prediction_type', 'rf_classification',...
    'output_type', 'detection',...
    'use_probs', 0,...
    'mask', mam_mask,...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 1024,...
    'incremental_results', 0);
toc;
%
display('Regression prediction');
tic;
[prediction_image_r] = predict_image(... % non-strict mode
    'image_in', mam_vdm, ... % the mandatory arguments
    'decomposition_args', decomposition_args, ...
    'predictor', rf_mat,...
    'prediction_type', 'rf_regression',...
    'output_type', 'detection',...
    'use_probs', 0,...
    'mask', mam_mask,...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 1024,...
    'incremental_results', 0);
toc;
%%
test_options.predict_all = 1;
dual_tree = compute_dual_tree(mam_vdm, 6);

[rows cols] = find(mask == 4);
[responses] = compute_filter_responses(mam_vdm, decomposition_args);
sampled_features = sample_image_features(responses, rows, cols, decomposition_args);
votes = zeros(size(rows,1),2);
for i_model = 1:20
    rf = u_load(['A:\D_2.6_Datasets\Nijmegen_21_1_2015_20_Model_' num2str(i_model) '.mat']);
    [~, votes_i] = classRF_predict(sampled_features, rf, test_options);
    votes = votes + votes_i;
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%

image_dir = 'c:\isbe\density\assure\procas318\images\';
full_mask_dir = 'c:\isbe\density\assure\procas318\full_masks\';
compressed_mask_dir = 'c:\isbe\density\assure\procas318\compressed_masks\';
label_mask_dir = 'c:\isbe\density\assure\procas318\label_masks\';
prediction_dir = 'C:\isbe\density\assure\procas318\predictions\detection\rf_classification\Nijmegen_21_1_2015_20_RF_MB\';
datalist_dir = 'c:\isbe\density\assure\procas318\data_lists\';
%
load([datalist_dir 'case_names.mat'], 'im_names');
%%
keep = true(num_cases,1);
case_scores_comp = zeros(num_cases,1);
case_scores_5000 = zeros(num_cases,1);
case_labels = false(num_cases,1);

scores_hist_control = zeros(201,1);
scores_hist_cancer = zeros(201,1);
hist_bins = linspace(0, 1, 201);

for i_case = 1:num_cases
    if isempty(im_names{i_case}) || ~exist([prediction_dir im_names{i_case} '_pred.mat'], 'file')
        keep(i_case) = 0;
    else
        full_mask = u_load([full_mask_dir im_names{i_case} '_f_mask.mat']);
        prediction_image = u_load([prediction_dir im_names{i_case} '_pred.mat']);

        case_scores_comp(i_case) = mean(prediction_image(full_mask==4));
        
        scores_hist_i = hist(prediction_image(full_mask==4), hist_bins);
        if strcmpi(im_names{i_case}(1:6), 'cancer')
            case_labels(i_case) = 1;
            scores_hist_cancer = scores_hist_cancer + scores_hist_i(:);
        else
            scores_hist_control = scores_hist_control + scores_hist_i(:);
        end
    end
    
    if isempty(case_score_names{i_case})
        keep(i_case) = 0;
    else
        score = load(case_score_names{i_case});
        case_scores_5000(i_case) = score.score_class_2;
    end
end

scores_sum_control = sum(scores_hist_control);
scores_sum_cancer = sum(scores_hist_cancer);

scores_hist_control = scores_hist_control / scores_sum_control;
scores_hist_cancer = scores_hist_cancer / scores_sum_cancer;
%%
[roc_pts_5000, auc_5000, ~, ~, auc_se_5000] = calculate_roc_curve(case_scores_5000(keep), case_labels(keep));
[roc_pts_comp, auc_comp, ~, ~, auc_se_comp] = calculate_roc_curve(case_scores_comp(keep), case_labels(keep));

figure; 
subplot(1,2,1); 
plot(roc_pts_5000(:,1), roc_pts_5000(:,2), '-x');
axis equal; axis([0 1 0 1]); 
title(['Emmanouil model: AUC = ' num2str(auc_5000, 2) ' \pm ' num2str(auc_se_5000,2)]);

subplot(1,2,2);
plot(roc_pts_comp(:,1), roc_pts_comp(:,2), '-x');
axis equal; axis([0 1 0 1]); 
title(['Matlab Nijmegen model: AUC = ' num2str(auc_comp, 2) ' \pm ' num2str(auc_se_comp,2)]);    
%%
figure; 
subplot(1,2,1); axis equal; hold all; 
title('Scatter plot of scores from EM and MB models');
plot(case_scores_5000(keep & ~case_labels), case_scores_comp(keep & ~case_labels)*100, 'bx');
plot(case_scores_5000(keep & case_labels), case_scores_comp(keep & case_labels)*100, 'rx');
plot([30 50], [30 50], 'k--');
xlabel('Mean scores of 5000 randomly selected pixels (EM model)');
ylabel('Mean scores of compressed region (MB model)');
axis([30 50 30 50]);

subplot(1,2,2); hold all; 
title('Bland-Altman plot of scores from EM and MB models');
plot(...
    (case_scores_5000(keep & ~case_labels)+case_scores_comp(keep & ~case_labels)*100)/2,...
    case_scores_5000(keep & ~case_labels) - case_scores_comp(keep & ~case_labels)*100, 'bx');
plot(...
    (case_scores_5000(keep & case_labels)+case_scores_comp(keep & case_labels)*100)/2,...
    case_scores_5000(keep & case_labels) - case_scores_comp(keep & case_labels)*100, 'rx');

plot([30 50], [0 0], 'k--');
xlabel('Mean MB and EM');
ylabel('EM - MB');
axis([30 50 -2.5 2.5]);
%%
figure;
hold all;
plot(hist_bins, scores_hist_control, 'b', 'linewidth', 2);
plot(hist_bins, scores_hist_cancer, 'r', 'linewidth', 2);
title('Distributions of indiviudal pixel scores');
legend({'Controls', 'Cancers'});
%%
mean_hist_bins = linspace(0.3,0.6,21);
mean_score_hist_control = hist(case_scores_comp(keep & ~case_labels), mean_hist_bins);
mean_score_hist_cancer = hist(case_scores_comp(keep & case_labels), mean_hist_bins);

figure; hold all;
plot(mean_hist_bins, mean_score_hist_control / num_controls, 'b', 'linewidth', 2);
plot(mean_hist_bins, mean_score_hist_cancer / num_cancers, 'r', 'linewidth', 2);
title('Histograms of case scores');
legend({'Controls', 'Cancers'});
%
grid_pts = linspace(0.3, 0.6, 101);
kdist_control = build_1d_kernel_distribution(case_scores_comp(keep & ~case_labels), grid_pts, 0);
kdist_cancer = build_1d_kernel_distribution(case_scores_comp(keep & case_labels), grid_pts, 0);

figure; hold all;
plot(kdist_control.x, kdist_control.D_f, 'b', 'linewidth', 2);
plot(kdist_cancer.x, kdist_cancer.D_f, 'r', 'linewidth', 2);
title('Kernel estimated distributions of case scores');
legend({'Controls', 'Cancers'});
%%
all_scores = zeros(scores_sum_control + scores_sum_cancer, 1, 'uint8');
all_labels = false(scores_sum_control + scores_sum_cancer, 1);
curr_sample = 0;
for i_case = 1:num_cases
    if ~isempty(im_names{i_case})
        full_mask = u_load([full_mask_dir im_names{i_case} '_f_mask.mat']);
        prediction_image = u_load([prediction_dir im_names{i_case} '_pred.mat']);

        scores_i = uint8(200*prediction_image(full_mask==4));
        idx = curr_sample + (1:length(scores_i));
        curr_sample = idx(end);
        
        all_scores(idx) = scores_i;
        if strcmpi(im_names{i_case}(1:6), 'cancer')
            all_labels(idx) = 1;
        end
    end
end
[roc_pts_all, auc_all, ~, ~, auc_se_all] = calculate_roc_curve(all_scores, all_labels, 0:200);
%%
scores_5000 = zeros(1268*5000, 1, 'uint8');
labels_5000 = false(1268*5000, 1);
case_scores_5000_new = zeros(num_cases,1);

curr_sample = 0;
for i_case = 1:num_cases
    if ~isempty(im_names{i_case})
        full_mask = u_load([full_mask_dir im_names{i_case} '_f_mask.mat']);
        prediction_image = u_load([prediction_dir im_names{i_case} '_pred.mat']);

        scores_i = prediction_image(full_mask==4);
        r_idx = randperm(length(scores_i), 5000);
        
        idx = curr_sample + (1:5000);
        curr_sample = idx(end);
        
        scores_5000(idx) = uint8(200*scores_i(r_idx));
        case_scores_5000_new(i_case) = mean(scores_i(r_idx));
        
        if strcmpi(im_names{i_case}(1:6), 'cancer')
            labels_5000(idx) = 1;
        end
    end
end
%%
[roc_pts_5000a, auc_5000a, ~, ~, auc_se_5000a] = calculate_roc_curve(scores_5000, labels_5000, 0:200);
figure; 
subplot(1,2,1);
plot(roc_pts_all(:,1), roc_pts_all(:,2), '-x');
axis equal; axis([0 1 0 1]); 
title(['Matlab model - all pixels pooled: AUC = ' num2str(auc_all, 2) ' \pm ' num2str(auc_se_all,2)]);

subplot(1,2,2); 
plot(roc_pts_5000a(:,1), roc_pts_5000a(:,2), '-x');
axis equal; axis([0 1 0 1]); 
title(['Matlab model - 5000 pixels per image: AUC = ' num2str(auc_5000a, 2) ' \pm ' num2str(auc_se_5000a,2)]);

%%
[roc_pts_5000n, auc_5000n, ~, ~, auc_se_5000n] = calculate_roc_curve(case_scores_5000_new(keep), case_labels(keep));

figure; 
subplot(1,2,1); 
plot(roc_pts_5000(:,1), roc_pts_5000(:,2), '-x');
axis equal; axis([0 1 0 1]); 
title(['Emmanouil model: AUC = ' num2str(auc_5000, 2) ' \pm ' num2str(auc_se_5000,2)]);

subplot(1,2,2);
plot(roc_pts_5000n(:,1), roc_pts_5000n(:,2), '-x');
axis equal; axis([0 1 0 1]); 
title(['Mean of random 5000 pixels: AUC = ' num2str(auc_5000n, 2) ' \pm ' num2str(auc_se_5000n,2)]);
%%
figure; 
subplot(1,2,1); axis equal; hold all; 
title('Scatter plot of scores from EM and MB models');
plot(case_scores_5000(keep & ~case_labels), case_scores_5000_new(keep & ~case_labels)*100, 'bx');
plot(case_scores_5000(keep & case_labels), case_scores_5000_new(keep & case_labels)*100, 'rx');
plot([30 50], [30 50], 'k--');
xlabel('Mean scores of 5000 randomly selected pixels (EM model)');
ylabel('Mean scores of compressed region (MB model)');
axis([30 50 30 50]);

subplot(1,2,2); hold all; 
title('Bland-Altman plot of scores from EM and MB models');
plot(...
    (case_scores_5000(keep & ~case_labels)+case_scores_5000_new(keep & ~case_labels)*100)/2,...
    case_scores_5000(keep & ~case_labels) - case_scores_5000_new(keep & ~case_labels)*100, 'bx');
plot(...
    (case_scores_5000(keep & case_labels)+case_scores_5000_new(keep & case_labels)*100)/2,...
    case_scores_5000(keep & case_labels) - case_scores_5000_new(keep & case_labels)*100, 'rx');

plot([30 50], [0 0], 'k--');
xlabel('Mean MB and EM');
ylabel('EM - MB');
axis([30 50 -2.5 2.5]);
%%
[max_scores, max_cases] = sort(case_scores_5000_new, 'descend');
%%
for i_im = 501:520%1:20
    i_case = max_cases(i_im);
    mam_vdm = u_load([image_dir im_names{i_case} '.mat']);
    mask = u_load([full_mask_dir im_names{i_case} '_f_mask.mat']);    
    prediction_image = u_load([prediction_dir im_names{i_case} '_pred.mat']);    
    b_inner = bwboundaries(mask == 4);
    
    pred_smooth = imfilter(prediction_image, fspecial('gauss', 27, 8));
    pred_smooth(mask~=4) = 0;
    glims = [min(mam_vdm(mask==4)) max(mam_vdm(mask==4))];
    
    figure; 
    subplot(1,2,1); 
    imgray(mam_vdm);
    plot(b_inner{1}(:,2), b_inner{1}(:,1), 'r--', 'linewidth', 2);
    caxis(glims);
    
    subplot(1,2,2);
    %imgray(prediction_image);
    %plot(b_inner{1}(:,2), b_inner{1}(:,1), 'r--', 'linewidth', 2);
    imgray(pred_smooth); caxis([0.2 0.6]);
    title(im_names{i_case});
end
%%
for i_im = 61:80
    i_case = max_cases(end - i_im - 1);
    mam_vdm = u_load([image_dir im_names{i_case} '.mat']);
    mask = u_load([full_mask_dir im_names{i_case} '_f_mask.mat']);
    prediction_image = u_load([prediction_dir im_names{i_case} '_pred.mat']);
    b_inner = bwboundaries(mask == 4);
    
    pred_smooth = imfilter(prediction_image, fspecial('gauss', 27, 8));
    pred_smooth(mask~=4) = 0;
    glims = [min(mam_vdm(mask==4)) max(mam_vdm(mask==4))];
    
    figure; 
    subplot(1,2,1); 
    imgray(mam_vdm);
    plot(b_inner{1}(:,2), b_inner{1}(:,1), 'r--', 'linewidth', 2);
    caxis(glims);
    
    subplot(1,2,2);
    %imgray(prediction_image);
    %plot(b_inner{1}(:,2), b_inner{1}(:,1), 'r--', 'linewidth', 2);
    imgray(pred_smooth); caxis([0.2 0.6]);
    title(im_names{i_case});
end
%%
[roc_pts_sort_t, auc_sort_t, ~, ~, auc_se_sort_t] = calculate_roc_curve(max_scores(1:634), sorted_labels(1:634));
[roc_pts_sort_b, auc_sort_b, ~, ~, auc_se_sort_b] = calculate_roc_curve(max_scores(634:1268), sorted_labels(634:1268));
figure; 
subplot(1,2,1);
plot(roc_pts_sort_b(:,1), roc_pts_sort_b(:,2), '-x'); hold on;
plot([0 1], [0 1], 'k');
axis equal; axis([0 1 0 1]); 
title(['Matlab model - cases sorted, bottom half: AUC = ' num2str(auc_sort_b, 2) ' \pm ' num2str(auc_se_sort_b,2)]);

subplot(1,2,2);
plot(roc_pts_sort_t(:,1), roc_pts_sort_t(:,2), '-x'); hold on;
plot([0 1], [0 1], 'k');
axis equal; axis([0 1 0 1]); 
title(['Matlab model - cases sorted, top half: AUC = ' num2str(auc_sort_t, 2) ' \pm ' num2str(auc_se_sort_t,2)]);
%%
vas_scores = xlsread('A:\PROCAS_CASE_CONTROL_SET\CC VAS and cumulus results.xls', 1, 'I2:I1267');
vas_labels = [true(317,1); false(949,1)];
[roc_pts_vas, auc_vas, ~, ~, auc_se_vas] = calculate_roc_curve(vas_scores, vas_labels);
[roc_pts_vdm, auc_vdm, ~, ~, auc_se_vdm] = calculate_roc_curve(g_means(keep), case_labels(keep));
figure; 
subplot(1,2,1);
plot(roc_pts_vas(:,1), roc_pts_vas(:,2), '-x'); hold on;
plot([0 1], [0 1], 'k');
axis equal; axis([0 1 0 1]); 
title(['VAS single reader LML: AUC = ' num2str(auc_vas, 2) ' \pm ' num2str(auc_se_vas,2)]);

subplot(1,2,2);
plot(roc_pts_vdm(:,1), roc_pts_vdm(:,2), '-x'); hold on;
plot([0 1], [0 1], 'k');
axis equal; axis([0 1 0 1]); 
title(['VDM mean of compressed: AUC = ' num2str(auc_vdm, 2) ' \pm ' num2str(auc_se_vdm,2)]);
%%
%%
for i_im = 101:120
    i_case = max_cases(i_im);
    mam_vdm = u_load([image_dir im_names{i_case} '.mat']);
    mask = u_load([full_mask_dir im_names{i_case} '_f_mask.mat']);    
    b_inner = bwboundaries(mask == 4);
    glims = [min(mam_vdm(mask==4)) max(mam_vdm(mask==4))];
    
    figure; 
    subplot(1,2,1); 
    imgray(mam_vdm);
    plot(b_inner{1}(:,2), b_inner{1}(:,1), 'r--', 'linewidth', 2);
    caxis(glims);
    title(im_names{i_case});
    
    i_case = max_cases(end - 1 - i_im);
    mam_vdm = u_load([image_dir im_names{i_case} '.mat']);
    mask = u_load([full_mask_dir im_names{i_case} '_f_mask.mat']);    
    b_inner = bwboundaries(mask == 4);
    glims = [min(mam_vdm(mask==4)) max(mam_vdm(mask==4))];
    
    subplot(1,2,2); 
    imgray(mam_vdm);
    plot(b_inner{1}(:,2), b_inner{1}(:,1), 'r--', 'linewidth', 2);
    caxis(glims);
    title(im_names{i_case});
end
%%
g_lims_all = zeros(num_cases,2);
g_means = zeros(num_cases,1);
for i_case = 1:num_cases
    if keep(i_case)
        mam_vdm = u_load([image_dir im_names{i_case} '.mat']);
        mask = u_load([full_mask_dir im_names{i_case} '_f_mask.mat']);    
        glims = [min(mam_vdm(mask==4)) max(mam_vdm(mask==4))];
        g_lims_all(i_case,:) = glims;
        g_means(i_case) = mean(mam_vdm(mask==4)); 
    end
end
glims_20_80 = [prctile(g_lims_all(keep,1), 20) prctile(g_lims_all(keep,2), 80)];
%%
[sorted_density sorted_density_idx] = sort(g_means, 'descend');
sorted_density_labels = case_labels(sorted_density_idx);
sorted_cancer_idx = sorted_density_idx(sorted_density_labels);
sorted_control_idx = sorted_density_idx(~sorted_density_labels);
sorted_density_scores = case_scores_comp(sorted_density_idx);

%%
num_cols = 6;
num_rows = 5;

curr_control = 1;
curr_cancer = 1;
%%
create_folder('C:\isbe\density\assure\tiled_density_ims\gray\');
create_folder('C:\isbe\density\assure\tiled_density_ims\jet\');
for i_fig = 1:21
    figure('windowstyle', 'docked', 'units', 'normalized',...
        'position', [0 0 1 1], 'menubar', 'none');
    
    idx_fig = [...
        sorted_control_idx(curr_control + (0:3*num_cols*num_rows/2-1));...
        sorted_cancer_idx(curr_cancer + (0:num_cols*num_rows/2-1))];
    
    g_lims_fig = g_lims_all(idx_fig,:);
    glims_20_80 = [prctile(g_lims_fig(:,1), 40) prctile(g_lims_fig(:,2), 60)];
    
    for i_row = 1:num_rows;
        for i_col = 1:(num_cols/2);

            i_case = sorted_control_idx(curr_control);
            curr_control = curr_control + 3;
            mam_vdm = u_load([image_dir im_names{i_case} '.mat']);
            mask = u_load([full_mask_dir im_names{i_case} '_f_mask.mat']);
            mam_vdm(mask==3) = nan;
            
            pos_x = (i_col - 1)/num_cols;
            pos_y = (i_row - 1)/num_rows;
            axes('units', 'normalized', 'position', [pos_x pos_y 1/num_cols 1/num_rows]);     
            imgray(rot90(mam_vdm)); axis off; 
            caxis(glims_20_80);
            %caxis(g_lims_all(i_case,:));         
        end
        for i_col = 1+(num_cols/2):num_cols;

            i_case = sorted_cancer_idx(curr_cancer);
            curr_cancer = curr_cancer + 1;
            mam_vdm = u_load([image_dir im_names{i_case} '.mat']);
            mask = u_load([full_mask_dir im_names{i_case} '_f_mask.mat']);
            mam_vdm(mask==3) = nan;
            
            pos_x = (i_col - 1)/num_cols;
            pos_y = (i_row - 1)/num_rows;
            axes('units', 'normalized', 'position', [pos_x pos_y 1/num_cols 1/num_rows]);     
            imgray(rot90(mam_vdm)); axis off;
            caxis(glims_20_80);
            %caxis(g_lims_all(i_case,:));
            
        end
    end
    exportfig(['C:\isbe\density\assure\tiled_density_ims\gray\im' zerostr(i_fig,2) '.png']);
    delete(['C:\isbe\density\assure\tiled_density_ims\gray\im' zerostr(i_fig,2) '.eps']);
    delete(['C:\isbe\density\assure\tiled_density_ims\gray\im' zerostr(i_fig,2) '.pdf']);
    colormap(jet(256));
    exportfig(['C:\isbe\density\assure\tiled_density_ims\jet\im' zerostr(i_fig,2) '.png']);
    delete(['C:\isbe\density\assure\tiled_density_ims\jet\im' zerostr(i_fig,2) '.eps']);
    delete(['C:\isbe\density\assure\tiled_density_ims\jet\im' zerostr(i_fig,2) '.pdf']);
end
%%
figure; hold on;
plot(g_means(keep & ~case_labels), case_scores_comp(keep & ~case_labels), 'bx');
plot(g_means(keep & case_labels), case_scores_comp(keep & case_labels), 'rx');
xlabel('Mean VDM of compressed region');
ylabel('Mean texture score of compressed region');
title('Scatter plot of density vs Manchester texture');

[roc_pts_sort_t, auc_sort_t, ~, ~, auc_se_sort_t] = calculate_roc_curve(sorted_density_scores(1:634), sorted_density_labels(1:634));
[roc_pts_sort_b, auc_sort_b, ~, ~, auc_se_sort_b] = calculate_roc_curve(sorted_density_scores(634:1268), sorted_density_labels(634:1268));
figure; 
subplot(1,2,1);
plot(roc_pts_sort_b(:,1), roc_pts_sort_b(:,2), '-x'); hold on;
plot([0 1], [0 1], 'k');
axis equal; axis([0 1 0 1]); 
title(['Matlab model - cases sorted by density, bottom half: AUC = ' num2str(auc_sort_b, 2) ' \pm ' num2str(auc_se_sort_b,2)]);

subplot(1,2,2);
plot(roc_pts_sort_t(:,1), roc_pts_sort_t(:,2), '-x'); hold on;
plot([0 1], [0 1], 'k');
axis equal; axis([0 1 0 1]); 
title(['Matlab model - cases sorted density, top half: AUC = ' num2str(auc_sort_t, 2) ' \pm ' num2str(auc_se_sort_t,2)]);