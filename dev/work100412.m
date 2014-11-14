%--------------------------------------------------------------------------
% Check our method for computing an error trace as we predict a test set
% from a random forest

X_test = load('C:\isbe\dev\classification\data\zip_data\zip.test');
X_train = load('C:\isbe\dev\classification\data\zip_data\zip.train');

forest_args.sampling_method = 'bootstrap_train_data';
forest_args.sampling_method_args = X_train;
forest_args.d = 16;
forest_args.tree_dir = 'C:\isbe\dev\classification\rf\rf_mb_test_data\';
forest_args.X_test = X_test(:,2:end);
forest_args.y_test = X_test(:,1);
forest_args.n_trees = 200;
         
[random_forest] = mb_random_forest_class_train(forest_args);
[y_fit votes error_trace] = mb_random_forest_class_predict(random_forest, X_test(:,2:end), X_test(:,1));
%%
%--------------------------------------------------------------------------
% Check that sampling 1x1 patches is the same as sampling 3x3 by and
% indexing the centre pixel (i.e. check we've got our indexing correct)
rand('twister', 5489);
sampling_method_args.bg_dir = 'M:\chen\data\smooth512x512_patches\train\';
sampling_method_args.save_path = [];
sampling_method_args.num_samples = 1e3;
sampling_method_args.width_range = [0 8];
sampling_method_args.contrast_range = [0 8];
sampling_method_args.win_size = 3;
sampling_method_args.num_levels = 6;
[X_test3 y_test3] = sample_training_data_pixelset(sampling_method_args);

rand('twister', 5489);
sampling_method_args.bg_dir = 'M:\chen\data\smooth512x512_patches\train\';
sampling_method_args.save_path = [];
sampling_method_args.num_samples = 1e3;
sampling_method_args.width_range = [0 8];
sampling_method_args.contrast_range = [0 8];
sampling_method_args.win_size = 3;
sampling_method_args.num_levels = 5;
[X_test3_5 y_test3_5] = sample_training_data_pixelset(sampling_method_args);

rand('twister', 5489);
sampling_method_args.bg_dir = 'M:\chen\data\smooth512x512_patches\train\';
sampling_method_args.save_path = [];
sampling_method_args.num_samples = 1e3;
sampling_method_args.width_range = [0 8];
sampling_method_args.contrast_range = [0 8];
sampling_method_args.win_size = 1;
sampling_method_args.num_levels = 6;

[X_test1 y_test1] = sample_training_data_pixelset(sampling_method_args);

X_test1a = X_test3(:,5:9:end);
display(max(abs(X_test1a(:) - X_test1(:))));
display(max(abs(y_test3(:) - y_test1(:))));

X_test3a = X_test3(:,[1:6*5*9 6*6*9+1:11*6*9]);
display(max(abs(X_test3a(:) - X_test3_5(:))));
display(max(abs(y_test3(:) - y_test3_5(:))));
%%
%--------------------------------------------------------------------------
% Sample a set of test data to compare our line detection random forests
sampling_method_args.bg_dir = 'M:\chen\data\smooth512x512_patches\train\';
sampling_method_args.save_path = [];
sampling_method_args.num_samples = 5e4;
sampling_method_args.width_range = [0 8];
sampling_method_args.contrast_range = [0 8];
sampling_method_args.win_size = 3;
sampling_method_args.bar_type = 'sin';

[X_test y_test] = sample_training_data_pixelset(sampling_method_args);
save C:\isbe\dev\classification\line_detection_results\test_data.mat X_test y_test
%%
% Compute forest error traces for 3x3 and 1x1 forests
forest_1 = u_load('M:\chen\data\rf_trees_182320\tree_combine\random_forest.mat');
forest_3 = u_load('M:\chen\data\rf_trees_182321\tree_combine\random_forest.mat');
%
[y_fit votes1 error_trace1] = mb_random_forest_class_predict(forest_1, X_test(:,5:9:end), y_test);
[y_fit votes3 error_trace3] = mb_random_forest_class_predict(forest_3, X_test, y_test);
%
[roc_pts1 auc1] = calculate_roc_curve(votes1(:,2)/400, y_test,(-1:400)/400);
[roc_pts3 auc3] = calculate_roc_curve(votes3(:,2)/400, y_test,(-1:400)/400);
%
f1 = figure; hold all;
axis equal; axis([0 1 0 1]);
plot(roc_pts1(:,1), roc_pts1(:,2), '-x');
plot(roc_pts3(:,1), roc_pts3(:,2), '-x');
legend({['1x1, A_z = ' num2str(auc1, 3)], ['3x3, A_z = ' num2str(auc3, 3)]}, 'location', 'southeast');
title('ROC curves with operating points shown for forests constructed using 1x1 and 3x3 patches - original test data');
xlabel('Specificity');
ylabel('Sensitivity');
saveas(f1, 'C:\isbe\dev\classification\line_detection_results\figures\original_roc_curves.bmp');
saveas(f1, 'C:\isbe\dev\classification\line_detection_results\figures\original_roc_curves.fig');
%
f2 = figure; hold all;
plot(1:400, error_trace1);
plot(1:400, error_trace3);
legend({'1x1' , '3x3'}, 'location', 'southeast');
plot([0 400], [max(error_trace1(:)) max(error_trace1(:))], 'k--');
plot([0 400], [max(error_trace3(:)) max(error_trace3(:))], 'k--');
title('ROC A_z for increasing number of trees in forests constructed using 1x1 and 3x3 patches - original test data');
xlabel('No. of trees');
ylabel('A_z');
saveas(f2, 'C:\isbe\dev\classification\line_detection_results\figures\original_roc_traces.bmp');
saveas(f2, 'C:\isbe\dev\classification\line_detection_results\figures\original_roc_traces.fig');
%%

%%
counts1 = hist(votes1(:,2)/400, (1:200)/200);
counts3 = hist(votes3(:,2)/400, (1:200)/200);
counts_real = hist(1-predict(1).probability_image(:), (1:200)/200);
%%
consensus_prob = (votes1(:,2) + votes3(:,2)) / 800;
%
bins = (0:200)/200;
[counts_real] = histc(prob_image(:), bins);

test_idx = zeros(sum(counts_real), 1);

curr_pt = 1;
for ii = 1:length(counts_real)
    if ii == length(counts_real)
        idx = find(consensus_prob == bins(ii));
    else
        idx = find(consensus_prob >= bins(ii) & (consensus_prob < bins(ii+1)));
    end
    n_pts = counts_real(ii);
    
    if n_pts
        test_idx(curr_pt:curr_pt+n_pts-1) = idx(ceil(length(idx)*rand(n_pts,1)));
    end
    curr_pt = curr_pt + n_pts;
end
test_idx_small = randsample(test_idx, 50000);
save C:\isbe\dev\classification\line_detection_results\test_data.mat X_test y_test test_idx test_idx_small
%%
forest_1 = u_load('M:\chen\data\rf_trees_182320\tree_combine\random_forest.mat');
forest_3 = u_load('M:\chen\data\rf_trees_182321\tree_combine\random_forest.mat');

[y_fit votes3 error_trace3] = mb_random_forest_class_predict(forest_3, X_test(test_idx_small,:), y_test(test_idx_small));
[roc_pts3 auc3] = calculate_roc_curve(votes3(:,2)/400, y_test(test_idx_small),(-1:400)/400);


[y_fit votes1 error_trace1] = mb_random_forest_class_predict(forest_1, X_test(test_idx_small,5:9:end), y_test(test_idx_small));
[roc_pts1 auc1] = calculate_roc_curve(votes1(:,2)/400, y_test(test_idx_small),(-1:400)/400);
%%
f1 = figure; hold all;
axis equal; axis([0 1 0 1]);
plot(roc_pts1(:,1), roc_pts1(:,2), '-x');
plot(roc_pts3(:,1), roc_pts3(:,2), '-x');
legend({['1x1, A_z = ' num2str(auc1, 3)], ['3x3, A_z = ' num2str(auc3, 3)]}, 'location', 'southeast');
title('ROC curves with operating points shown for forests constructed using 1x1 and 3x3 patches - modified test data');
xlabel('1 - Specificity');
ylabel('Sensitivity');
saveas(f1, 'C:\isbe\dev\classification\line_detection_results\figures\modified_roc_curves.bmp');
saveas(f1, 'C:\isbe\dev\classification\line_detection_results\figures\modified_roc_curves.fig');
%
f2 = figure; hold all;
plot(1:400, error_trace1);
plot(1:400, error_trace3);
legend({'1x1' , '3x3'}, 'location', 'southeast');
plot([0 400], [max(error_trace1(:)) max(error_trace1(:))], 'k--');
plot([0 400], [max(error_trace3(:)) max(error_trace3(:))], 'k--');
title('ROC A_z for increasing number of trees in forests constructed using 1x1 and 3x3 patches - modified test data');
xlabel('No. of trees');
ylabel('A_z');
saveas(f2, 'C:\isbe\dev\classification\line_detection_results\figures\modified_roc_traces.bmp');
saveas(f2, 'C:\isbe\dev\classification\line_detection_results\figures\modified_roc_traces.fig');
%%
% Repeat process for other forest parameters
%
load C:\isbe\dev\classification\line_detection_results\test_data.mat X_test y_test test_idx_small
%%
display('Results for 3x3 patches');
% 1. 3x3, 3 levels, all bgs (182326)
lev = 3;
forest = u_load('M:\chen\data\rf_trees_182326\tree_combine\random_forest.mat');
forest.tree_dir = 'M:\chen\data\rf_trees_182326\tree_combine\';
data_cols = [1:6*lev*9 6*6*9+1:(lev+6)*6*9];

[y_fit votes] = mb_random_forest_class_predict(forest, X_test(test_idx_small,data_cols)); clear y_fit;
[roc_pts auc] = calculate_roc_curve(votes(:,2)/200, y_test(test_idx_small),(-1:200)/200); clear votes;
save C:\isbe\dev\classification\line_detection_results\roc_3x3_3_all_modified.mat roc_pts auc;
display(['level 3 modified roc = ', num2str(auc)]);

[y_fit votes] = mb_random_forest_class_predict(forest, X_test(:,data_cols)); clear y_fit;
[roc_pts auc] = calculate_roc_curve(votes(:,2)/200, y_test(:),(-1:200)/200); clear votes;
save C:\isbe\dev\classification\line_detection_results\roc_3x3_3_all_original.mat roc_pts auc;
display(['level 3 original roc = ', num2str(auc)]);
%
% 2. 3x3, 4 levels, all bgs (182327)
lev = 4;
forest = u_load('M:\chen\data\rf_trees_182327\tree_combine\random_forest.mat');
forest.tree_dir = 'M:\chen\data\rf_trees_182327\tree_combine\';
data_cols = [1:6*lev*9 6*6*9+1:(lev+6)*6*9];

[y_fit votes] = mb_random_forest_class_predict(forest, X_test(test_idx_small,data_cols)); clear y_fit;
[roc_pts auc] = calculate_roc_curve(votes(:,2)/200, y_test(test_idx_small),(-1:200)/200); clear votes;
display(['level 4 modified roc = ', num2str(auc)]);
save C:\isbe\dev\classification\line_detection_results\roc_3x3_4_all_modified.mat roc_pts auc;

[y_fit votes] = mb_random_forest_class_predict(forest, X_test(:,data_cols)); clear y_fit;
[roc_pts auc] = calculate_roc_curve(votes(:,2)/200, y_test(:),(-1:200)/200); clear votes;
save C:\isbe\dev\classification\line_detection_results\roc_3x3_4_all_original.mat roc_pts auc;
display(['level 4 original roc = ', num2str(auc)]);

% 3. 3x3, 5 levels, all bgs (182328)
lev = 5;
forest = u_load('M:\chen\data\rf_trees_182328\tree_combine\random_forest.mat');
forest.tree_dir = 'M:\chen\data\rf_trees_182328\tree_combine\';
data_cols = [1:6*lev*9 6*6*9+1:(lev+6)*6*9];

[y_fit votes] = mb_random_forest_class_predict(forest, X_test(test_idx_small,data_cols)); clear y_fit;
[roc_pts auc] = calculate_roc_curve(votes(:,2)/200, y_test(test_idx_small),(-1:200)/200); clear votes;
save C:\isbe\dev\classification\line_detection_results\roc_3x3_5_all_modified.mat roc_pts auc;
display(['level 5 modified roc = ', num2str(auc)]);

[y_fit votes] = mb_random_forest_class_predict(forest, X_test(:,data_cols)); clear y_fit;
[roc_pts auc] = calculate_roc_curve(votes(:,2)/200, y_test(:),(-1:200)/200); clear votes;
save C:\isbe\dev\classification\line_detection_results\roc_3x3_5_all_original.mat roc_pts auc;
display(['level 5 original roc = ', num2str(auc)]);

%
% 7. 3x3, 6 levels, all bgs (182272)
lev = 6;
forest = u_load('M:\chen\data\rf_trees_182270\tree_combine\random_forest.mat');
forest.tree_dir = 'M:\chen\data\rf_trees_182270\tree_combine\';
data_cols = 1:12*6*9;

[y_fit votes] = mb_random_forest_class_predict(forest, X_test(test_idx_small,data_cols)); clear y_fit;
[roc_pts auc] = calculate_roc_curve(votes(:,2)/200, y_test(test_idx_small),(-1:200)/200); clear votes;
save C:\isbe\dev\classification\line_detection_results\roc_3x3_6_all_modified.mat roc_pts auc;
display(['level 6 modified roc = ', num2str(auc)]);

[y_fit votes] = mb_random_forest_class_predict(forest, X_test(:,data_cols)); clear y_fit;
[roc_pts auc] = calculate_roc_curve(votes(:,2)/200, y_test(:),(-1:200)/200); clear votes;
save C:\isbe\dev\classification\line_detection_results\roc_3x3_6_all_original.mat roc_pts auc;
display(['level 6 original roc = ', num2str(auc)]);

display('Results for 1x1 patches');
%
% 4. 1x1, 3 levels, all bgs (182329)
lev = 3;
forest = u_load('M:\chen\data\rf_trees_182329\tree_combine\random_forest.mat');
forest.tree_dir = 'M:\chen\data\rf_trees_182329\tree_combine\';
data_cols = [5:9:6*lev*9 6*6*9+6:9:(lev+6)*6*9];

[y_fit votes] = mb_random_forest_class_predict(forest, X_test(test_idx_small,data_cols)); clear y_fit;
[roc_pts auc] = calculate_roc_curve(votes(:,2)/200, y_test(test_idx_small),(-1:200)/200); clear votes;
save C:\isbe\dev\classification\line_detection_results\roc_1x1_3_all_modified.mat roc_pts auc;
display(['level 3 modified roc = ', num2str(auc)]);

[y_fit votes] = mb_random_forest_class_predict(forest, X_test(:,data_cols)); clear y_fit;
[roc_pts auc] = calculate_roc_curve(votes(:,2)/200, y_test(:),(-1:200)/200); clear votes;
save C:\isbe\dev\classification\line_detection_results\roc_1x1_3_all_original.mat roc_pts auc;
display(['level 3 original roc = ', num2str(auc)]);

% 5. 1x1, 4  levels, all bgs (182330)
lev = 4;
forest = u_load('M:\chen\data\rf_trees_182330\tree_combine\random_forest.mat');
forest.tree_dir = 'M:\chen\data\rf_trees_182330\tree_combine\';
data_cols = [5:9:6*lev*9 6*6*9+6:9:(lev+6)*6*9];

[y_fit votes] = mb_random_forest_class_predict(forest, X_test(test_idx_small,data_cols)); clear y_fit;
[roc_pts auc] = calculate_roc_curve(votes(:,2)/200, y_test(test_idx_small),(-1:200)/200); clear votes;
save C:\isbe\dev\classification\line_detection_results\roc_1x1_4_all_modified.mat roc_pts auc;
display(['level 4 modified roc = ', num2str(auc)]);

[y_fit votes] = mb_random_forest_class_predict(forest, X_test(:,data_cols)); clear y_fit;
[roc_pts auc] = calculate_roc_curve(votes(:,2)/200, y_test(:),(-1:200)/200); clear votes;
save C:\isbe\dev\classification\line_detection_results\roc_1x1_4_all_original.mat roc_pts auc;
display(['level 4 original roc = ', num2str(auc)]);

% 6. 1x1, 5 levels, all bgs (182331)
lev = 5;
forest = u_load('M:\chen\data\rf_trees_182331\tree_combine\random_forest.mat');
forest.tree_dir = 'M:\chen\data\rf_trees_182331\tree_combine\';
data_cols = [5:9:6*lev*9 6*6*9+6:9:(lev+6)*6*9];

[y_fit votes] = mb_random_forest_class_predict(forest, X_test(test_idx_small,data_cols)); clear y_fit;
[roc_pts auc] = calculate_roc_curve(votes(:,2)/200, y_test(test_idx_small),(-1:200)/200); clear votes;
save C:\isbe\dev\classification\line_detection_results\roc_1x1_5_all_modified.mat roc_pts auc;
display(['level 5 modified roc = ', num2str(auc)]);

[y_fit votes] = mb_random_forest_class_predict(forest, X_test(:,data_cols)); clear y_fit;
[roc_pts auc] = calculate_roc_curve(votes(:,2)/200, y_test(:),(-1:200)/200); clear votes;
save C:\isbe\dev\classification\line_detection_results\roc_1x1_5_all_original.mat roc_pts auc;
display(['level 5 original roc = ', num2str(auc)]);
%
% 7. 1x1, 6 levels, all bgs (182272)
lev = 6;
forest = u_load('M:\chen\data\rf_trees_182272\tree_combine\random_forest.mat');
forest.tree_dir = 'M:\chen\data\rf_trees_182272\tree_combine\';
data_cols = [5:9:6*lev*9 6*6*9+6:9:(lev+6)*6*9];

[y_fit votes] = mb_random_forest_class_predict(forest, X_test(test_idx_small,data_cols)); clear y_fit;
[roc_pts auc] = calculate_roc_curve(votes(:,2)/200, y_test(test_idx_small),(-1:200)/200); clear votes;
save C:\isbe\dev\classification\line_detection_results\roc_1x1_6_all_modified.mat roc_pts auc;
display(['level 6 modified roc = ', num2str(auc)]);

[y_fit votes] = mb_random_forest_class_predict(forest, X_test(:,data_cols)); clear y_fit;
[roc_pts auc] = calculate_roc_curve(votes(:,2)/200, y_test(:),(-1:200)/200); clear votes;
save C:\isbe\dev\classification\line_detection_results\roc_1x1_6_all_original.mat roc_pts auc;
display(['level 6 original roc = ', num2str(auc)]);

%%
f = figure; hold all; axis equal; axis([0 1 0 1]);

load C:\isbe\dev\classification\line_detection_results\roc_1x1_3_all_original.mat
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{1} = ['1x1, 3 levels, Az = ', num2str(auc, 3)];

load C:\isbe\dev\classification\line_detection_results\roc_1x1_4_all_original.mat
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{2} = ['1x1, 4 levels, Az = ', num2str(auc, 3)];

load C:\isbe\dev\classification\line_detection_results\roc_1x1_5_all_original.mat
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{3} = ['1x1, 5 levels, Az = ', num2str(auc, 3)];

load C:\isbe\dev\classification\line_detection_results\roc_1x1_6_all_original.mat
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{4} = ['1x1, 6 levels, Az = ', num2str(auc, 3)];

title('ROC curves, for forests constructed from 1x1 patches, varying no. of levels - original test data')
legend(legend_text, 'location', 'southeast');
xlabel('1 - Specificity');
ylabel('Sensitivity');

saveas(f, 'C:\isbe\dev\classification\line_detection_results\figures\original_roc_curves_1x1_levels.fig');
saveas(f, 'C:\isbe\dev\classification\line_detection_results\figures\original_roc_curves_1x1_levels.bmp'); clear f;
%
f = figure; hold all; axis equal; axis([0 1 0 1]);

load C:\isbe\dev\classification\line_detection_results\roc_1x1_3_all_modified.mat
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{1} = ['1x1, 3 levels, Az = ', num2str(auc, 3)];

load C:\isbe\dev\classification\line_detection_results\roc_1x1_4_all_modified.mat
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{2} = ['1x1, 4 levels, Az = ', num2str(auc, 3)];

load C:\isbe\dev\classification\line_detection_results\roc_1x1_5_all_modified.mat
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{3} = ['1x1, 5 levels, Az = ', num2str(auc, 3)];

load C:\isbe\dev\classification\line_detection_results\roc_1x1_6_all_modified.mat
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{4} = ['1x1, 6 levels, Az = ', num2str(auc, 3)];

title('ROC curves, for forests constructed from 1x1 patches, varying no. of levels - modified test data')
legend(legend_text, 'location', 'southeast');
xlabel('1 - Specificity');
ylabel('Sensitivity');

saveas(f, 'C:\isbe\dev\classification\line_detection_results\figures\modified_roc_curves_1x1_levels.fig');
saveas(f, 'C:\isbe\dev\classification\line_detection_results\figures\modified_roc_curves_1x1_levels.bmp'); clear f;
%%
f = figure; hold all; axis equal; axis([0 1 0 1]);

load C:\isbe\dev\classification\line_detection_results\roc_3x3_3_all_original.mat
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{1} = ['3x3, 3 levels, Az = ', num2str(auc, 3)];

load C:\isbe\dev\classification\line_detection_results\roc_3x3_4_all_original.mat
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{2} = ['3x3, 4 levels, Az = ', num2str(auc, 3)];

load C:\isbe\dev\classification\line_detection_results\roc_3x3_5_all_original.mat
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{3} = ['3x3, 5 levels, Az = ', num2str(auc, 3)];

load C:\isbe\dev\classification\line_detection_results\roc_3x3_6_all_original.mat
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{4} = ['3x3, 6 levels, Az = ', num2str(auc, 3)];

title('ROC curves, for forests constructed from 3x3 patches, varying no. of levels - original test data')
legend(legend_text, 'location', 'southeast');
xlabel('1 - Specificity');
ylabel('Sensitivity');

saveas(f, 'C:\isbe\dev\classification\line_detection_results\figures\original_roc_curves_3x3_levels.fig');
saveas(f, 'C:\isbe\dev\classification\line_detection_results\figures\original_roc_curves_3x3_levels.bmp'); 
axis([0 .3 0 .3]);
saveas(f, 'C:\isbe\dev\classification\line_detection_results\figures\original_roc_curves_3x3_levels_zoom.bmp'); clear f;
%
f = figure; hold all; axis equal; axis([0 1 0 1]);

load C:\isbe\dev\classification\line_detection_results\roc_3x3_3_all_modified.mat
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{1} = ['3x3, 3 levels, Az = ', num2str(auc, 3)];

load C:\isbe\dev\classification\line_detection_results\roc_3x3_4_all_modified.mat
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{2} = ['3x3, 4 levels, Az = ', num2str(auc, 3)];

load C:\isbe\dev\classification\line_detection_results\roc_3x3_5_all_modified.mat
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{3} = ['3x3, 5 levels, Az = ', num2str(auc, 3)];

load C:\isbe\dev\classification\line_detection_results\roc_3x3_6_all_modified.mat
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{4} = ['3x3, 6 levels, Az = ', num2str(auc, 3)];

title('ROC curves, for forests constructed from 3x3 patches, varying no. of levels - modified test data')
legend(legend_text, 'location', 'southeast');
xlabel('1 - Specificity');
ylabel('Sensitivity');

saveas(f, 'C:\isbe\dev\classification\line_detection_results\figures\modified_roc_curves_3x3_levels.fig');
saveas(f, 'C:\isbe\dev\classification\line_detection_results\figures\modified_roc_curves_3x3_levels.bmp'); 
axis([0 .3 0 .3]);
saveas(f, 'C:\isbe\dev\classification\line_detection_results\figures\modified_roc_curves_3x3_levels_zoom.bmp'); clear f;
%%
%--------------------------------------------------------------------------
test_dir = 'C:\isbe\dev\classification\data\testimage_contrast0to8_multibars_sin\test_images\';
prob_dir = 'C:\isbe\dev\classification\data\testimage_contrast0to8_multibars_sin\probability_images\';
[roc_pts_3x3_6, auc_3x3_6] = compute_roc_image_set_lines(test_dir, [prob_dir '182270\'], 'centre_line');
[roc_pts_1x1_6, auc_1x1_6] = compute_roc_image_set_lines(test_dir, [prob_dir '182272\'], 'centre_line');

figure; hold all; axis equal; axis([0 1 0 1]);
plot(roc_pts_1x1_6(:,1), roc_pts_1x1_6(:,2), '-x');
plot(roc_pts_3x3_6(:,1), roc_pts_3x3_6(:,2), '-x');
legend({['1x1, A_z = ' num2str(auc_1x1_6, 3)], ['3x3, A_z = ' num2str(auc_3x3_6, 3)]}, 'location', 'southeast');
title('ROC curves with operating points shown for forests constructed using 1x1 and 3x3 patches - image test data');
xlabel('1 - Specificity');
ylabel('Sensitivity');
%%
f = figure; hold all; axis equal; axis([0 1 0 1]);

[roc_pts, auc] = compute_roc_image_set_lines(test_dir, [prob_dir '182326\'], 'centre_line');
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{1} = ['3x3, 3 levels, Az = ', num2str(auc, 3)];

[roc_pts, auc] = compute_roc_image_set_lines(test_dir, [prob_dir '182327\'], 'centre_line');
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{2} = ['3x3, 4 levels, Az = ', num2str(auc, 3)];

[roc_pts, auc] = compute_roc_image_set_lines(test_dir, [prob_dir '182328\'], 'centre_line');
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{3} = ['3x3, 5 levels, Az = ', num2str(auc, 3)];

[roc_pts, auc] = compute_roc_image_set_lines(test_dir, [prob_dir '182270\'], 'centre_line');
plot(roc_pts(:,1), roc_pts(:,2), '-x');
legend_text{4} = ['3x3, 6 levels, Az = ', num2str(auc, 3)];

title('ROC curves, for forests constructed from 3x3 patches, varying no. of levels - modified test data')
legend(legend_text, 'location', 'southeast');
xlabel('1 - Specificity');
ylabel('Sensitivity');

saveas(f, 'C:\isbe\dev\classification\line_detection_results\figures\image_roc_curves_3x3_levels.fig');
saveas(f, 'C:\isbe\dev\classification\line_detection_results\figures\image_roc_curves_3x3_levels.bmp'); clear f;
%%
mb_combine_rfs(...
    'rf_dir', 'M:\chen\data\rf_trees_182329\',...
    'tree_dir', 'M:\chen\data\rf_trees_182329\tree_combine\',...
    'tree_root',[],...
    'replace_tree_root', 'M:\chen\data\',...
    'copy_trees', 1,...% the optional arguments
    'delete_trees', 0,...
    'save_path', 'M:\chen\data\rf_trees_182329\tree_combine\');