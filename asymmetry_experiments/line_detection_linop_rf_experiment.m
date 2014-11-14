%Combine together random forests from SAN
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% 1x1 windows
mb_combine_rfs('rf_dir', '\\isbe-san1\mberks\dev\classification\rf\',...
    'tree_dir', 'C:\isbe\dev\classification\rf\rf_linop_05\', 'delete_trees', 0,...
    'replace_tree_root', '\\isbe-san1\mberks\dev\',...
    'save_path', 'C:\isbe\dev\classification\rf\rf_linop_05\');
%%
mb_combine_rfs('rf_dir', '\\isbe-san1\mberks\dev\classification\rf\rf_linop_1_5_8\',...
    'tree_dir', '\\isbe-san1\mberks\dev\classification\rf\rf_linop_1_5_8\', 'delete_trees', 1,...
    'replace_tree_root', '\\isbe-san1\mberks\dev\',...
    'save_path', '\\isbe-san1\mberks\dev\classification\rf\rf_linop_1_5_8\');
%%
random_forest = u_load('C:\isbe\dev\classification\rf\rf_linop_1_5_8\random_forest.mat');
param_dir = 'rf_linop_1_5_8';

test_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin';
prob_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin\probability_images\';
mkdir([prob_dir param_dir]);

for ii = 1:100
    test_image = load([test_dir '\image' zerostr(ii,3) '.mat']);

    %Levels 1, min wave length 1, Onf 0.65
    prob_im = classify_image_linop(...
        'image_in', test_image.image,...
        'forest', random_forest,...
        'forest_type', 'isbe',...
        'num_levels', 5,...
        'win_size', 1);
    save([prob_dir param_dir '\probability_image' zerostr(ii,3) '.mat'], 'prob_im');

end
%%
test_dir = 'C:\isbe\dev\classification\data\test_28\';
prob_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin\probability_images\';
results_dir = 'C:\isbe\dev\classification\line_detection_results\';

[roc_a, auc_a, tp_counts_a, fp_counts_a, t_counts_a, f_counts_a] =...
    compute_roc_image_set_lines(test_dir, [prob_dir 'rf_linop_1_5_8\'], 'centre_not_all');

[roc_b, auc_b, tp_counts_b, fp_counts_b, t_counts_b, f_counts_b] =...
    compute_roc_image_set_lines(test_dir, [prob_dir 'rf_linop_1_5_8\'], 'centre_line');

[roc_c, auc_c, tp_counts_c, fp_counts_c, t_counts_c, f_counts_c] =...
    compute_roc_image_set_lines(test_dir, [prob_dir 'rf_linop_1_5_8\'], 'all_line');

save([results_dir 'rf_linop_1_5_8.mat'], '*_a', '*_b', '*_c');
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% 3x3 windows
mb_combine_rfs('rf_dir', '\\isbe-san1\mberks\dev\classification\rf\',...
    'tree_dir', 'C:\isbe\dev\classification\rf\rf_linop_3_5_8\', 'delete_trees', 0,...
    'replace_tree_root', '\\isbe-san1\mberks\dev\',...
    'save_path', 'C:\isbe\dev\classification\rf\rf_linop_3_5_8\');
%%
mb_combine_rfs('rf_dir', '\\isbe-san1\mberks\dev\classification\rf\',...
    'tree_dir', '\\isbe-san1\mberks\dev\classification\rf\rf_linop_3_5_8\', 'delete_trees', 1,...
    'replace_tree_root', '\\isbe-san1\mberks\dev\',...
    'save_path', '\\isbe-san1\mberks\dev\classification\rf\rf_linop_3_5_8\');
%%
random_forest = u_load('C:\isbe\dev\classification\rf\rf_linop_1_5_8\random_forest.mat');
param_dir = 'rf_linop_1_5_8';

test_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin';
prob_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin\probability_images\';
mkdir([prob_dir param_dir]);

for ii = 1:100
    test_image = load([test_dir '\image' zerostr(ii,3) '.mat']);

    %Levels 1, min wave length 1, Onf 0.65
    prob_im = classify_image_linop(...
        'image_in', test_image.image,...
        'forest', random_forest,...
        'forest_type', 'isbe',...
        'num_levels', 5,...
        'win_size', 1);
    save([prob_dir param_dir '\probability_image' zerostr(ii,3) '.mat'], 'prob_im');

end
%%
test_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin';
prob_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin\probability_images\';
results_dir = 'C:\isbe\dev\classification\line_detection_results\';

[roc_a, auc_a, tp_counts_a, fp_counts_a, t_counts_a, f_counts_a] =...
    compute_roc_image_set_lines(test_dir, [prob_dir 'rf_linop_1_5_8\'], 'centre_not_all');

[roc_b, auc_b, tp_counts_b, fp_counts_b, t_counts_b, f_counts_b] =...
    compute_roc_image_set_lines(test_dir, [prob_dir 'rf_linop_1_5_8\'], 'centre_line');

[roc_c, auc_c, tp_counts_c, fp_counts_c, t_counts_c, f_counts_c] =...
    compute_roc_image_set_lines(test_dir, [prob_dir 'rf_linop_1_5_8\'], 'all_line');

save([results_dir 'rf_linop_1_5_8.mat'], '*_a', '*_b', '*_c');