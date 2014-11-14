%Dual tree experiment script - add more commments

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
test_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin'; %ZC changes these
prob_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin\probability_images\';
results_dir = 'C:\isbe\dev\classification\line_detection_results\';
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

param_dir = '182392';
mb_name = 'dt_3_5_all_pm';       
[roc_a, auc_a, tp_counts_a, fp_counts_a, t_counts_a, f_counts_a] =...
    compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all'); %#ok
[roc_b, auc_b, tp_counts_b, fp_counts_b, t_counts_b, f_counts_b] =...
    compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_line'); %#ok
[roc_c, auc_c, tp_counts_c, fp_counts_c, t_counts_c, f_counts_c] =...
    compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all'); %#ok
save([results_dir mb_name '.mat'], '*_a', '*_b', '*_c');
plot(roc_b(:,1), roc_b(:,2), '.'); legend_text{1} = ['3x3, 5 levels, phase/mag, A_z = ', num2str(auc_b, 3)];

param_dir = '182396';
mb_name = 'dt_1_5_all_pm';       
[roc_a, auc_a, tp_counts_a, fp_counts_a, t_counts_a, f_counts_a] =...
    compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all'); %#ok
[roc_b, auc_b, tp_counts_b, fp_counts_b, t_counts_b, f_counts_b] =...
    compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_line'); %#ok
[roc_c, auc_c, tp_counts_c, fp_counts_c, t_counts_c, f_counts_c] =...
    compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all'); %#ok
save([results_dir mb_name '.mat'], '*_a', '*_b', '*_c');

param_dir = '182417';
mb_name = 'dt_3_5_all_rr';       
[roc_a, auc_a, tp_counts_a, fp_counts_a, t_counts_a, f_counts_a] =...
    compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all'); %#ok
[roc_b, auc_b, tp_counts_b, fp_counts_b, t_counts_b, f_counts_b] =...
    compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_line'); %#ok
[roc_c, auc_c, tp_counts_c, fp_counts_c, t_counts_c, f_counts_c] =...
    compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all'); %#ok
save([results_dir mb_name '.mat'], '*_a', '*_b', '*_c');

param_dir = '182418';
mb_name = 'dt_3_5_all_mm';       
[roc_a, auc_a, tp_counts_a, fp_counts_a, t_counts_a, f_counts_a] =...
    compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all'); %#ok
[roc_b, auc_b, tp_counts_b, fp_counts_b, t_counts_b, f_counts_b] =...
    compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_line'); %#ok
[roc_c, auc_c, tp_counts_c, fp_counts_c, t_counts_c, f_counts_c] =...
    compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all'); %#ok
save([results_dir mb_name '.mat'], '*_a', '*_b', '*_c');

param_dir = '182419';
mb_name = 'dt_3_5_all_pp';       
[roc_a, auc_a, tp_counts_a, fp_counts_a, t_counts_a, f_counts_a] =...
    compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all'); %#ok
[roc_b, auc_b, tp_counts_b, fp_counts_b, t_counts_b, f_counts_b] =...
    compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_line'); %#ok
[roc_c, auc_c, tp_counts_c, fp_counts_c, t_counts_c, f_counts_c] =...
    compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all'); %#ok
save([results_dir mb_name '.mat'], '*_a', '*_b', '*_c');

%%
test_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin'; %ZC changes these
prob_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin\probability_images\';
results_dir = 'C:\isbe\dev\classification\line_detection_results\';

figure(...
    'WindowStyle', 'normal', ...
    'units', 'pixels',...
    'position', [10 10 800 800]);
hold all; axis equal; axis([0 1 0 1]);
legend_text = cell(5,1);

%
load([results_dir 'dt_3_5_all_pm.mat'], '*roc_b', '*auc_b');
plot(roc_b(:,1), roc_b(:,2), '.'); legend_text{1} = ['3x3, 5 levels, phase/mag, A_z = ', num2str(auc_b, 3)];
%
load([results_dir 'dt_1_5_all_pm.mat'], '*roc_b', '*auc_b');
plot(roc_b(:,1), roc_b(:,2), '.'); legend_text{2} = ['1x1, 5 levels, phase/mag, A_z = ', num2str(auc_b, 3)];
%
load([results_dir 'dt_3_5_all_rr.mat'], '*roc_b', '*auc_b');
plot(roc_b(:,1), roc_b(:,2), '.'); legend_text{3} = ['3x3, 5 levels, real, A_z = ', num2str(auc_b, 3)];
%
load([results_dir 'dt_3_5_all_mm.mat'], '*roc_b', '*auc_b');
plot(roc_b(:,1), roc_b(:,2), '.'); legend_text{4} = ['3x3, 5 levels, mag, A_z = ', num2str(auc_b, 3)];
%
load([results_dir 'dt_3_5_all_pp.mat'], '*roc_b', '*auc_b');
plot(roc_b(:,1), roc_b(:,2), '.'); legend_text{5} = ['3x3, 5 levels, phase, A_z = ', num2str(auc_b, 3)];


title('ROC curves for a dual-tree based line detector')
legend(legend_text, 'location', 'southeast');