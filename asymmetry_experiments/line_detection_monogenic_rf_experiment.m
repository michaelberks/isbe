%Monogenic/RF experiment script - add more commments

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
test_dir = 'M:\chen\data\testimage_contrast1to8_exprnd_sin\'; %ZC changes these
prob_dir = 'M:\chen\data\testimage_contrast1to8_exprnd_sin\probability_image\';
results_dir = 'C:\isbe\dev\classification\line_detection_results\';
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

for win_size = [1 3]
    for level = 3:6
        param_dir = ['monogenic_rf_trees_W' num2str(win_size) 'L' num2str(level)];       
        [roc_a, auc_a, tp_counts_a, fp_counts_a, t_counts_a, f_counts_a] =...
            compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all'); %#ok
        [roc_b, auc_b, tp_counts_b, fp_counts_b, t_counts_b, f_counts_b] =...
            compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_line'); %#ok
        [roc_c, auc_c, tp_counts_c, fp_counts_c, t_counts_c, f_counts_c] =...
            compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all'); %#ok
        save([results_dir param_dir '.mat'], '*_a', '*_b', '*_c');
    end
end