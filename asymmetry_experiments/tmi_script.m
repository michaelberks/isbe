%Code for generating results and figures for paper submission to MIA/TMI
results_dir = 'C:\isbe\asymmetry_project\data\line_detection_rfs\results\';

mkdir C:\isbe\asymmetry_project\data\line_detection_rfs\results\

rf_id = {'191616', '191630'};
line_type = {'circle', 'line'};
for lt = 1:2
    for id = 1:2
        
        prob_dir = ['C:\isbe\asymmetry_project\data\synthetic_lines\' line_type{lt} 's512\results\' rf_id{id} '\'];
        test_dir = ['C:\isbe\asymmetry_project\data\synthetic_lines\' line_type{lt} 's512\'];
        
        save_name = [line_type{lt} 's512_' rf_id{id} '.mat'];
        
        [roc_a, auc_a, tp_counts_a, fp_counts_a, t_counts_a, f_counts_a] =...
            compute_roc_image_set_lines(test_dir, prob_dir, 'centre_not_all'); %#ok
        [roc_b, auc_b, tp_counts_b, fp_counts_b, t_counts_b, f_counts_b] =...
            compute_roc_image_set_lines(test_dir, prob_dir, 'centre_line'); %#ok
        [roc_c, auc_c, tp_counts_c, fp_counts_c, t_counts_c, f_counts_c] =...
            compute_roc_image_set_lines(test_dir, prob_dir, 'centre_not_all'); %#ok
        save([results_dir save_name], '*_a', '*_b', '*_c');
    end
end
%%
rf_id = {'191616', '191630'};
line_type = {'line', 'circle'};
results_dir = 'C:\isbe\asymmetry_project\data\line_detection_rfs\results\';

leg = {'straight lines', 'curves'};

figure; axis equal; axis([0 1 0 1]); hold all;
legend_txt = [];
for lt = 1:2
    for id = 1:2
        
        save_name = [line_type{lt} 's512_' rf_id{id} '.mat'];
        load([results_dir save_name]);
        display(['Lines: ' line_type{lt} ', training: ', rf_id{id} ', AUC = ', num2str(auc_b)]);
        legend_txt{end+1,1} = ['Training: ' leg{id} ', testing: ', leg{lt} ', AUC = ', num2str(auc_b)];
        plot(roc_b(:,1), roc_b(:,2), 'linewidth', 2);
    end
end
legend(legend_txt, 'location', 'southeast');
%%
rf_id = {'191616', '191630', '191656', '191658'};
line_type = {'line'};
results_dir = 'C:\isbe\asymmetry_project\data\line_detection_rfs\results\';

figure; axis equal; axis([0 1 0 1]); hold all;
legend_txt = [];
for lt = 1
    for id = 1:4
        
        save_name = [line_type{lt} 's512_' rf_id{id} '.mat'];
        load([results_dir save_name]);
        display(['Lines: ' line_type{lt} ', training: ', rf_id{id} ', AUC = ', num2str(auc_b)]);
        legend_txt{end+1,1} = ['Lines: ' rf_id{id} ', training: ', line_type{lt} ', AUC = ', num2str(auc_b)];
        plot(roc_b(:,1), roc_b(:,2));
    end
end
legend(legend_txt, 'location', 'southeast');