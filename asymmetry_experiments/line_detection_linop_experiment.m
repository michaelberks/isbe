%Linop experiment script - add more commments

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
test_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\';
prob_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\';
results_dir = 'C:\isbe\asymmetry_project\data\line_detection_rfs\results\';
%
% for angles = [8 12]
% %     figure; hold all; axis equal; axis([0 1 0 1]);
% %     legend_text = cell(5,1);
%     for level = 4:6
% 
%         param_dir = ['linop_octave_' zerostr(level,2) '_' zerostr(angles,2)];
%         mkdir([prob_dir param_dir]);
% 
%         for ii = 1:100
%             load([test_dir '\image' zerostr(ii,3) '.mat']);
%     
%             %Levels 1, min wave length 1, Onf 0.65
%             prob_im = line_operator_octave(image, angles, level, 'degrees', 0);
%             save([prob_dir param_dir '\probability_image' zerostr(ii,3) '.mat'], 'prob_im');
%     
%         end
%         
%         [roc_a, auc_a, tp_counts_a, fp_counts_a, t_counts_a, f_counts_a] =...
%             compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all');
%         [roc_b, auc_b, tp_counts_b, fp_counts_b, t_counts_b, f_counts_b] =...
%             compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_line');
%         [roc_c, auc_c, tp_counts_c, fp_counts_c, t_counts_c, f_counts_c] =...
%             compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all');
%         save([results_dir param_dir '.mat'], '*_a', '*_b', '*_c');
% %         plot(roc_b(:,1), roc_b(:,2), '.'); legend_text{1} = [num2str(level) 'scale(s), A_z = ', num2str(auc_b, 3)];
%             
%     end
% %     legend(legend_text, 'location', 'southeast');
% %     title(['ROC curves for a linop line detector (octave scales), minimum wavelength set to ' num2str(wave)])
% end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
for angles = [8]
%     figure; hold all; axis equal; axis([0 1 0 1]);
%     legend_text = cell(5,1);
    for level = [53]

        param_dir = ['linop_odd_' zerostr(level,2) '_' zerostr(angles,2)];
        mkdir([prob_dir param_dir]);

        for ii = 1:100
            load([test_dir '\image' zerostr(ii,3) '.mat']);
    
            %Levels 1, min wave length 1, Onf 0.65
            prob_im = line_operator_conv(test_image, angles, 3, level, 'degrees');
            save([prob_dir param_dir '\probability_image' zerostr(ii,3) '.mat'], 'prob_im');
    
        end
        
        [roc_a, auc_a, tp_counts_a, fp_counts_a, t_counts_a, f_counts_a] =...
            compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all');
        [roc_b, auc_b, tp_counts_b, fp_counts_b, t_counts_b, f_counts_b] =...
            compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_line');
        [roc_c, auc_c, tp_counts_c, fp_counts_c, t_counts_c, f_counts_c] =...
            compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all');
        save([results_dir param_dir '.mat'], '*_a', '*_b', '*_c');
%         plot(roc_b(:,1), roc_b(:,2), '.'); legend_text{1} = [num2str(level) 'scale(s), A_z = ', num2str(auc_b, 3)];
            
    end
%     legend(legend_text, 'location', 'southeast');
%     title(['ROC curves for a linop line detector (octave scales), minimum wavelength set to ' num2str(wave)])
end
%%

test_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin'; %ZC changes these
prob_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin\probability_images\';
results_dir = 'C:\isbe\dev\classification\line_detection_results\';

for angles = [12]
    figure(...
        'WindowStyle', 'normal', ...
        'units', 'pixels',...
        'position', [10 10 800 800]);
    hold all; axis equal; axis([0 1 0 1]);
    legend_text = cell(3,1);
    for level = 4:6

        param_dir = ['linop_octave_' zerostr(level,2) '_' zerostr(angles,2)];
        load([results_dir param_dir '.mat'], '*roc_b', '*auc_b');
        plot(roc_b(:,1), roc_b(:,2), '.'); legend_text{level-3} = [num2str(level) 'scale(s), A_z = ', num2str(auc_b, 3)];

    end
    legend(legend_text, 'location', 'southeast');
    title(['ROC curves for a linop, number of angles set to ' num2str(angles)])
end
%%
test_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin'; %ZC changes these
prob_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin\probability_images\';
results_dir = 'C:\isbe\dev\classification\line_detection_results\';

for angles = [8 12]
    figure(...
        'WindowStyle', 'normal', ...
        'units', 'pixels',...
        'position', [10 10 800 800]);
    hold all; axis equal; axis([0 1 0 1]);
    legend_text = cell(5,1);
    ii = 1;
    for level = [21 29 37 53 85]

        param_dir = ['linop_odd_' zerostr(level,2) '_' zerostr(angles,2)];
        load([results_dir param_dir '.mat'], '*roc_b', '*auc_b');
        plot(roc_b(:,1), roc_b(:,2), '.'); legend_text{ii} = ['scales 3 to ' num2str(level) ', A_z = ', num2str(auc_b, 3)];
        ii = ii+1;
    end
    legend(legend_text, 'location', 'southeast');
    title(['ROC curves for a linop, number of angles set to ' num2str(angles)])
end
