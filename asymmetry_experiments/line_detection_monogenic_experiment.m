%Monogenic experiment script - add more commments

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% test_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin'; %ZC changes these
% prob_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin\probability_images\';
% results_dir = 'C:\isbe\dev\classification\line_detection_results\';
test_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\';
prob_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\';
results_dir = 'C:\isbe\asymmetry_project\data\line_detection_rfs\results\';
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%1. Varying number of scales and wavelengths
for wave = 8%[8 16]
    %figure; hold all; axis equal; axis([0 1 0 1]);
    %legend_text = cell(5,1);
    for level = 2%1:4

        param_dir = ['monogenic_0p65_' zerostr(wave,2) '_' zerostr(level,2)];
        mkdir([prob_dir param_dir]);

        for ii = 1:100
            load([test_dir '\image' zerostr(ii,3) '.mat']);
    
            %Levels 1, min wave length 1, Onf 0.65
            prob_im = monogenic_phase_cong(test_image, level, wave, 2, 0.65);
            save([prob_dir param_dir '\probability_image' zerostr(ii,3) '.mat'], 'prob_im');
    
        end
        
        [roc_a, auc_a, tp_counts_a, fp_counts_a, t_counts_a, f_counts_a] =...
            compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all');
        [roc_b, auc_b, tp_counts_b, fp_counts_b, t_counts_b, f_counts_b] =...
            compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_line');
        [roc_c, auc_c, tp_counts_c, fp_counts_c, t_counts_c, f_counts_c] =...
            compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\'], 'centre_not_all');
        save([results_dir param_dir '.mat'], '*_a', '*_b', '*_c');
        %plot(roc_b(:,1), roc_b(:,2), '.'); legend_text{level} = [num2str(level) 'scale(s), A_z = ', num2str(auc_b, 3)];
            
    end
    %legend(legend_text, 'location', 'southeast');
    %title(['ROC curves for a monogenic signal based line detector, minimum wavelength set to ' num2str(wave)])
end
%%
test_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin'; %ZC changes these
prob_dir = 'C:\isbe\dev\classification\data\testimage_contrast1to8_exprnd_sin\probability_images\';
results_dir = 'C:\isbe\dev\classification\line_detection_results\';
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%1. Varying number of scales and wavelengths
for wave = [8 16]
    figure(...
        'WindowStyle', 'normal', ...
        'units', 'pixels',...
        'position', [10 10 800 800]);
    hold all; axis equal; axis([0 1 0 1]);
    legend_text = cell(4,1);
    for level = 1:4

        param_dir = ['monogenic_0p65_' zerostr(wave,2) '_' zerostr(level,2)];
        load([results_dir param_dir '.mat'], '*roc_b', '*auc_b');
        plot(roc_b(:,1), roc_b(:,2), '.'); legend_text{level} = [num2str(level) 'scale(s), A_z = ', num2str(auc_b, 3)];
            
    end
    legend(legend_text, 'location', 'southeast');
    title(['ROC curves for a monogenic signal based line detector, minimum wavelength set to ' num2str(wave)])
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% let's actually look at some line prob maps to check this is working
%% ok...
if 0
    for ii = 1:5
        test = load([test_dir '\image' zerostr(ii,3) '.mat']);

        %Levels 3, min wave length 3, Onf 0.65
        load([prob_dir 'monogenic_0p65_02_03\probability_image' zerostr(ii,3) '.mat'], 'prob_im');
        figure; 
        subplot(1,2,1); imagesc(test.image); axis image; colormap(gray(256));
        subplot(1,2,2); imagesc(prob_im); axis image; colormap(gray(256));

        %Levels 4, min wave length 4, Onf 0.65
        load([prob_dir 'monogenic_0p65_02_04\probability_image' zerostr(ii,3) '.mat'], 'prob_im');
        figure; 
        subplot(1,2,1); imagesc(test.image); axis image; colormap(gray(256));
        subplot(1,2,2); imagesc(prob_im); axis image; colormap(gray(256));

        %Levels 5, min wave length 5, Onf 0.65
        load([prob_dir 'monogenic_0p65_02_05\probability_image' zerostr(ii,3) '.mat'], 'prob_im');
        figure; 
        subplot(1,2,1); imagesc(test.image); axis image; colormap(gray(256));
        subplot(1,2,2); imagesc(prob_im); axis image; colormap(gray(256));


    end
end
%%
test_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\';
prob_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\';
results_dir = 'C:\isbe\asymmetry_project\data\line_detection_rfs\results\';
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%1. Varying number of scales and wavelengths
for wave = 8%[8 16]
    %figure; hold all; axis equal; axis([0 1 0 1]);
    %legend_text = cell(5,1);
    for level = 2%1:4

        param_dir = ['monogenic_0p65_' zerostr(wave,2) '_' zerostr(level,2) '_ori'];
        mkdir([prob_dir param_dir]);

        for ii = 1:100
            load([test_dir '\image' zerostr(ii,3) '.mat']);
    
            %Levels 1, min wave length 1, Onf 0.65
            [d d orientation_image] = monogenic_phase_cong(test_image, level, wave, 2, 0.65);
            save([prob_dir param_dir '\probability_image' zerostr(ii,3) '.mat'], 'orientation_image');
    
        end
        
        [ori_err_mono] = compute_image_orientation_errors(test_dir, [prob_dir param_dir '\'], 'centre_line');
        save([results_dir param_dir '.mat'], 'ori_err_mono');
            
    end
end