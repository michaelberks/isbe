%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%----------------------- BMVC paper script --------------------------------
%--------------------------------------------------------------------------

%Script to generate figures for BMVC 2010 paper
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%1. Create ROC curves of different line detection methods
test_dir = 'M:\chen\data\testimage_contrast1to8_exprnd_sin\';
prob_dir = 'M:\chen\data\testimage_contrast1to8_exprnd_sin\probability_image\';
results_dir = 'M:\chen\data\testimage_contrast1to8_exprnd_sin\line_detection_results\';

f = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [50 50 800 800],...
    'PaperPositionMode','auto'); 
hold all; axis equal; axis([0 1 0 1]);
color_order = get(gca, 'colororder');
set(gca, 'colororder', color_order([1 3 2 4:end],:));
legend_text = cell(5,1);

% i) Dual-tree with random forest
% Load in saved ROC curves
load([results_dir 'dt_3_5_all_pm.mat'], '*roc_b', '*auc_b', 't_counts_b', 'f_counts_b');

%Compute errors on ROC curve
n_pos = sum(t_counts_b);
n_neg = sum(f_counts_b);
q2 = sum( ...
    (roc_b(2:end,1)-roc_b(1:end-1,1)) .* ...
    (roc_b(1:end-1,2).^2 + roc_b(1:end-1,2).*(roc_b(2:end,2)-roc_b(1:end-1,2)) +...
    (roc_b(2:end,2)-roc_b(1:end-1,2)).*(roc_b(2:end,2)-roc_b(1:end-1,2))/3) );
q1 = sum( ...
    (roc_b(2:end,2)-roc_b(1:end-1,2)) .* ...
    ((1-roc_b(2:end,1)).^2 + (1-roc_b(2:end,1)).*(roc_b(2:end,1)-roc_b(1:end-1,1)) +...
    (roc_b(2:end,1)-roc_b(1:end-1,1)).*(roc_b(2:end,1)-roc_b(1:end-1,1))/3) );
auc_se = sqrt( (auc_b*(1-auc_b) + (n_pos-1)*(q1-auc_b^2) + (n_neg-1)*(q2-auc_b^2)) / (n_neg*n_pos) );
ci_95 = 1.96*auc_se;

%Plot curve
plot(roc_b(:,1), roc_b(:,2), '-', 'LineWidth', 2);
%Store legend labels
legend_text{1} = ['\fontsize{16}  DT-CWT/RF, {\bfA_z = ', num2str(auc_b, 3) '} \pm ' num2str(ci_95,2)];
%Compute sensitivity for some specificity and display result
sens = interp1(roc_b(:,1), roc_b(:,2), 0.05);
display(['DT sensitivity at 95% specificity = ', num2str(sens)]);

%
% ii) Linop with random forest
% Load in saved ROC curves
load([results_dir 'rf_linop_1_5_8.mat'], '*roc_b', '*auc_b', 't_counts_b', 'f_counts_b');

%Compute errors on ROC curve
n_pos = sum(t_counts_b);
n_neg = sum(f_counts_b);
q2 = sum( ...
    (roc_b(2:end,1)-roc_b(1:end-1,1)) .* ...
    (roc_b(1:end-1,2).^2 + roc_b(1:end-1,2).*(roc_b(2:end,2)-roc_b(1:end-1,2)) +...
    (roc_b(2:end,2)-roc_b(1:end-1,2)).*(roc_b(2:end,2)-roc_b(1:end-1,2))/3) );
q1 = sum( ...
    (roc_b(2:end,2)-roc_b(1:end-1,2)) .* ...
    ((1-roc_b(2:end,1)).^2 + (1-roc_b(2:end,1)).*(roc_b(2:end,1)-roc_b(1:end-1,1)) +...
    (roc_b(2:end,1)-roc_b(1:end-1,1)).*(roc_b(2:end,1)-roc_b(1:end-1,1))/3) );
auc_se = sqrt( (auc_b*(1-auc_b) + (n_pos-1)*(q1-auc_b^2) + (n_neg-1)*(q2-auc_b^2)) / (n_neg*n_pos) );
ci_95 = 1.96*auc_se;

%Plot curve
plot(roc_b(:,1), roc_b(:,2), '-', 'LineWidth', 2);
%Store legend labels
legend_text{2} = ['\fontsize{16}  Linop/RF, {\bfA_z = ', num2str(auc_b, 3) '} \pm ' num2str(ci_95,2)];

%Compute sensitivity for some specificity and display result
sens = interp1(roc_b(:,1), roc_b(:,2), 0.05);
display(['Linop/rf sensitivity at 95% specificity = ', num2str(sens)]);

%
% iii) Monogenic with random forest
load([results_dir 'monogenic_rf_trees_W3L5.mat'], '*roc_b', '*auc_b', 't_counts_b', 'f_counts_b');

%Compute errors on ROC curve
n_pos = sum(t_counts_b);
n_neg = sum(f_counts_b);
q2 = sum( ...
    (roc_b(2:end,1)-roc_b(1:end-1,1)) .* ...
    (roc_b(1:end-1,2).^2 + roc_b(1:end-1,2).*(roc_b(2:end,2)-roc_b(1:end-1,2)) +...
    (roc_b(2:end,2)-roc_b(1:end-1,2)).*(roc_b(2:end,2)-roc_b(1:end-1,2))/3) );
q1 = sum( ...
    (roc_b(2:end,2)-roc_b(1:end-1,2)) .* ...
    ((1-roc_b(2:end,1)).^2 + (1-roc_b(2:end,1)).*(roc_b(2:end,1)-roc_b(1:end-1,1)) +...
    (roc_b(2:end,1)-roc_b(1:end-1,1)).*(roc_b(2:end,1)-roc_b(1:end-1,1))/3) );
auc_se = sqrt( (auc_b*(1-auc_b) + (n_pos-1)*(q1-auc_b^2) + (n_neg-1)*(q2-auc_b^2)) / (n_neg*n_pos) );
ci_95 = 1.96*auc_se;

%Plot curve
plot(roc_b(:,1), roc_b(:,2), '-', 'LineWidth', 2);
%Store legend labels
legend_text{3} = ['\fontsize{16}  Monogenic/RF, {\bfA_z = ', num2str(auc_b, 3) '} \pm ' num2str(ci_95,2)];

%Compute sensitivity for some specificity and display result
sens = interp1(roc_b(:,1), roc_b(:,2), 0.05);
display(['Monogenic/rf sensitivity at 95% specificity = ', num2str(sens)]);

%
% iv) Monogenic, analytic method
load([results_dir 'monogenic_0p65_08_02.mat'], '*roc_b', '*auc_b', 't_counts_b', 'f_counts_b');

%Compute errors on ROC curve
n_pos = sum(t_counts_b);
n_neg = sum(f_counts_b);
q2 = sum( ...
    (roc_b(2:end,1)-roc_b(1:end-1,1)) .* ...
    (roc_b(1:end-1,2).^2 + roc_b(1:end-1,2).*(roc_b(2:end,2)-roc_b(1:end-1,2)) +...
    (roc_b(2:end,2)-roc_b(1:end-1,2)).*(roc_b(2:end,2)-roc_b(1:end-1,2))/3) );
q1 = sum( ...
    (roc_b(2:end,2)-roc_b(1:end-1,2)) .* ...
    ((1-roc_b(2:end,1)).^2 + (1-roc_b(2:end,1)).*(roc_b(2:end,1)-roc_b(1:end-1,1)) +...
    (roc_b(2:end,1)-roc_b(1:end-1,1)).*(roc_b(2:end,1)-roc_b(1:end-1,1))/3) );
auc_se = sqrt( (auc_b*(1-auc_b) + (n_pos-1)*(q1-auc_b^2) + (n_neg-1)*(q2-auc_b^2)) / (n_neg*n_pos) );
ci_95 = 1.96*auc_se;

%Plot curve
plot(roc_b(:,1), roc_b(:,2), '-', 'LineWidth', 2);
%Store legend labels
legend_text{4} = ['\fontsize{16}  Monogenic, {\bfA_z = ', num2str(auc_b, 3) '} \pm ' num2str(ci_95,2)];
%Compute sensitivity for some specificity and display result
sens = interp1(roc_b(:,1), roc_b(:,2), 0.05);
display(['Monogenic sensitivity at 95% specificity = ', num2str(sens)]);

%
% v) Linop, analytic method
load([results_dir 'linop_octave_05_08.mat'], '*roc_b', '*auc_b', 't_counts_b', 'f_counts_b');

%Compute errors on ROC curve
n_pos = sum(t_counts_b);
n_neg = sum(f_counts_b);
q2 = sum( ...
    (roc_b(2:end,1)-roc_b(1:end-1,1)) .* ...
    (roc_b(1:end-1,2).^2 + roc_b(1:end-1,2).*(roc_b(2:end,2)-roc_b(1:end-1,2)) +...
    (roc_b(2:end,2)-roc_b(1:end-1,2)).*(roc_b(2:end,2)-roc_b(1:end-1,2))/3) );
q1 = sum( ...
    (roc_b(2:end,2)-roc_b(1:end-1,2)) .* ...
    ((1-roc_b(2:end,1)).^2 + (1-roc_b(2:end,1)).*(roc_b(2:end,1)-roc_b(1:end-1,1)) +...
    (roc_b(2:end,1)-roc_b(1:end-1,1)).*(roc_b(2:end,1)-roc_b(1:end-1,1))/3) );
auc_se = sqrt( (auc_b*(1-auc_b) + (n_pos-1)*(q1-auc_b^2) + (n_neg-1)*(q2-auc_b^2)) / (n_neg*n_pos) );
ci_95 = 1.96*auc_se;

%Plot curve
plot(roc_b(:,1), roc_b(:,2), '-', 'LineWidth', 2);
%Store legend labels
legend_text{5} = ['\fontsize{16}  Linop, {\bfA_z = ', num2str(auc_b, 3) '} \pm ' num2str(ci_95,2)];

%Compute sensitivity for some specificity and display result
sens = interp1(roc_b(1:80,1), roc_b(1:80,2), 0.05);
display(['Linop sensitivity at 95% specificity = ', num2str(sens)]);

%
% Label title and axes of graph
title('\fontsize{18} \bf ROC curves for line detection')
xlabel('\fontsize{16} \bf False positive rate');
ylabel('\fontsize{16} \bf True positive rate');
legend(legend_text, 'location', 'southeast');

%Save graph
print('-dtiff', '-noui', '-painters', f, '-r300', 'M:\asymmetry_project\Paper submissions\bmvc2010\figures\line_detection_roc.tif');
clear f;
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%2. Generate ROC curves for spicule classification
decomp = {'DT', 'Monogenic', 'Linop'};
decomp2 = {'DT-CWT', 'Monogenic', 'Linop'};
for win_size = [1 3] %loop through different window sizes
    for level = 5:6 %loop through different levels
        f = figure(...
            'windowstyle', 'normal',...
            'Units', 'pixels',...
            'position', [50 50 800 800],...
            'PaperPositionMode','auto'); 
        hold all; axis equal; axis([0 1 0 1]);
        leg_text = cell(1,1);

        j = 1;
        for dd = 1:3 %loop through different decomposition types
        
        
            pooled_test_labels = [];
            pooled_test_scores = [];
            fold_aucs = zeros(10,1);

            for ii = 1:10 %loop through different folds of cross-fold evaluation
                
                %load in data
                load(...
                    ['M:\chen\data\spicules\' decomp{dd} ...
                     '\rf_spic_' decomp{dd} '_' zerostr(win_size,1) '_' zerostr(level,1) '_all_', zerostr(ii,2) '_results.mat']);
                fold_accuracy = mean(test_labels == str2double(test_predictions));
                display(['accuracy of fold ', num2str(ii), ' = ', num2str(fold_accuracy)]);

                test_scores = test_votes(:,2) / sum(test_votes(1,:));

                %Pool labels and scores for this fold
                pooled_test_labels = [pooled_test_labels; test_labels]; %#ok
                pooled_test_scores = [pooled_test_scores; test_scores]; %#ok

            end 
            
            % Comptue ROC statistics for pooled data
            [pooled_roc_pts pooled_auc dum dum auc_se] = calculate_roc_curve(pooled_test_scores,pooled_test_labels,linspace(-0.0001,1.0001,100));
            ci_95 = 1.96*auc_se;
            
            %Plot curve
            plot(pooled_roc_pts(:,1), pooled_roc_pts(:,2), '-', 'linewidth', 2);

            %Save legend label
            leg_text{j} = ['\fontsize{16} ' decomp2{dd} ', {\bfA_{z}= ' num2str(pooled_auc,3) '} \pm ' num2str(ci_95,2)];
            j = j+1;
        end
        
        %Label title and axes
        xlabel('\fontsize{16} \bf False positive rate');
        ylabel('\fontsize{16} \bf True positive rate');
        title({...
            '\fontsize{18} \bf ROC curves for pooled cross-fold classification - ';
            ['spicule vs non-spicule, {\fontname{times} \it w} = ' num2str(win_size) ', {\fontname{times} \it s} = ' num2str(level)]});
        legend(leg_text, 'location', 'southeast');
        
        %Save the figure
        f_name = ['M:\asymmetry_project\Paper submissions\bmvc2010\figures\rf_spic_' zerostr(win_size,1) '_' zerostr(level,1) '.tif'];
        print('-dtiff', '-noui', '-painters', f, '-r300', f_name);
        clear f;
    end  
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% 3. Generate monogenic and linop line probabilities for a real mass - Chen
% has alreday created maps for all the random forest methods
load('M:\chen\data\line_detection_mammo\mass028.mat')
figure; imagesc(mass_roi); axis image; colormap(gray(256));
write_im_from_colormap(mass_roi, 'M:\asymmetry_project\Paper submissions\bmvc2010\figures\mass028.bmp', gray(256));
prob_im_linop = line_operator_octave(mass_roi, 8, 5, 'degrees', 0);
prob_im_mono = monogenic_phase_cong(mass_roi, 3, 4, 2, 0.65);
write_im_from_colormap(prob_im_linop, 'M:\asymmetry_project\Paper submissions\bmvc2010\figures\mass028_linop.bmp', gray(256));
write_im_from_colormap(prob_im_mono, 'M:\asymmetry_project\Paper submissions\bmvc2010\figures\mass028_monogenic.bmp', gray(256));
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% 4. Save probability maps (as bmp images) for a real test image
load('M:\chen\data\testimage_contrast1to8_exprnd_sin\image005.mat');      
prob_im_dt_3_5 = u_load('M:\chen\data\testimage_contrast1to8_exprnd_sin\probability_image\182392\probability_image005.mat');
prob_im_dt_3_mag = u_load('M:\chen\data\testimage_contrast1to8_exprnd_sin\probability_image\182418\probability_image005.mat');
prob_im_dt_3_pha = u_load('M:\chen\data\testimage_contrast1to8_exprnd_sin\probability_image\182419\probability_image005.mat');
prob_im_linop_rf = u_load('M:\chen\data\testimage_contrast1to8_exprnd_sin\probability_image\rf_linop_1_5_8\probability_image005.mat');
prob_im_linop = u_load('M:\chen\data\testimage_contrast1to8_exprnd_sin\probability_image\linop_octave_05_08\probability_image005.mat');
prob_im_monogenic = u_load('M:\chen\data\testimage_contrast1to8_exprnd_sin\probability_image\monogenic_0p65_04_03\probability_image005.mat');
prob_im_monogenic_rf = u_load('M:\chen\data\testimage_contrast1to8_exprnd_sin\probability_image\monogenic_rf_trees_W3L4\probability_image005.mat');

figure; imagesc(prob_im_dt_3_5); axis image; colormap(gray(256));
figure; imagesc(prob_im_dt_3_mag); axis image; colormap(gray(256));
figure; imagesc(prob_im_dt_3_pha); axis image; colormap(gray(256));
figure; imagesc(prob_im_linop_rf); axis image; colormap(gray(256));
figure; imagesc(prob_im_linop); axis image; colormap(gray(256));
figure; imagesc(prob_im_monogenic); axis image; colormap(gray(256));
figure; imagesc(prob_im_monogenic_rf); axis image; colormap(gray(256));

% Have to create a new figure for using only max of DT as don't yet have
% results for new test images (I know this model is trained on the old data
% - but it doesn't matter just for display purposes)
load('\\isbe-matrix\mammography\chen\data\rf_trees_182267\tree_combine\random_forest.mat')
random_forest.tree_dir = 'M:\chen\data\rf_trees_182267\tree_combine\';
prob_im_dt_max = classify_image('image_in', image, 'forest', random_forest, 'do_max', 1);
%
write_im_from_colormap(image, 'M:\asymmetry_project\Paper submissions\bmvc2010\figures\lines_test.bmp', gray(256));
write_im_from_colormap(prob_im_dt_3_5, 'M:\asymmetry_project\Paper submissions\bmvc2010\figures\lines_dt_pm.bmp', gray(256));
write_im_from_colormap(prob_im_dt_3_mag, 'M:\asymmetry_project\Paper submissions\bmvc2010\figures\lines_dt_mm.bmp', gray(256));
write_im_from_colormap(prob_im_dt_max, 'M:\asymmetry_project\Paper submissions\bmvc2010\figures\lines_dt_max.bmp', gray(256));
write_im_from_colormap(prob_im_dt_3_pha, 'M:\asymmetry_project\Paper submissions\bmvc2010\figures\lines_dt_pp.bmp', gray(256));
write_im_from_colormap(prob_im_linop_rf, 'M:\asymmetry_project\Paper submissions\bmvc2010\figures\lines_linop.bmp', gray(256));
write_im_from_colormap(prob_im_linop, 'M:\asymmetry_project\Paper submissions\bmvc2010\figures\lines_linop_rf.bmp', gray(256));
write_im_from_colormap(prob_im_monogenic, 'M:\asymmetry_project\Paper submissions\bmvc2010\figures\lines_monogenic.bmp', gray(256));
write_im_from_colormap(prob_im_monogenic_rf, 'M:\asymmetry_project\Paper submissions\bmvc2010\figures\lines_monogenic_rf.bmp', gray(256));