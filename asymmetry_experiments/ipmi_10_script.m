%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%----------------------- IPMI paper script --------------------------------
%--------------------------------------------------------------------------

%Script to generate figures for IPMI 2011 paper
%--------------------------------------------------------------------------

test_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\';
prob_dir = 'Z:\data\synthetic_lines\lines512\results\';
results_dir = 'C:\isbe\asymmetry_project\experiments\line_detection\';
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%COmpute line and orientation maps for analytic g2d methods
param_dir = 'g2d_scales_16';

mkdir([prob_dir param_dir '\lines\']);
mkdir([prob_dir param_dir '\orientations\']);

for ii = 11:100
    display(['Processing image ', num2str(ii) ' of 100']);
    %Load in test image
    test_im = u_load([test_dir '\image' zerostr(ii,3) '.mat']);
    
    %
    [line_strength, orientation_map, scale_map] = gaussian_2nd_derivative_line(test_im, [1 2 4 8]);
    [grad_strength, grad_orientation] = gaussian_1st_derivative_gradient(test_im, 10);

    discard_map = ...
        ((abs(mb_mod(orientation_map - grad_orientation + pi/2, pi)) < pi/6) &...
        (grad_strength > 25)) | ...
        (line_strength > 0);
    
    line_map = abs(line_strength);
    line_map(discard_map) = 0;
    line_map = line_map ./ max(line_map(:));
    
    if ii < 11
        figure; colormap(gray(256));
        subplot(1,2,1); imagesc(test_im); axis image;
        subplot(1,2,2); imagesc(line_map); axis image;
    end
        
    save([prob_dir param_dir '\lines\image' zerostr(ii,3) '_line.mat'], 'line_map');
    save([prob_dir param_dir '\orientations\image' zerostr(ii,3) '_ori.mat'], 'orientation_map');

end
%%
% Compute line and orientation maps for analytic linop
param_dir = 'linop_odd_53_08';
mkdir([prob_dir param_dir '\lines\']);
mkdir([prob_dir param_dir '\orientations\']);

for ii = 1:100
    test_im = u_load([test_dir '\image' zerostr(ii,3) '.mat']);
    [line_map orientation_map] = line_operator_conv(test_im, 8, 3, 53, 'degrees');
    save([prob_dir param_dir '\lines\image' zerostr(ii,3) '_line.mat'], 'line_map');
    save([prob_dir param_dir '\orientations\image' zerostr(ii,3) '_ori.mat'], 'orientation_map');

end
%%
%Compute line detection and orientation regression errors for all the
%various methods

%Line detection RFs first
param_dir = {'191905', '191960', '191961', '233141', '192038'};

for ii = 4
    %Compute ROC scores
    [roc_b, auc_b, tp_counts_b, fp_counts_b, t_counts_b, f_counts_b] =...
        compute_roc_image_set_lines(test_dir, [prob_dir param_dir{ii} '\'], 'centre_line');
    %Save results
    save([results_dir param_dir{ii} '_roc_results.mat'], '*_b');
end

%%
param_dir = 'g2d_scales_16';
[roc_b, auc_b, tp_counts_b, fp_counts_b, t_counts_b, f_counts_b] =...
    compute_roc_image_set_lines(test_dir, [prob_dir param_dir '\lines\'], 'centre_line');
save([results_dir param_dir '_roc_results.mat'], '*_b');
%%
%Now do the same for orientations
param_dir = {'213209', '191958', '191959', '233142'};
for ii = 4
    %Compute orientation errors
    [orientation_errors] = compute_image_orientation_errors(...
        test_dir, [prob_dir param_dir{ii}], 'centre_line');
    
    %Save results
    save([results_dir param_dir{ii} '_ori_results.mat'], 'orientation_errors');
end

param_dir = 'g2d_scales_16';
[orientation_errors] = compute_image_orientation_errors(...
     test_dir, [prob_dir param_dir '\orientations\'], 'centre_line');
save([results_dir param_dir '_ori_results.mat'], 'orientation_errors');

%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
map_types = {'f1', 'f2', 'ni', 'np'};
for ii = 1:4
    [rf.thresh rf.percentiles] =...
        mass_maps_thresh('k_stellate_maps_rf', 'dist_ranges', [60 90 120 150 180], 'sigma', 8, 'map_type', map_types{ii});
    [rf.hits rf.num_maxima rf.maxima_ranks rf.maxima_scores] =...
        mass_maps_maxima('k_stellate_maps_rf', rf.thresh, 'dist_ranges', [60 90 120 150 180], 'sigma', 8, 'map_type', map_types{ii});
    save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rf_hits_' map_types{ii} '.mat'], 'rf');
    %
    [wrf.thresh wrf.percentiles] =...
        mass_maps_thresh('k_stellate_maps_wrf', 'dist_ranges', [60 90 120 150 180], 'sigma', 8, 'map_type', map_types{ii});
    [wrf.hits wrf.num_maxima wrf.maxima_ranks wrf.maxima_scores] =...
        mass_maps_maxima('k_stellate_maps_wrf', wrf.thresh, 'dist_ranges', [60 90 120 150 180], 'sigma', 8, 'map_type', map_types{ii});
    save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrf_hits_' map_types{ii} '.mat'], 'wrf');
    %
    [g2.thresh g2.percentiles] =...
        mass_maps_thresh('k_stellate_maps_g2', 'dist_ranges', [60 90 120 150 180], 'sigma', 8, 'map_type', map_types{ii}, 'data_type', '2004_screening_processed/abnormals');
    [g2.hits g2.num_maxima g2.maxima_ranks g2.maxima_scores] =...
        mass_maps_maxima('k_stellate_maps_g2', g2.thresh, 'dist_ranges', [60 90 120 150 180], 'sigma', 8, 'map_type', map_types{ii}, 'data_type', '2004_screening_processed/abnormals');
    save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2_hits_' map_types{ii} '.mat'], 'g2');
end
%%
map_types = {'f1', 'f2', 'ni', 'np'};
for ii = 1:4
    
    load(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rf_hits_' map_types{ii} '.mat'], 'rf');
    display([map_types{ii} ' hits for rf:']);
    display(num2str(nansum(rf.hits)));
    %
    load(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrf_hits_' map_types{ii} '.mat'], 'wrf');
    display([map_types{ii} ' hits for wrf:']);
    display(num2str(nansum(wrf.hits)));
    %
    load(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2_hits_' map_types{ii} '.mat'], 'g2');
    display([map_types{ii} ' hits for g2:']);
    display(num2str(nansum(g2.hits)));
    
end
%%
map_types = {'f1', 'f2'};
for ii = 2
    [rf.thresh rf.percentiles] =...
        mass_maps_thresh('k_stellate_maps_rf_scaled', 'sigma', 8, 'map_type', map_types{ii});
    [rf.hits rf.num_maxima rf.maxima_ranks rf.maxima_scores] =...
        mass_maps_maxima('k_stellate_maps_rf_scaled', rf.thresh, 'sigma', 8, 'map_type', map_types{ii});
    save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rfs_hits_' map_types{ii} '.mat'], 'rf');
    %
    [wrf.thresh wrf.percentiles] =...
        mass_maps_thresh('k_stellate_maps_wrf_scaled', 'sigma', 8, 'map_type', map_types{ii});
    [wrf.hits wrf.num_maxima wrf.maxima_ranks wrf.maxima_scores] =...
        mass_maps_maxima('k_stellate_maps_wrf_scaled', wrf.thresh, 'sigma', 8, 'map_type', map_types{ii});
    save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrfs_hits_' map_types{ii} '.mat'], 'wrf');
    %
    [g2.thresh g2.percentiles] =...
        mass_maps_thresh('k_stellate_maps_g2_scaled', 'sigma', 8, 'map_type', map_types{ii}, 'data_type', '2004_screening_processed/abnormals');
    [g2.hits g2.num_maxima g2.maxima_ranks g2.maxima_scores] =...
        mass_maps_maxima('k_stellate_maps_g2_scaled', g2.thresh, 'sigma', 8, 'map_type', map_types{ii}, 'data_type', '2004_screening_processed/abnormals');
    save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2s_hits_' map_types{ii} '.mat'], 'g2');
end
%%
map_types = {'f1', 'f2'};
mam_type = {'normals', 'abnormals'};
for ii = 1:2
    for jj = 2
        [rf.tp rf.fp rf.fp_pixels] = mass_maps_roc(...
            'k_stellate_maps_rf_scaled', ['2004_screening/' mam_type{jj}], 'map_type', map_types{ii});

        save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\rfs_froc_'...
            mam_type{jj} '_' map_types{ii} '.mat'], 'rf');
        
        if ii == 2, continue, end
        
        [wrf.tp wrf.fp wrf.fp_pixels] = mass_maps_roc(...
            'k_stellate_maps_wrf_scaled', ['2004_screening/' mam_type{jj}], 'map_type', map_types{ii});

        save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\wrfs_froc_'...
            mam_type{jj} '_' map_types{ii} '.mat'], 'wrf');
        
        [g2.tp g2.fp g2.fp_pixels] = mass_maps_roc(...
            'k_stellate_maps_g2_scaled', ['2004_screening_processed/' mam_type{jj}], 'map_type', map_types{ii});

        save(['C:\isbe\asymmetry_project\experiments\karssemeijer_maps\g2s_froc_'...
            mam_type{jj} '_' map_types{ii} '.mat'], 'g2');
    end
end
%%
clear; pack;
for win_size = [1 3]
    for level = [5 4]  
        for ii = 1:10
            spicule_classification_crossfold(ii, 10, level, win_size, 'all', 'g2d');
            pack;
        end
    end
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Generate ROC curves for line detection methods
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%1. Create ROC curves of different line detection methods
results_dir = 'C:\isbe\asymmetry_project\experiments\line_detection\';

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
load([results_dir '191905_roc_results.mat'], '*roc_b', '*auc_b', 't_counts_b', 'f_counts_b');

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
legend_text{1} = ['\fontsize{20}  DT-CWT/RF, {\bfA_z = ', num2str(auc_b, 3) '} \pm ' num2str(ci_95,2)];
%Compute sensitivity for some specificity and display result
sens = interp1(roc_b(:,1), roc_b(:,2), 0.1);
display(['DT sensitivity at 95% specificity = ', num2str(sens)]);

%
% ii) Linop with random forest
% Load in saved ROC curves
load([results_dir '191961_roc_results.mat'], '*roc_b', '*auc_b', 't_counts_b', 'f_counts_b');

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
legend_text{2} = ['\fontsize{20}  Linop/RF, {\bfA_z = ', num2str(auc_b, 3) '} \pm ' num2str(ci_95,2)];

%Compute sensitivity for some specificity and display result
sens = interp1(roc_b(:,1), roc_b(:,2), 0.1);
display(['Linop/rf sensitivity at 95% specificity = ', num2str(sens)]);

%
% iii) Gaussian with random forest
% Load in saved ROC curves
load([results_dir '233141_roc_results.mat'], '*roc_b', '*auc_b', 't_counts_b', 'f_counts_b');

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
legend_text{3} = ['\fontsize{20}  Gaussian/RF, {\bfA_z = ', num2str(auc_b, 3) '} \pm ' num2str(ci_95,2)];

%Compute sensitivity for some specificity and display result
sens = interp1(roc_b(:,1), roc_b(:,2), 0.1);
display(['Gaussian/rf sensitivity at 95% specificity = ', num2str(sens)]);

%
% iv) Monogenic with random forest
load([results_dir '191960_roc_results.mat'], '*roc_b', '*auc_b', 't_counts_b', 'f_counts_b');

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
legend_text{4} = ['\fontsize{20}  Monogenic/RF, {\bfA_z = ', num2str(auc_b, 3) '} \pm ' num2str(ci_95,2)];

%Compute sensitivity for some specificity and display result
sens = interp1(roc_b(:,1), roc_b(:,2), 0.1);
display(['Monogenic/rf sensitivity at 95% specificity = ', num2str(sens)]);

%
% v) Monogenic, analytic method
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
legend_text{5} = ['\fontsize{20}  Monogenic, {\bfA_z = ', num2str(auc_b, 3) '} \pm ' num2str(ci_95,2)];
%Compute sensitivity for some specificity and display result
sens = interp1(roc_b(:,1), roc_b(:,2), 0.1);
display(['Monogenic sensitivity at 95% specificity = ', num2str(sens)]);

%
% vi) Gaussian, analytic method
load([results_dir 'g2d_scales_16_roc_results.mat'], '*roc_b', '*auc_b', 't_counts_b', 'f_counts_b');

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
legend_text{6} = ['\fontsize{20}  Gaussian, {\bfA_z = ', num2str(auc_b, 3) '} \pm ' num2str(ci_95,2)];

%Compute sensitivity for some specificity and display result
sens = interp1(roc_b(1:80,1), roc_b(1:80,2), 0.1);
display(['Linop sensitivity at 95% specificity = ', num2str(sens)]);

%
% vii) Linop, analytic method
load([results_dir 'linop_odd_53_08.mat'], '*roc_b', '*auc_b', 't_counts_b', 'f_counts_b');

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
legend_text{7} = ['\fontsize{20}  Linop, {\bfA_z = ', num2str(auc_b, 3) '} \pm ' num2str(ci_95,2)];

%Compute sensitivity for some specificity and display result
sens = interp1(roc_b(1:80,1), roc_b(1:80,2), 0.1);
display(['Linop sensitivity at 95% specificity = ', num2str(sens)]);

%
% Label title and axes of graph
title('\fontsize{24} \bf ROC curves for line detection')
xlabel('\fontsize{20} \bf False positive rate');
ylabel('\fontsize{20} \bf True positive rate');
legend(legend_text, 'location', 'southeast');

%Save graph
print('-dtiff', '-noui', '-painters', f, '-r300', 'M:\asymmetry_project\Paper submissions\ipmi2011\figures\line_detection_roc.tif');
clear f;
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%2. Generate ROC curves for spicule classification
decomp = {'DT', 'Monogenic', 'g2d', 'Linop'};
decomp2 = {'DT-CWT', 'Monogenic', 'Gaussian', 'Linop'};
for win_size = 1 %loop through different window sizes
    for level = 6 %loop through different levels
        f = figure(...
            'windowstyle', 'normal',...
            'Units', 'pixels',...
            'position', [50 50 800 800],...
            'PaperPositionMode','auto'); 
        hold all; axis equal; axis([0 1 0 1]);
        leg_text = cell(1,1);

        j = 1;
        for dd = 1:4 %loop through different decomposition types
        
        
            pooled_test_labels = [];
            pooled_test_scores = [];
            fold_aucs = zeros(10,1);

            for ii = 1:10 %loop through different folds of cross-fold evaluation
                
                %load in data
                load(...
                    ['I:\asymmetry_project\experiments\spicule_classification\' ...
                     decomp{dd} '\rf_spic_' decomp{dd} '_' zerostr(win_size,1) '_' zerostr(level,1) '_all_', zerostr(ii,2) '_results.mat']);
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
            leg_text{j} = ['\fontsize{20} ' decomp2{dd} ', {\bfA_{z}= ' num2str(pooled_auc,3) '} \pm ' num2str(ci_95,2)];
            j = j+1;
        end
        
        %Label title and axes
        xlabel('\fontsize{20} \bf False positive rate');
        ylabel('\fontsize{20} \bf True positive rate');
        title({...
            '\fontsize{24} \bf ROC curves for pooled cross-fold classification - ';
            ['spicule vs non-spicule, {\fontname{times} \it w} = ' num2str(win_size) ', {\fontname{times} \it s} = ' num2str(level)]});
        legend(leg_text, 'location', 'southeast');
        
        %Save the figure
        f_name = ['M:\asymmetry_project\Paper submissions\ipmi2011\figures\rf_spic_' zerostr(win_size,1) '_' zerostr(level,1) '.tif'];
        print('-dtiff', '-noui', '-painters', f, '-r300', f_name);
        clear f;
    end  
end
%%
%--------------------------------------------------------------------------
%Make images of the line detection results
line_im = u_load([test_dir 'image003.mat']);
write_im_from_colormap(line_im, 'M:\asymmetry_project\Paper submissions\ipmi2011\figures\line512_003.bmp');

param_dir = {'191905', '191960', '191961', '233141'};

for ii = 1:4
    im_list = dir([prob_dir param_dir{ii} '\*.mat']);
    line_map = load_uint8([prob_dir param_dir{ii} '\' im_list(3).name]);
    figure; imagesc(line_map); axis image; colormap(gray(256));
    write_im_from_colormap(line_map, ['M:\asymmetry_project\Paper submissions\ipmi2011\figures\line512_003_rf_' param_dir{ii} '.bmp'], gray(256));
end
%%
line_map = u_load('C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\g2d_scales_16\lines\image003_line.mat');
write_im_from_colormap(line_map, 'M:\asymmetry_project\Paper submissions\ipmi2011\figures\line512_003_gaussian.bmp', gray(256));
figure; imagesc(line_map); axis image; colormap(gray(256));

line_map = u_load('C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\monogenic_0p65_08_02\probability_image003.mat');
write_im_from_colormap(line_map, 'M:\asymmetry_project\Paper submissions\ipmi2011\figures\line512_003_monogenic.bmp', gray(256));
figure; imagesc(line_map); axis image; colormap(gray(256));

line_map = u_load('C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\linop_odd_53_08\lines\probability_image003.mat');
write_im_from_colormap(line_map, 'M:\asymmetry_project\Paper submissions\ipmi2011\figures\line512_003_lines.bmp', gray(256));
figure; imagesc(line_map); axis image; colormap(gray(256));
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Now do the same for orientations
warning('off', 'load_uint8:missing_variables');

wstyle = 'docked';
label_type = 'all_line';

f1 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [50 50 1000 800],...
    'PaperPositionMode','auto');
a1 = axes; hold all; 
title('\fontsize{24} \bf  CDF of orientation estimation errors'); 
xlabel('\fontsize{20} \bf Orientation error (degrees)');
% figure('windowstyle', wstyle); a2 = axes; hold all; 
% title('CDF of dispersion magnitudes');
% xlabel('Dispersion magnitudes');
% figure('windowstyle', wstyle); a3 = axes; hold all; 
% title('CDF of orientation errors weighted by dispersion magnitudes');
% xlabel('Weighted orientation errors');
% figure('windowstyle', wstyle); a4 = axes; hold all; 
% title('Mean orientation error for the Nth percentile of dispersion magnitudes');
% ylabel('Mean orientation error')
% xlabel('Percentile of samples sorted by dispersion magnitude');

ori_dir = {'213209', '191959', '233142', '191958', 'monogenic_0p65_08_02_ori', 'g2d_scales_16\orientations', 'linop_odd_53_08\orientations'};
line_dir = {'191905', '191961', '233141', '191960', 'monogenic_0p65_08_02', 'g2d_scales_16\lines', 'linop_odd_53_08\lines'};
test_label = {'DT-CWT/RF', 'Linop/RF', 'Gaussian/RF', 'Monogenic/RF', 'Monogenic', 'Gaussian', 'Linop'};

ori_leg = cell(length(ori_dir)+1,1);
ori_leg{1} = '{\fontsize{20} \bf Mean (Median) Angular Errors} \bf \color{white} z_z';
plot(a1, 0, 0, 'w');
% mag_leg = cell(length(ori_dir),1);
% com_leg = cell(length(ori_dir),1);
% pct_leg = cell(length(ori_dir),1);
%
for ii = 1:length(ori_dir)    
    %Save results

    errors = compute_image_orientation_errors(...
        test_dir, [prob_dir ori_dir{ii} '\'], label_type);
%     errors(:,2) = get_image_classifications(...
%         test_dir, [prob_dir line_dir{ii} '\'], label_type);
    
    ori_errors = sort(abs(errors(:,1)));
%     mag_errors = sort(errors(:,2));
%     com_errors = sort(prod(abs(errors),2) / mean(mag_errors));
%     errors = sortrows(errors, -2);
    
    ori_cdf = zeros(101,1);
%     mag_cdf = zeros(101,1);
%     com_cdf = zeros(101,1);
%     mean_pct = zeros(100,1);
    for jj = 1:100
        x = ceil(jj*size(ori_errors, 1)/100);
        ori_cdf(jj+1) = ori_errors(x,1);
%         mag_cdf(jj+1) = mag_errors(x,1);
%         com_cdf(jj+1) = com_errors(x,1);
%         mean_pct(jj) = mean(abs(errors(1:x,1)));
    end
    
    plot(a1, ori_cdf, (0:100)/100, 'linewidth', 2);
%     plot(a2, mag_cdf, (0:100)/100, 'linewidth', 2);
%     plot(a3, com_cdf, (0:100)/100, 'linewidth', 2);
%     plot(a4, 1:100, mean_pct, 'linewidth', 2);
    
    ori_leg{ii+1} = ...
        ['\fontsize{24} ' test_label{ii} ': ' sprintf('%3.1f', mean(ori_errors)) ' (' sprintf('%3.1f', median(ori_errors)) ') \bf \color{white} z_z'];
%     mag_leg{ii} = ...
%         ['\fontsize{20} ' test_label{ii} ': (mean, median) = (' num2str(mean(mag_errors),4) ', ' num2str(median(mag_errors),4) ')'];
%     com_leg{ii} = ...
%         ['\fontsize{20} ' test_label{ii} ': (mean, median) = (' num2str(mean(com_errors),4) ', ' num2str(median(com_errors),4) ')'];
%     pct_leg{ii} = ['\fontsize{20} ' test_label{ii} ': mean at 50th pcntile = ' num2str(mean_pct(50),4)];
end

legend(a1, ori_leg, 'location', 'southeast');
% legend(a2, mag_leg, 'location', 'southeast');
% legend(a3, com_leg, 'location', 'southeast');
% legend(a4, pct_leg, 'location', 'southeast');
%%
%Save graph
print('-dtiff', '-noui', '-painters', f1, '-r300', 'M:\asymmetry_project\Paper submissions\ipmi2011\figures\orientation_estimation_cdf.tif');
clear f;
%%
for ii =1:length(ori_dir)
    ori_leg{ii} = ['\fontsize{20} ' ori_leg{ii}];
    mag_leg{ii} = ['\fontsize{20} ' mag_leg{ii}];
    com_leg{ii} = ['\fontsize{20} ' com_leg{ii}];
    pct_leg{ii} = ['\fontsize{20} ' pct_leg{ii}];
end
legend(a1, ori_leg, 'location', 'southeast');
legend(a2, mag_leg, 'location', 'southeast');
legend(a3, com_leg, 'location', 'southeast');
legend(a4, pct_leg, 'location', 'southeast');

%%

f = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [50 50 800 800],...
    'PaperPositionMode','auto'); 

ori_map = load_uint8('C:\isbe\asymmetry_project\data\orientation_maps\rf_thin\2004_screening_processed\mass_roi\039RCC_roi.mat');
mass = load_uint8('C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\mass_roi\039RCC_roi.mat');

ori_small = ori_map(145:656, 145:656);
mass_small = mass(145:656, 145:656);

display_orientation(f, mass_small, ori_small, 4);
print('-dtiff', '-noui', '-painters', f, '-r300', 'M:\asymmetry_project\Paper submissions\ipmi2011\figures\mass_ori_quiver.tif');
clear f;
