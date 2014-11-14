rand('twister', 5489);
sampling_method_args.bg_dir = '\\isbe-san1\mberks\dev\classification\data\smooth512x512_patches\train\';
sampling_method_args.save_path = [];
sampling_method_args.num_samples = 1e4;
sampling_method_args.width_range = [0 8];
sampling_method_args.contrast_range = [0 8];
sampling_method_args.win_size = 3;
sampling_method_args.num_levels = 3;
sampling_method_args.min_wavelength = 4;
[X_test3 y_test3] = sample_training_data_monogenic(sampling_method_args);
%

rand('twister', 5489);
sampling_method_args.bg_dir = '\\isbe-san1\mberks\dev\classification\data\smooth512x512_patches\train\';
sampling_method_args.save_path = [];
sampling_method_args.num_samples = 1e4;
sampling_method_args.width_range = [0 8];
sampling_method_args.contrast_range = [0 8];
sampling_method_args.win_size = 1;
sampling_method_args.num_levels = 3;
sampling_method_args.min_wavelength = 4;
[X_test1 y_test1] = sample_training_data_monogenic(sampling_method_args);

%
X_test1a = X_test3(:,5:9:end);
display(max(abs(X_test1a(:) - X_test1(:))));
display(max(abs(y_test3(:) - y_test1(:))));
%%
for ii = 1:10
    test_image = load(['M:\chen\data\testimage_contrast1to8_exprnd_sin\image' zerostr(ii,3) '.mat']);
    write_im_from_colormap(test_image.image, ['M:\asymmetry_project\presentation_for_cjt\test_image' zerostr(ii,3) '_enhanced.jpg'], gray(256));
    imwrite(uint8(test_image.image), ['M:\asymmetry_project\presentation_for_cjt\test_image' zerostr(ii,3) '.jpg']);
end
%%
test_dir = 'M:\chen\data\testimage_contrast1to8_exprnd_sin\'; %ZC changes these
prob_dir = 'M:\chen\data\testimage_contrast1to8_exprnd_sin\probability_image\';
results_dir = 'C:\isbe\dev\classification\line_detection_results\';

figure(...
    'WindowStyle', 'normal', ...
    'units', 'pixels',...
    'position', [10 10 800 800]);
hold all; axis equal; axis([0 1 0 1]);
legend_text = cell(5,1);

%
load([results_dir 'dt_3_5_all_pm.mat'], '*roc_a', '*auc_a');
plot(roc_a(:,1), roc_a(:,2), '-x', 'LineWidth', 2); legend_text{1} = ['DT-CWT/RF, 3x3, 5 levels, phase/mag, A_z = ', num2str(auc_a, 3)];
%
load([results_dir 'rf_linop_1_5_8.mat'], '*roc_a', '*auc_a');
plot(roc_a(:,1), roc_a(:,2), '-x', 'LineWidth', 2); legend_text{2} = ['Linop/RF, 1x1, 5 levels, phase/mag, A_z = ', num2str(auc_a, 3)];
%
load([results_dir 'monogenic_rf_trees_W3L5.mat'], '*roc_a', '*auc_a');
plot(roc_a(:,1), roc_a(:,2), '-x', 'LineWidth', 2); legend_text{3} = ['Monogenic/RF, 3x3, 5 levels, phase/mag, A_z = ', num2str(auc_a, 3)];
%
load([results_dir 'linop_octave_05_08.mat'], '*roc_a', '*auc_a');
plot(roc_a(:,1), roc_a(:,2), '-x', 'LineWidth', 2); legend_text{4} = ['Linop 5 levels, 8 angles, A_z = ', num2str(auc_a, 3)];
%
load([results_dir 'monogenic_0p65_08_02.mat'], '*roc_a', '*auc_a');
plot(roc_a(:,1), roc_a(:,2), '-x', 'LineWidth', 2); legend_text{5} = ['Monogenic, 2 levels, min wave length 8, real, A_z = ', num2str(auc_a, 3)];
%

title('Comparison of ROC curves for line detectors')
legend(legend_text, 'location', 'southeast');
%%
ex = linspace(0,pi, 100);
figure(...
    'WindowStyle', 'normal', ...
    'units', 'pixels',...
    'position', [10 10 800 800]);
hold all; axis([0 pi 0 1]);
legend_text = cell(11,1);

for squash = 0:10
    base_wave = sin(ex);
    for ii = 3:2:630; 
        base_wave = base_wave + squash*sin(ii*ex)/(ii*10); 
    end
    base_wave = base_wave / max(base_wave(:));
    
    plot(ex, base_wave, 'LineWidth', 2);
    legend_text{squash+1} = ['Box factor = ', num2str(squash/10,2)];
end
legend(legend_text, 'location', 'south');