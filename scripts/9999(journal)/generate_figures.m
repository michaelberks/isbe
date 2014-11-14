%--------------------------------------------------------------------------
% Figures
%--------------------------------------------------------------------------
fig_dir = 'C:\isbe\matlab_code\trunk\papers\9999 journal (orientation)\figs\';
%%
mkdir([fig_dir 'fibre']);
%%
%--------------------------------------------------------------------------
% Response across orientation for each filter type
g2d_responses_a = zeros(180,6);
gabor_responses_a = zeros(180,6);
dt_responses_a = zeros(180,6);
for i_ori = 1:180
    
    im = create_gauss_bar(1,1,i_ori,32,32,16,16);
    
    [g2d_responses] = compute_gaussian_2nd_derivatives(im, 2);
    g2d_responses_a(i_ori,:) = ...
        squeeze( abs(steer_gaussian_2nd_derivatives(g2d_responses(16,16,1,:), 6)) )';
    
    [gabor_responses] = compute_gabor_responses(im, 2, 6);
    gabor_responses_a(i_ori,:) = ...
        squeeze( real(gabor_responses(16,16,1,:)) )';
    
    dt = compute_dual_tree(im, 2, 0);
    dt_responses_a(i_ori,:) = ...
        squeeze( abs(sample_dt_data(dt, 16, 16, 'levels', 2, 'win_size', 1)) )';
end

g2d_responses_a = g2d_responses_a / max(g2d_responses_a(:));
gabor_responses_a = gabor_responses_a / max(gabor_responses_a(:));
dt_responses_a = dt_responses_a / max(dt_responses_a(:));

figure; 
plot(g2d_responses_a, 'linewidth', 3);
title('Response for 6 steered orientations of Gaussian 2nd derivative filters',...
    'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
ylabel('Absolute response',...
    'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
xlabel('Line orientation (degrees)',...
    'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
exportfig([fig_dir 'filtering\g2d_oriented_response.png']);

figure; 
plot(gabor_responses_a, 'linewidth', 3);
title('Response for 6 oriented of Gabor filters',...
    'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
ylabel('Absolute response',...
    'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
xlabel('Line orientation (degrees)',...
    'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
exportfig([fig_dir 'filtering\gabor_oriented_response.png']);

figure; 
plot(dt_responses_a, 'linewidth', 3);
title('Response for the 6 oriented sub-bands in the DT-CWT',...
    'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
ylabel('Absolute response',...
    'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
xlabel('Line orientation (degrees)',...
    'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
exportfig([fig_dir 'filtering\dt_oriented_response.png']);
%%
% Response across orientation for each filter type
g2d_responses_a = zeros(32,4);
gabor_responses_a = zeros(32,4);
dt_responses_a = zeros(32,4);
for i_width = 1:32
    
    im = create_gauss_bar(i_width/2,1,90,32,32,16,16);
    
    [g2d_responses] = compute_gaussian_2nd_derivatives(im, [1 2 4 8]);
    g2d_responses_a(i_width,:) = ...
        squeeze( abs( steer_gaussian_2nd_derivatives(g2d_responses(16,16,:,:), 1)) )';

    [gabor_responses] = compute_gabor_responses(im, [1 2 4 8], 1);
    gabor_responses_a(i_width,:) = ...
        squeeze( real(gabor_responses(16,16,:,1)) )';
    
    dt = compute_dual_tree(im, 4, 0);
    temp = squeeze( abs(sample_dt_data(dt, 16, 16, 'levels', 1:4, 'win_size', 1)) );
    dt_responses_a(i_width,:) = temp(3,:);
        
end
figure; 
plot(g2d_responses_a);
title('Response for increasing Gaussian 2nd derivative filters',...
    'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
ylabel('Absolute response');
xlabel('Line orientation (degrees)');
legend({'\sigma = 1', '\sigma = 2', '\sigma = 4', '\sigma = 8'}, 'location', 'se');
exportfig([fig_dir 'filtering\g2d_scale_response.png']);

figure; 
plot(gabor_responses_a);
title('Response for increasing scale Gabor filters');
ylabel('Absolute response');
xlabel('Line orientation (degrees)');
legend({'\sigma = 1', '\sigma = 2', '\sigma = 4', '\sigma = 8'}, 'location', 'se');
exportfig([fig_dir 'filtering\gabor_scale_response.png']);

figure; 
plot(dt_responses_a);
title('Response at each level in the DT-CWT');
ylabel('Absolute response');
xlabel('Line orientation (degrees)');
legend({'Level = 1', 'Level = 2', 'Level = 4', 'Level = 8'}, 'location', 'se');
exportfig([fig_dir 'filtering\dt_scale_response.png']);
%%
%3D version - this will be slow to compute, so if possible, load in
%precomputed values. Also, axes labels and titles look out of place on 3d
%plots so they're omitted
load('C:\isbe\asymmetry_project\experiments\filter_responses_surfaces.mat');
oris = 4:4:180;
widths = (1:32)/2;

figure; a1 = gca; hold all;
figure; a2 = gca; hold all;
figure; a3 = gca; hold all;

for i_scale = 2:5
    for i_band = 1:6    
        axes(a1); surface(widths, oris, g2d_responses_a(:,:,i_band,i_scale), 'edgecolor','none','facecolor', 'interp');%'edgecolor', colors(i_band), 
        axes(a2); surface(widths, oris, gabor_responses_a(:,:,i_band,i_scale), 'edgecolor','none','facecolor', 'interp');
        axes(a3); surface(widths, oris, dt_responses_a(:,:,i_band,i_scale), 'edgecolor','none','facecolor', 'interp');
    end
end
set(a1, 'view', [-29 52], 'zlim', [0 0.5]);
set(a2, 'view', [-29 52], 'zlim', [0 0.7]);
set(a3, 'view', [-29 52], 'zlim', [0 1.2]);

axes(a1);
exportfig([fig_dir 'filtering\g2d_3d_response.png']);
axes(a2); 
exportfig([fig_dir 'filtering\gabor_3d_response.png']);
axes(a3); 
exportfig([fig_dir 'filtering\dt_3d_response.png']);

figure;
for i_scale = 2:5
    surface(widths, oris, mono_responses_a(:,:,i_scale), 'edgecolor','none','facecolor', 'interp');
end
set(gca, 'view', [-29 52], 'zlim', [0 0.6]);
exportfig([fig_dir 'filtering\mono_3d_response.png']);
%--------------------------------------------------------------------------
%%
% Percentage of oriented bands correctly predicted: retinograms: Gabor
retroot = [asymmetryroot 'data\retinograms\DRIVE_clean\training\'];
load([asymmetryroot '/experiments/DRIVE/gabor_orientation_prediction/v_responses.mat'], 'v_responses');
load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations');

sample_idx = randperm(length(gt_orientations));
sample_idx = sample_idx(1:2e4);

v_responses = reshape(v_responses(sample_idx,:), [], 1, 5, 6);
gt_orientations = gt_orientations(sample_idx,:);

true_oris = angle(gt_orientations) / 2;
clear gt_orientations;
pack;

[~, predicted_bands, predicted_scales] = max_response_line(v_responses);
clear v_responses;

%Orientation
num_bins = 50;
ori_bins = linspace(-pi/2, pi/2, num_bins);
counts = zeros(6, num_bins);
for i_band = 1:6
    counts(i_band,:) = hist(true_oris(predicted_bands == i_band), ori_bins);
end

figure;
plot(ori_bins, bsxfun(@rdivide, counts, sum(counts))', 'linewidth', 3);
set(gca, 'ylim', [0 1], 'xlim', [-pi/2 pi/2], 'fontsize', 18);
xlabel('Orientation of vessel - $$\phi(p)$$',...
    'fontsize', 18, 'Interpreter', 'Latex');
ylabel('Percentage of pixels where $$\hat{\theta}_{Re}(p) = \theta$$',...
    'fontsize', 18, 'fontsize', 18, 'Interpreter', 'Latex');
title({'Which oriented sub-band produces maximum response?'; 'Retinogram DRIVE data'},...
    'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
legend({...
    '$$\theta = 0$$',...
    '$$\theta = ^{\pi}/_{6}$$',...
    '$$\theta = ^{2\pi}/_{6}$$',...
    '$$\theta = ^{3\pi}/_{6}$$',...
    '$$\theta = ^{4\pi}/_{6}$$',...
    '$$\theta = ^{5\pi}/_{6}$$'}, 'fontsize', 18, 'Interpreter', 'Latex');
exportfig([fig_dir 'retina\ret_vessels_gabor_ori_subbands.png']);

%Scale
load([retroot '/width_maps/gt/all_gt_widths.mat'], 'all_gt_widths');
all_gt_widths = all_gt_widths(sample_idx,:);
counts = zeros(5, 16);
for i_scale = 1:5
    counts(i_scale,:) = hist(all_gt_widths(predicted_scales == i_scale), 1:16);
end
figure; plot(1:16, bsxfun(@rdivide, counts, sum(counts))', 'linewidth', 3); 
xlabel('Width of vessel - $$w(p)$$',...
    'fontsize', 18, 'Interpreter', 'Latex');
ylabel('Percentage of pixels where $$\hat{\sigma}_{G}(p) = \sigma$$',...
    'fontsize', 18, 'Interpreter', 'Latex');
title({'Which scale produces maximum response?'; 'Retinogram DRIVE data'},...
    'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
legend({'$$\sigma = 1$$', '$$\sigma = 2$$', '$$\sigma = 4$$', '$$\sigma = 8$$', '$$\sigma = 16$$'},...
    'fontsize', 18, 'Interpreter', 'Latex');
exportfig([fig_dir 'retina\ret_vessels_gabor_scales.png']);
%%
%--------------------------------------------------------------------------
% Percentage of oriented bands correctly predicted: retinograms: DT
retroot = [asymmetryroot 'data\retinograms\DRIVE_clean\training\'];
load([asymmetryroot '\experiments\DRIVE\dt_orientation_prediction\v_responses.mat'], 'v_responses');
load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations');

sample_idx = randperm(length(gt_orientations));
sample_idx = sample_idx(1:2e4);

v_responses = reshape(v_responses(sample_idx,:), [], 1, 6, 5);
v_responses = permute(v_responses, [1 2 4 3]);
gt_orientations = gt_orientations(sample_idx,:);
true_oris = angle(gt_orientations) / 2;
clear gt_orientations;
pack;

[~, predicted_bands, predicted_scales] = max_response_line(v_responses);
clear v_responses;

num_bins = 50;
ori_bins = linspace(-pi/2, pi/2, num_bins);
counts = zeros(6, num_bins);
for i_band = 1:6
    counts(i_band,:) = hist(true_oris(predicted_bands == i_band), ori_bins);
end

figure; plot(ori_bins, bsxfun(@rdivide, counts, sum(counts))');
xlabel('Orientation of vessel - $$\phi(p)$$',...
    'fontsize', 18, 'Interpreter', 'Latex');
ylabel('Percentage of pixels where $$\hat{\theta}_{Re}(p) = \theta$$',...
    'fontsize', 18, 'fontsize', 18, 'Interpreter', 'Latex');
title({'Which oriented sub-band produces maximum response?'; 'Retinogram DRIVE data'},...
    'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
legend({...
    '$$\theta = 0$$',...
    '$$\theta = ^{\pi}/_{6}$$',...
    '$$\theta = ^{2\pi}/_{6}$$',...
    '$$\theta = ^{3\pi}/_{6}$$',...
    '$$\theta = ^{4\pi}/_{6}$$',...
    '$$\theta = ^{5\pi}/_{6}$$'}, 'fontsize', 18, 'Interpreter', 'Latex');
exportfig([fig_dir 'retina\ret_vessels_dt_ori_subbands.png']);
%
load([retroot '/width_maps/gt/all_gt_widths.mat'], 'all_gt_widths');
counts = zeros(5, 16);
width_bins = linspace(min(all_gt_widths), max(all_gt_widths), num_bins);
for i_scale = 1:5
    counts(i_scale,:) = hist(all_gt_widths(predicted_scales == i_scale), 1:16);
end

figure; plot(1:16, bsxfun(@rdivide, counts, sum(counts))'); 
xlabel('Orientation of vessel');
ylabel('Percentage of pixels with maximum response');
title('Which filter scale produces maximum response?');
legend({'Level = 1', 'Level = 2', 'Level = 3', 'Level = 4', 'Level = 5'});
exportfig([fig_dir 'retina\ret_vessels_dt_scales.png']);
%--------------------------------------------------------------------------
%%
% Percentage of oriented bands correctly predicted: Retinograms (DRIVE): G2d
retroot = [asymmetryroot 'data\retinograms\DRIVE_clean\training\'];
load([asymmetryroot '\experiments\DRIVE\g2d_orientation_prediction\v_responses.mat'], 'v_responses');

sample_idx = randperm(size(v_responses, 1));
sample_idx = sample_idx(1:2e4);

v_responses = reshape(v_responses(sample_idx,:), [], 1, 5, 3);
[~, predicted_theta, predicted_scales] = gaussian_2nd_derivative_line(v_responses);
clear v_responses;
%
load([retroot '/width_maps/gt/all_gt_widths.mat'], 'all_gt_widths');
counts = zeros(5, 16);
for i_scale = 1:5
    counts(i_scale,:) = hist(all_gt_widths(predicted_scales == i_scale), 1:16);
end

figure; plot(1:16, bsxfun(@rdivide, counts, sum(counts))', 'linewidth', 3); 
xlabel('Width of vessel - $$w(p)$$',...
    'fontsize', 18, 'Interpreter', 'Latex');
ylabel('Percentage of pixels where $$\hat{\sigma}_{G}(p) = \sigma$$',...
    'fontsize', 18, 'Interpreter', 'Latex');
title({'Which scale produces maximum response?'; 'Retinogram DRIVE data'},...
    'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
legend({'$$\sigma = 1$$', '$$\sigma = 2$$', '$$\sigma = 4$$', '$$\sigma = 8$$', '$$\sigma = 16$$'},...
    'fontsize', 18, 'Interpreter', 'Latex');
exportfig([fig_dir 'retina\ret_vessels_g2d_scales.png']);
%--------------------------------------------------------------------------
%%
% Percentage of scale bands correctly predicted: synthetic lines with
% increasing noise: G2d
exp_dir = 'C:\isbe\asymmetry_project\experiments\synthetic_lines\line_on_edge\';
for noise = [0 1 2 3]
    data_dir = [exp_dir 'test\rician_' num2str(noise) '\1\'];

    %Load and reformat data
    responses_g2d = u_load([data_dir 'responses_g2d.mat']);
    responses_g2d = reshape(responses_g2d, [], 9, 4, 3);
    responses_g2d = responses_g2d(:, 5, :, :);
    
    %Load ground truth data
    load([data_dir 'true_labels.mat'], 'true_oris', 'true_widths');
    
    [~, predicted_theta, predicted_scales] = gaussian_2nd_derivative_line(responses_g2d);
    clear responses_g2d;
    
    counts = zeros(4, 16);
    width_bins = linspace(1, 8, 16);
    for i_scale = 1:4
        counts(i_scale,:) = hist(true_widths(predicted_scales == 2^(i_scale-1)), width_bins);
    end

    figure; plot(width_bins, bsxfun(@rdivide, counts, sum(counts))', 'linewidth', 2); 
    xlabel('Width of vessel - $$w(p)$$', 'fontsize', 18, 'Interpreter', 'Latex');
    ylabel('Percentage of pixels where $$\hat{\sigma}_{G}(p) = \sigma$$', 'fontsize', 18, 'Interpreter', 'Latex');
    title({'Which scale produces maximum response?'; ['Noise level, \lambda = ' num2str(noise)]},...
        'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
    legend({'$$\sigma = 1$$', '$$\sigma = 2$$', '$$\sigma = 4$$', '$$\sigma = 8$$'}, 'fontsize', 18, 'Interpreter', 'Latex');
    exportfig([fig_dir 'synthetic\syn_lines_g2d_scales_' num2str(noise) '.png']);
end
%--------------------------------------------------------------------------
%%
% Percentage of oriented bands correctly predicted: synthetic lines with
% increasing noise: Gabor
exp_dir = 'C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\';
for noise = [0 1 2 3]
    data_dir = [exp_dir 'test\rician_' num2str(noise) '\1\'];

    %Load and reformat data
    responses_gabor = u_load([data_dir 'responses_gabor.mat']);
    responses_gabor = reshape(responses_gabor, [], 9, 4, 6);
    responses_gabor = real(responses_gabor(:, 5, :, :));
    
    %Load ground truth data
    load([data_dir 'true_labels.mat'], 'true_oris', 'true_widths');
    if ~noise
        true_oris = true_oris/2;
    end
    
    [~, predicted_bands, predicted_scales] = max_response_line(responses_gabor);
    clear responses_g2d;
    
    %Orientation
    num_bins = 50;
    ori_bins = linspace(0, pi, num_bins);
    counts = zeros(6, num_bins);
    for i_band = 1:6
        counts(i_band,:) = hist(true_oris(predicted_bands == i_band), ori_bins);
    end

    figure;
    plot(ori_bins, bsxfun(@rdivide, counts, sum(counts))', 'linewidth', 2);
    set(gca, 'ylim', [0 1], 'xlim', [0 pi], 'fontsize', 18);
    xlabel('Orientation of vessel - $$\phi(p)$$', 'fontsize', 18, 'fontsize', 18, 'Interpreter', 'Latex');
    ylabel('Percentage of pixels where $$\hat{\theta}_{Re}(p) = \theta$$', 'fontsize', 18, 'fontsize', 18, 'Interpreter', 'Latex');
    title({'Which oriented sub-band produces maximum response?'; ['Noise level, \lambda = ' num2str(noise)]},...
        'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
    legend({...
        '$$\theta = 0$$',...
        '$$\theta = ^{\pi}/_{6}$$',...
        '$$\theta = ^{2\pi}/_{6}$$',...
        '$$\theta = ^{3\pi}/_{6}$$',...
        '$$\theta = ^{4\pi}/_{6}$$',...
        '$$\theta = ^{5\pi}/_{6}$$'}, 'fontsize', 18, 'Interpreter', 'Latex');
    exportfig([fig_dir 'synthetic\syn_lines_gabor_ori_subbands_' num2str(noise) '.png']);
    
    counts = zeros(4, 16);
    width_bins = linspace(1, 8, 16);
    for i_scale = 1:4
        counts(i_scale,:) = hist(true_widths(predicted_scales == i_scale), width_bins);
    end

    figure; plot(width_bins, bsxfun(@rdivide, counts, sum(counts))', 'linewidth', 2); 
    xlabel('Width of vessel - $$w(p)$$', 'fontsize', 18, 'Interpreter', 'Latex');
    ylabel('Percentage of pixels where $$\hat{\sigma}_{Re}(p) = \sigma$$', 'fontsize', 18, 'Interpreter', 'Latex');
    title({'Which scale produces maximum response?'; ['Noise level, \lambda = ' num2str(noise)]},...
        'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
    legend({'$$\sigma = 1$$', '$$\sigma = 2$$', '$$\sigma = 4$$', '$$\sigma = 8$$'}, 'fontsize', 18, 'Interpreter', 'Latex');
    exportfig([fig_dir 'synthetic\syn_lines_gabor_scales_' num2str(noise) '.png']);
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
%Compare distributions of responses for line and background pixels -
%synthetic lines, g2d
exp_dir = 'C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\';

%Analytic
for noise = [0 1 2 3]
    data_dir = [exp_dir 'test\rician_' num2str(noise) '\1\'];

    %Load and reformat data
    responses_g2d = u_load([data_dir 'responses_g2d.mat']);
    responses_g2d = reshape(responses_g2d, [], 9, 4, 3);
    responses_g2d = -responses_g2d(:, 5, :, :);
    [predicted_fg_responses] = ...
        gaussian_2nd_derivative_line(responses_g2d(:,:,:,:));
    clear responses_g2d;
    
    responses_g2d = u_load([data_dir 'bg_responses_g2d.mat']);
    responses_g2d = reshape(responses_g2d, [], 9, 4, 3);
    responses_g2d = -responses_g2d(:, 5, :, :);
    [predicted_bg_responses] = ...
        gaussian_2nd_derivative_line(responses_g2d(:,:,:,:));
    clear responses_g2d;
    
    min_val = min(min(predicted_fg_responses), min(predicted_bg_responses));
    max_val = max(max(predicted_fg_responses), max(predicted_bg_responses));
    
    val_bins = linspace(min_val, max_val, 100);
    
    fg_counts = hist(predicted_fg_responses, val_bins);
    bg_counts = hist(predicted_bg_responses, val_bins);

    figure; hold on;
    plot(val_bins, fg_counts, 'r', 'linewidth', 2);
    plot(val_bins, bg_counts, 'b', 'linewidth', 2);
    xlabel('Steered filter response at maximal scale - $$\hat{I}_{G}(p)$$',...
        'fontsize', 18, 'Interpreter', 'Latex');
    ylabel('Percentage of pixels', 'fontsize', 18, 'Interpreter', 'Latex');
    title({'Distribution of \it G" \rm filter responses for synthetic lines'; ['Noise level, \lambda = ' num2str(noise)]},...
        'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
    legend({'Foreground pixels', 'Background pixels'},...
        'location', 'northeast', 'fontsize', 18, 'Interpreter', 'Latex');
    
    exportfig([fig_dir 'synthetic\syn_lines_g2d_responses_cdf_' num2str(noise) '.png']);
end
%--------------------------------------------------------------------------
%%
%Compare distributions of responses for line and background pixels -
%synthetic lines, Gabor
exp_dir = 'C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\';

%Analytic
for noise = [0 1 2 3]
    data_dir = [exp_dir 'test\rician_' num2str(noise) '\1\'];

    %Load and reformat data
    responses_gabor = u_load([data_dir 'responses_gabor.mat']);
    responses_gabor = reshape(responses_gabor, [], 9, 4, 6);
    responses_gabor = real(responses_gabor(:, 5, :, :));
    [predicted_fg_responses] = ...
        max_response_line(responses_gabor);
    clear responses_gabor;
    
    responses_gabor = u_load([data_dir 'bg_responses_gabor.mat']);
    responses_gabor = reshape(responses_gabor, [], 9, 4, 6);
    responses_gabor = real(responses_gabor(:, 5, :, :));
    [predicted_bg_responses] = ...
        max_response_line(responses_gabor(:,:,:,:));
    clear responses_gabor;
    
    min_val = min(min(predicted_fg_responses), min(predicted_bg_responses));
    max_val = max(max(predicted_fg_responses), max(predicted_bg_responses));
    
    val_bins = linspace(min_val, max_val, 100);
    
    fg_counts = hist(predicted_fg_responses, val_bins) / numel(predicted_fg_responses);
    bg_counts = hist(predicted_bg_responses, val_bins) / numel(predicted_bg_responses);

    figure; hold on;
    plot(val_bins, fg_counts, 'r', 'linewidth', 2);
    plot(val_bins, bg_counts, 'b', 'linewidth', 2);
    xlabel('Filter response at maximal scale and orientation - $$\hat{I}_{Re}(p)$$',...
        'fontsize', 18, 'Interpreter', 'Latex');
    ylabel('Percentage of pixels', 'fontsize', 18, 'Interpreter', 'Latex');
    title({'Distribution of Gabor filter responses for synthetic lines'; ['Noise level, \lambda = ' num2str(noise)]},...
        'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
    legend({'Foreground pixels', 'Background pixels'},...
        'location', 'northeast', 'fontsize', 18, 'Interpreter', 'Latex');
    
    exportfig([fig_dir 'synthetic\syn_lines_gabor_responses_cdf_' num2str(noise) '.png']);
end  
%--------------------------------------------------------------------------
%%
%Compare distributions of responses for line and background pixels -
%retina DRIVE, g2d
exp_dir = 'C:\isbe\asymmetry_project\experiments\DRIVE\comparing_representations\';
data_dir = [exp_dir 'test\1\'];

%Load and reformat data
responses_g2d = u_load([data_dir 'responses_g2d.mat']);
responses_g2d = reshape(responses_g2d, [], 9, 5, 3);
responses_g2d = responses_g2d(:, 5, :, :);
[predicted_fg_responses] = ...
    gaussian_2nd_derivative_line(responses_g2d(:,:,:,:));
clear responses_g2d;

responses_g2d = u_load([data_dir 'bg_responses_g2d.mat']);
responses_g2d = reshape(responses_g2d, [], 9, 5, 3);
responses_g2d = responses_g2d(:, 5, :, :);
[predicted_bg_responses] = ...
    gaussian_2nd_derivative_line(responses_g2d(:,:,:,:));
clear responses_g2d;

min_val = min(min(predicted_fg_responses), min(predicted_bg_responses));
max_val = max(max(predicted_fg_responses), max(predicted_bg_responses));

val_bins = linspace(min_val, max_val, 100);

fg_counts = hist(predicted_fg_responses, val_bins);
bg_counts = hist(predicted_bg_responses, val_bins);

figure; hold on;
plot(val_bins, fg_counts, 'r', 'linewidth', 2);
plot(val_bins, bg_counts, 'b', 'linewidth', 2);
xlabel('Steered filter response at maximal scale - $$\hat{I}_{G}(p)$$',...
    'fontsize', 18, 'Interpreter', 'Latex');
ylabel('Percentage of pixels', 'fontsize', 18, 'Interpreter', 'Latex');
title({'Distribution of \it G" \rm filter responses for synthetic lines'; 'Retinograms DRIVE data'},...
    'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
legend({'Foreground pixels', 'Background pixels'},...
    'location', 'northeast', 'fontsize', 18, 'Interpreter', 'Latex');

exportfig([fig_dir 'retina\ret_vessels_g2d_responses_cdf.png']);
%--------------------------------------------------------------------------
%%
%Compare distributions of responses for line and background pixels -
%retina DRIVE, Gabor
exp_dir = 'C:\isbe\asymmetry_project\experiments\DRIVE\comparing_representations\';
data_dir = [exp_dir 'test\1\'];

%Load and reformat data
responses_gabor = u_load([data_dir 'responses_gabor.mat']);

sample_idx = randperm(size(responses_gabor, 1));
sample_idx = sample_idx(1:2e4);

responses_gabor = reshape(responses_gabor(sample_idx,:), [], 9, 5, 6);
responses_gabor = -real(responses_gabor(:, 5, :, :));
[predicted_fg_responses] = ...
    max_response_line(responses_gabor);
clear responses_gabor;

responses_gabor = u_load([data_dir 'bg_responses_gabor.mat']);
responses_gabor = reshape(responses_gabor(sample_idx,:), [], 9, 5, 6);
responses_gabor = -real(responses_gabor(:, 5, :, :));
[predicted_bg_responses] = ...
    max_response_line(responses_gabor(:,:,:,:));
clear responses_gabor;

min_val = min(min(predicted_fg_responses), min(predicted_bg_responses));
max_val = max(max(predicted_fg_responses), max(predicted_bg_responses));

val_bins = linspace(min_val, max_val, 100);

fg_counts = hist(predicted_fg_responses, val_bins) / numel(predicted_fg_responses);
bg_counts = hist(predicted_bg_responses, val_bins) / numel(predicted_bg_responses);

figure; hold on;
plot(val_bins, fg_counts, 'r', 'linewidth', 2);
plot(val_bins, bg_counts, 'b', 'linewidth', 2);
xlabel('Filter response at maximal scale and orientation - $$\hat{I}_{Re}(p)$$',...
    'fontsize', 18, 'Interpreter', 'Latex');
ylabel('Percentage of pixels', 'fontsize', 18, 'Interpreter', 'Latex');
title({'Distribution of Gabor filter responses for synthetic lines'; 'Retinograms DRIVE data'},...
    'fontsize', 18, 'Interpreter', 'tex', 'fontname', 'times');
legend({'Foreground pixels', 'Background pixels'},...
    'location', 'northeast', 'fontsize', 18, 'Interpreter', 'Latex');

exportfig([fig_dir 'retina\ret_vessels_gabor_responses_cdf.png']);
%--------------------------------------------------------------------------
%%
%Prediction error of orientation for G2d: synthetic lines with increasing
%noise
exp_dir = 'C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\';

%Analytic
figure; a1 = gca; hold all;
for noise = [0 1 2 3]
    data_dir = [exp_dir 'test\rician_' num2str(noise) '\1\'];

    %Load and reformat data
    responses_g2d = u_load([data_dir 'responses_g2d.mat']);
    responses_g2d = reshape(responses_g2d, [], 9, 4, 3);
    responses_g2d = responses_g2d(:, 5, :, :);

    %Load ground truth data
    load([data_dir 'true_labels.mat'], 'true_oris', 'true_widths');
    complex_oris = complex(cos(2*true_oris), sin(2*true_oris));
       
    [~, predicted_theta, predicted_scales] = ...
        gaussian_2nd_derivative_line(responses_g2d);

    predicted_ori = complex(cos(2*predicted_theta), sin(2*predicted_theta));
    [~, err_stats] = ori_error(complex_oris, predicted_ori);

    plot(a1, err_stats.abs_percentiles, 1:100, 'linewidth', 3);
    %plot(a1, [0 err_stats.abs_percentiles(75)], [75 75], 'k--', 'linewidth', 1);
    %plot(a1, [err_stats.abs_percentiles(75) err_stats.abs_percentiles(75)], [0 75], 'k--', 'linewidth', 1);        
end
xlabel(a1, 'Absolute error in orientation prediction',...
    'fontsize', 18, 'fontname', 'times');
ylabel(a1, 'Percentage of pixels', 'fontsize', 18,...
    'fontname', 'times');
title(a1, 'G" cumulative distribution of orientation prediction errors',...
    'fontsize', 18, 'fontname', 'times');
legend(a1, {'$$\lambda = 0$$', '$$\lambda = 1$$', '$$\lambda = 2$$', '$$\lambda = 3$$'},...
    'location', 'southeast', 'fontsize', 18, 'Interpreter', 'Latex');
exportfig([fig_dir 'synthetic\syn_lines_g2d_noise_cdf.png']);
%%
%Random forest regressed
figure; a1 = gca; hold all;
for noise = [0 1 2 3]
    %Load ground truth data
    load([exp_dir 'test\rician_' num2str(noise) '\1\true_labels.mat'], 'true_oris', 'true_widths');
    complex_oris = complex(cos(2*true_oris), sin(2*true_oris));

    predicted_ori = u_load([exp_dir 'results\rician_' num2str(noise) '\1\orientation\g2d\orig\3\results.mat']);
    [~, err_stats] = ori_error(complex_oris, predicted_ori);

    plot(a1, err_stats.abs_percentiles, 1:100, 'linewidth', 3);
    %plot(a1, [0 err_stats.abs_percentiles(75)], [75 75], 'k--', 'linewidth', 1);
    %plot(a1, [err_stats.abs_percentiles(75) err_stats.abs_percentiles(75)], [0 75], 'k--', 'linewidth', 1);
        
end
xlabel(a1, 'Absolute error in orientation prediction',...
    'fontsize', 18, 'fontname', 'times');
ylabel(a1, 'Percentage of pixels',...
    'fontsize', 18, 'fontname', 'times');
title(a1, 'G" cumulative distribution of orientation prediction errors',...
    'fontsize', 18, 'fontname', 'times');
legend(a1, {'$$\lambda = 0$$', '$$\lambda = 1$$', '$$\lambda = 2$$', '$$\lambda = 3$$'},...
    'location', 'southeast', 'fontsize', 18, 'Interpreter', 'Latex');
exportfig([fig_dir 'synthetic\syn_lines_g2d_RF_noise_cdf.png']);
%--------------------------------------------------------------------------
%%
% Prediction of G2d for retinograms, all levels
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];
load('C:\isbe\asymmetry_project\experiments\DRIVE\g2d_orientation_prediction\v_responses.mat', 'v_responses');
load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
load([retroot '/width_maps/gt/all_gt_widths.mat'], 'all_gt_widths');

v_responses = reshape(v_responses, [], 1, 5, 3);
true_oris = angle(gt_orientations) / 2;
%
figure; hold all; a1 = gca;
for i_scale = [6 1:5]
    if i_scale < 6
        [~, predicted_theta] = gaussian_2nd_derivative_line(v_responses(:,:,i_scale,:));
    else
        [~, predicted_theta] = gaussian_2nd_derivative_line(v_responses(:,:,1:4,:));
    end
    predicted_ori = complex(cos(2*predicted_theta), sin(2*predicted_theta));
    [~, err_stats] = ori_error(predicted_ori, true_oris);
    display(['Scale: ' num2str(i_scale)]);
    display(err_stats);
    plot(a1, err_stats.abs_percentiles, 1:100, 'linewidth', 3);
end
xlabel(a1, 'Absolute error in orientation prediction', 'fontsize', 18, 'fontname', 'times');
ylabel(a1, 'Percentage of pixels', 'fontsize', 18, 'fontname', 'times');
legend(a1, {'All scales', '$$\sigma = 1$$', '$$\sigma = 2$$', '$$\sigma = 4$$', '$$\sigma = 8$$', '$$\sigma = 16$$'},...
    'location', 'southeast', 'fontsize', 18, 'Interpreter', 'Latex');
title(a1, {'G" cumulative distribution of orientation error'}, ...
    'fontsize', 18, 'fontname', 'times');

exportfig([fig_dir 'retina\ret_vessels_g2d_scales_cdf.png']);
%%
%G2d RF regressed prediction of retinogram orientation - DRIVE
retroot = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\training\';
pred_dir = [retroot 'predictions\orientation\rf_regression\'];
fov_mask_dir = [retroot 'fov_masks\'];
fg_mask_dir = [retroot 'vessel_masks\'];

rf_codes(6,:) = {'24069'};
rf_codes(1,:) = {'30469'};
rf_codes(2,:) = {'30470'};
rf_codes(3,:) = {'30471'};
rf_codes(4,:) = {'30472'};
rf_codes(5,:) = {'30473'};

figure; hold all; a1 = gca;
for i_scale = [6 1:5]
    [~, ~, err_stats] =...
        compute_image_orientation_errors([pred_dir rf_codes{i_scale,1} '\'], fg_mask_dir,...
        'gt_orientations', gt_orientations, 'fov_mask_dir', fov_mask_dir);
    plot(a1, err_stats.abs_percentiles, 1:100, 'linewidth', 3);
end

xlabel(a1, 'Absolute error in orientation prediction', 'fontsize', 18, 'fontname', 'times');
ylabel(a1, 'Percentage of pixels', 'fontsize', 18, 'fontname', 'times');
legend(a1, {'All scales', '$$\sigma = 1$$', '$$\sigma = 2$$', '$$\sigma = 4$$', '$$\sigma = 8$$', '$$\sigma = 16$$'},...
    'location', 'southeast', 'fontsize', 18, 'Interpreter', 'Latex');
title(a1, {'G" cumulative distribution of orientation error'}, ...
    'fontsize', 18, 'fontname', 'times');
exportfig([fig_dir 'retina\ret_vessels_g2d_RF_scales_cdf.png']);
%%
%G2d RF regressed prediction of retinogram orientation, all levels - synthetic data
exp_dir = 'C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\';
noise = 2;
%Analytic
figure; a1 = gca; hold all;
for i_scale = [5 1:4]
    if i_scale == 5
        scales = 1:4;
    else
        scales = i_scale;
    end
    data_dir = [exp_dir 'test\rician_' num2str(noise) '\1\'];

    %Load and reformat data
    responses_g2d = u_load([data_dir 'responses_g2d.mat']);
    responses_g2d = reshape(responses_g2d, [], 9, 4, 3);
    responses_g2d = responses_g2d(:, 5, scales, :);

    %Load ground truth data
    load([data_dir 'true_labels.mat'], 'true_oris', 'true_widths');
    complex_oris = complex(cos(2*true_oris), sin(2*true_oris));
       
    [~, predicted_theta] = ...
        gaussian_2nd_derivative_line(responses_g2d);

    predicted_ori = complex(cos(2*predicted_theta), sin(2*predicted_theta));
    [~, err_stats] = ori_error(complex_oris, predicted_ori);

    plot(a1, err_stats.abs_percentiles, 1:100, 'linewidth', 3);
    %plot(a1, [0 err_stats.abs_percentiles(75)], [75 75], 'k--', 'linewidth', 1);
    %plot(a1, [err_stats.abs_percentiles(75) err_stats.abs_percentiles(75)], [0 75], 'k--', 'linewidth', 1);        
end
xlabel(a1, 'Absolute error in orientation prediction', 'fontsize', 18, 'fontname', 'times');
ylabel(a1, 'Percentage of pixels', 'fontsize', 18, 'fontname', 'times');
legend(a1, {'All scales', '$$\sigma = 1$$', '$$\sigma = 2$$', '$$\sigma = 4$$', '$$\sigma = 8$$'},...
    'location', 'southeast', 'fontsize', 18, 'Interpreter', 'Latex');
title(a1, {'G" cumulative distribution of orientation error'}, ...
    'fontsize', 18, 'fontname', 'times');
exportfig([fig_dir 'synthetic\syn_lines_g2d_scales_cdf.png']);
%%
%Random forest regressed
figure; a1 = gca; hold all;
for i_scale = [5 1:4]
    %Load ground truth data
    load([exp_dir 'test\rician_' num2str(noise) '\1\true_labels.mat'], 'true_oris', 'true_widths');
    complex_oris = complex(cos(2*true_oris), sin(2*true_oris));

    if i_scale == 5
        predicted_ori = u_load([exp_dir 'results\rician_' num2str(noise) '\1\orientation\g2d\orig\3\results.mat']);
    else
        predicted_ori = u_load([exp_dir 'results\rician_' num2str(noise) '\1\orientation\g2d\levels\' num2str(i_scale) '\3\results.mat']);
    end
    [~, err_stats] = ori_error(complex_oris, predicted_ori);

    plot(a1, err_stats.abs_percentiles, 1:100, 'linewidth', 3);
    %plot(a1, [0 err_stats.abs_percentiles(75)], [75 75], 'k--', 'linewidth', 1);
    %plot(a1, [err_stats.abs_percentiles(75) err_stats.abs_percentiles(75)], [0 75], 'k--', 'linewidth', 1);
        
end
xlabel(a1, 'Absolute error in orientation prediction', 'fontsize', 18, 'fontname', 'times');
ylabel(a1, 'Percentage of pixels', 'fontsize', 18, 'fontname', 'times');
legend(a1, {'All scales', '$$\sigma = 1$$', '$$\sigma = 2$$', '$$\sigma = 4$$', '$$\sigma = 8$$'},...
    'location', 'southeast', 'fontsize', 18, 'Interpreter', 'Latex');
title(a1, {'G" cumulative distribution of orientation error'}, ...
    'fontsize', 18, 'fontname', 'times');
exportfig([fig_dir 'synthetic\syn_lines_g2d_RF_scales_cdf.png']);
%%
%G2d RF regressed prediction of retinogram orientation, all levels error vs
%width - retinogram DRIVE data
retroot = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\training\';
fov_mask_dir = [retroot 'fov_masks\'];
fg_mask_dir = [retroot 'vessel_masks\'];
pred_dir = [retroot 'predictions\orientation\rf_regression\'];
label_dir = [retroot 'orientations\'];

load([label_dir 'gt\all_gt_orientations.mat'], 'gt_orientations');
load([retroot 'width_maps\gt\all_gt_widths.mat'], 'all_gt_widths');
[uni_widths,~, width_idx] = unique(all_gt_widths);
rf_codes = {...
    '30469';...
	'30470';...
	'30471';...
	'30472';...
	'30473';...
	'24069'};

figure; hold all; a1 = gca;
for ii = [6 1:5]

    [orientation_errors] =...
        compute_image_orientation_errors([pred_dir rf_codes{ii,1} '\'], fg_mask_dir,...
        'gt_orientations', gt_orientations, 'fov_mask_dir', fov_mask_dir);
    valid_idx = ~isnan(orientation_errors);
    [width_centres, smoothed_abs_errs] =...
        kernel_smoother(all_gt_widths(valid_idx), 180*abs(orientation_errors(valid_idx))/pi );
    plot(a1, width_centres, smoothed_abs_errs, 'linewidth', 2);
end
set(gca, 'xlim', [1 8], 'ylim', [0 45])
xlabel(a1, 'Vessel width - $$w(p)$$', 'fontsize', 18, 'Interpreter', 'Latex');
ylabel(a1, 'Absolute error (degrees)', 'fontsize', 18, 'fontname', 'times');
legend(a1, {'All scales', '$$\sigma = 1$$', '$$\sigma = 2$$', '$$\sigma = 4$$', '$$\sigma = 8$$', '$$\sigma = 16$$'},...
    'location', 'southeast', 'fontsize', 18, 'Interpreter', 'Latex');
title(a1, {'Estimate of MAE in orientation prediction as a function of vessel width'}, ...
    'fontsize', 18, 'fontname', 'times');
exportfig([fig_dir 'retina\ret_vessels_g2d_RF_scales_v_width.pdf']);
%%
%G2d analytic prediction of retinogram orientation, all levels error vs
%width - retinogram DRIVE data
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];
load('C:\isbe\asymmetry_project\experiments\DRIVE\g2d_orientation_prediction\v_responses.mat', 'v_responses');
load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
load([retroot '/width_maps/gt/all_gt_widths.mat'], 'all_gt_widths');

v_responses = reshape(v_responses, [], 1, 5, 3);
true_oris = angle(gt_orientations) / 2;
%
figure; hold all; a1 = gca;
for i_scale = [6 1:5]
    if i_scale < 6
        [~, predicted_theta] = gaussian_2nd_derivative_line(v_responses(:,:,i_scale,:));
    else
        [~, predicted_theta] = gaussian_2nd_derivative_line(v_responses(:,:,1:4,:));
    end
    predicted_ori = complex(cos(2*predicted_theta), sin(2*predicted_theta));
    [orientation_errors] = ori_error(predicted_ori, true_oris);
    valid_idx = ~isnan(orientation_errors);
    [width_centres, smoothed_abs_errs] =...
        kernel_smoother(all_gt_widths(valid_idx), 180*abs(orientation_errors(valid_idx))/pi );
    plot(a1, width_centres, smoothed_abs_errs, 'linewidth', 2);
    
    %figure; boxplot(abs(orientation_errors), width_idx);
end
set(gca, 'xlim', [1 8], 'ylim', [0 45])
xlabel(a1, 'Line width - $$w(p)$$', 'fontsize', 18, 'Interpreter', 'Latex');
ylabel(a1, 'Absolute error (degrees)', 'fontsize', 18, 'fontname', 'times');
legend(a1, {'All scales', '$$\sigma = 1$$', '$$\sigma = 2$$', '$$\sigma = 4$$', '$$\sigma = 8$$', '$$\sigma = 16$$'},...
    'location', 'southeast', 'fontsize', 18, 'Interpreter', 'Latex');
title(a1, {'Estimate of MAE in orientation prediction as a function of vessel width'}, ...
    'fontsize', 18, 'fontname', 'times');
exportfig([fig_dir 'retina\ret_vessels_g2d_scales_v_width.pdf']);
%%
%G2d RF regressed prediction of retinogram orientation, all levels estimate of error vs width - synthetic data
exp_dir = 'C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\';
noise = 2;
%Analytic
figure; a1 = gca; hold all;
for i_scale = [5 1:4]
    if i_scale == 5
        scales = 1:4;
    else
        scales = i_scale;
    end
    data_dir = [exp_dir 'test\rician_' num2str(noise) '\1\'];

    %Load and reformat data
    responses_g2d = u_load([data_dir 'responses_g2d.mat']);
    responses_g2d = reshape(responses_g2d, [], 9, 4, 3);
    responses_g2d = responses_g2d(:, 5, scales, :);

    %Load ground truth data
    load([data_dir 'true_labels.mat'], 'true_oris', 'true_widths');
    complex_oris = complex(cos(2*true_oris), sin(2*true_oris));
       
    [~, predicted_theta] = ...
        gaussian_2nd_derivative_line(responses_g2d);

    predicted_ori = complex(cos(2*predicted_theta), sin(2*predicted_theta));
    [orientation_errors] = ori_error(predicted_ori, true_oris);
    valid_idx = ~isnan(orientation_errors);
    [width_centres, smoothed_abs_errs] =...
        kernel_smoother(true_widths(valid_idx), 180*abs(orientation_errors(valid_idx))/pi );
    plot(a1, width_centres, smoothed_abs_errs, 'linewidth', 2);       
end
set(gca, 'xlim', [1 8], 'ylim', [0 45])
xlabel(a1, 'Line width - $$w(p)$$', 'fontsize', 18, 'Interpreter', 'Latex');
ylabel(a1, 'Absolute error (degrees)', 'fontsize', 18, 'fontname', 'times');
legend(a1, {'All scales', '$$\sigma = 1$$', '$$\sigma = 2$$', '$$\sigma = 4$$', '$$\sigma = 8$$'},...
    'location', 'southeast', 'fontsize', 18, 'Interpreter', 'Latex');
title(a1, {'Estimate of MAE in orientation prediction as a function of line width'}, ...
    'fontsize', 18, 'fontname', 'times');
exportfig([fig_dir 'synthetic\syn_lines_g2d_scales_v_width.png']);
%%
%Something else
retroot = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\training\';
fov_mask_dir = [retroot 'fov_masks\'];
fg_mask_dir = [retroot 'vessel_masks\'];
pred_dir = [retroot 'predictions\orientation\rf_regression\'];
label_dir = [retroot 'orientations\'];

load([label_dir 'gt\all_gt_orientations.mat'], 'gt_orientations');
load([retroot 'width_maps\gt\all_gt_widths.mat'], 'all_gt_widths');
[uni_widths,~, width_idx] = unique(all_gt_widths);

figure; hold all; a1 = gca;

[orientation_errors] =...
    compute_image_orientation_errors([pred_dir '24308\'], fg_mask_dir,...
    'gt_orientations', gt_orientations, 'fov_mask_dir', fov_mask_dir);
valid_idx = ~isnan(orientation_errors);
[width_centres, smoothed_abs_errs] =...
    kernel_smoother(all_gt_widths(valid_idx), 180*abs(orientation_errors(valid_idx))/pi );
plot(a1, width_centres, smoothed_abs_errs, 'linewidth', 2);

clear orientation_errors;
load('C:\isbe\asymmetry_project\experiments\DRIVE\gabor_orientation_prediction\v_responses.mat', 'v_responses');

v_responses = reshape(v_responses, [], 1, 5, 6);
true_oris = angle(gt_orientations) / 2;

rad_angles = pi/2 + pi*(0:5)'/6;
complex_angles = complex(cos(2*rad_angles), sin(2*rad_angles));
[~, max_band] = max_response_line(v_responses(:,:,1:4,:));

predicted_ori = complex_angles(max_band);
[orientation_errors] = ori_error(predicted_ori, true_oris);
valid_idx = ~isnan(orientation_errors);
[width_centres, smoothed_abs_errs] =...
    kernel_smoother(all_gt_widths(valid_idx), 180*abs(orientation_errors(valid_idx))/pi );
plot(a1, width_centres, smoothed_abs_errs, 'linewidth', 2);

set(gca, 'xlim', [1 8], 'ylim', [0 30])
xlabel(a1, 'Vessel width - $$w(p)$$', 'fontsize', 18, 'Interpreter', 'Latex');
ylabel(a1, 'Absolute error (degrees)', 'fontsize', 18, 'fontname', 'times');
legend(a1, {'RF predictions', 'Analytic predictions'},...
    'location', 'southeast', 'fontsize', 18, 'Interpreter', 'Latex');
title(a1, {'Estimate of MAE in orientation prediction as a function of vessel width'}, ...
    'fontsize', 18, 'fontname', 'times');
%exportfig([fig_dir 'retina\ret_vessels_gabor_scales_v_width.pdf']);
%%
%Random forest regressed
noise = 2;
figure; a1 = gca; hold all;
for i_scale = [5 1:4]
    %Load ground truth data
    load([exp_dir 'test\rician_' num2str(noise) '\1\true_labels.mat'], 'true_oris', 'true_widths');
    complex_oris = complex(cos(2*true_oris), sin(2*true_oris));

    if i_scale == 5
        predicted_ori = u_load([exp_dir 'results\rician_' num2str(noise) '\1\orientation\g2d\orig\3\results.mat']);
    else
        predicted_ori = u_load([exp_dir 'results\rician_' num2str(noise) '\1\orientation\g2d\levels\' num2str(i_scale) '\3\results.mat']);
    end
    [~, err_stats] = ori_error(complex_oris, predicted_ori);

    [orientation_errors] = ori_error(predicted_ori, true_oris);
    valid_idx = ~isnan(orientation_errors);
    [width_centres, smoothed_abs_errs] =...
        kernel_smoother(true_widths(valid_idx), 180*abs(orientation_errors(valid_idx))/pi );
    plot(a1, width_centres, smoothed_abs_errs, 'linewidth', 2);  
        
end
set(gca, 'xlim', [1 8], 'ylim', [0 45])
xlabel(a1, 'Line width - $$w(p)$$', 'fontsize', 18, 'Interpreter', 'Latex');
ylabel(a1, 'Absolute error (degrees)', 'fontsize', 18, 'fontname', 'times');
legend(a1, {'All scales', '$$\sigma = 1$$', '$$\sigma = 2$$', '$$\sigma = 4$$', '$$\sigma = 8$$'},...
    'location', 'southeast', 'fontsize', 18, 'Interpreter', 'Latex');
title(a1, {'Estimate of MAE in orientation prediction as a function of line width'}, ...
    'fontsize', 18, 'fontname', 'times');
exportfig([fig_dir 'synthetic\syn_lines_g2d_RF_scales_v_width.png']);
%%
%Something else
retroot = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\';
fov_mask_dir = [retroot 'fov_masks\'];
fg_mask_dir = [retroot 'vessel_masks\'];
pred_dir = [retroot 'predictions\orientation\rf_regression\'];
label_dir = [retroot 'orientations\'];

load([label_dir 'gt\all_gt_orientations.mat'], 'gt_orientations');

rf_codes = {...
    '88171';...
	'24308';...
	'42813';...
	'13295'};

figure; hold all; a1 = gca;

for ii = 1:4

    [orientation_errors predicted_orientations] =...
        compute_image_orientation_errors([pred_dir rf_codes{ii} '\'], fg_mask_dir,...
        'gt_orientations', gt_orientations, 'fov_mask_dir', fov_mask_dir);
    valid_idx = ~isnan(orientation_errors);
    ori_variances = abs(predicted_orientations);
    
    [var_centres, smoothed_abs_errs] =...
        kernel_smoother(ori_variances(valid_idx), 180*abs(orientation_errors(valid_idx))/pi );
    plot(a1, var_centres, smoothed_abs_errs, 'linewidth', 2);
    
end
legend({'Gaussian + RF', 'Gabor + RF',  'DT-CWT + RF', 'Monogenic + RF'},...
    'location', 'northeast','fontsize', 18, 'fontname', 'times');
xlabel(a1, 'Angular variance of predicted orientation - $$|t_{est}(p)|$$', 'fontsize', 18, 'Interpreter', 'Latex');
ylabel(a1, 'Absolute error (degrees)', 'fontsize', 18, 'fontname', 'times');

title(a1, {'Estimate of MAE in orientation prediction'; 'as a function of predicted angular variance'}, ...
    'fontsize', 18, 'fontname', 'times');
exportfig([fig_dir 'retina\drive_RF_err_v_var.png']);

%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
retroot = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\predictions\detection\rf_classification\';
%Segmentation maps
prediction_image = u_load([retroot '88122\02_test_pred.mat']);
write_im_from_colormap(prediction_image, [fig_dir 'retina\02_segmentation_gabor.png'], gray(256));
write_im_from_colormap(1-prediction_image, [fig_dir 'retina\02_segmentation_gabor_inv.png'], gray(256));
write_im_from_colormap(prediction_image > .75, [fig_dir 'retina\02_segmentation_gabor_thresh.png'], gray(256));
%%
%Segmentation maps
prediction_image = u_load([retroot '9160\02_test_pred.mat']);
write_im_from_colormap(prediction_image, [fig_dir 'retina\02_segmentation_g2d.png'], gray(256));
write_im_from_colormap(1-prediction_image, [fig_dir 'retina\02_segmentation_g2d_inv.png'], gray(256));
write_im_from_colormap(prediction_image > .75, [fig_dir 'retina\02_segmentation_g2d_thresh.png'], gray(256));
prediction_image = u_load([retroot '12393\02_test_pred.mat']);
write_im_from_colormap(prediction_image, [fig_dir 'retina\02_segmentation_g12d.png'], gray(256));
write_im_from_colormap(1-prediction_image, [fig_dir 'retina\02_segmentation_g12d_inv.png'], gray(256));
write_im_from_colormap(prediction_image > .75, [fig_dir 'retina\02_segmentation_g12d_thresh.png'], gray(256));

%Example at optic disk
prediction_image = u_load([retroot '9160\02_test_pred.mat']);
write_im_from_colormap(1-prediction_image(210:340, 375:505), [fig_dir 'retina\02_optic_g2d_inv.png'], gray(256));
prediction_image = u_load([retroot '12393\02_test_pred.mat']);
write_im_from_colormap(1-prediction_image(210:340, 375:505), [fig_dir 'retina\02_optic_g12d_inv.png'], gray(256));
ret = u_load('C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\images\02_test.mat');
imwrite(ret(210:340, 375:505,:), [fig_dir 'retina\02_optic.png']);
%%
%Best/worst result in DRIVE test set
%Segmentation maps
rf_codes = cell(4,2);
rf_codes( 1,:) = {'9162',  'dt'};
rf_codes( 2,:) = {'88122', 'gabor'};
rf_codes( 3,:) = {'88126', 'gh2da'};
rf_codes( 4,:) = {'12394', 'mono'};
retroot = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\predictions\detection\rf_classification\';
for i_im = 8%[8 19]
    ret = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\images\' zerostr(i_im,2) '_test.mat']);
    imwrite(ret, [fig_dir 'retina\' zerostr(i_im,2) '_DRIVE_ret.png']);
    for i_rf = 1:4
        prediction_image = u_load([retroot rf_codes{i_rf,1} '\' zerostr(i_im,2) '_test_pred.mat']);
        write_im_from_colormap(prediction_image, [fig_dir 'retina\' zerostr(i_im,2) '_DRIVE_segmentation_' rf_codes{i_rf,2} '.png'], gray(256));
        write_im_from_colormap(1-prediction_image, [fig_dir 'retina\' zerostr(i_im,2) '_DRIVE_segmentation_' rf_codes{i_rf,2} '_inv.png'], gray(256));
    end
end
%%
% Biggest difference resampling makes
pred_dir = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\training\predictions\detection\rf_classification\';
p1 = u_load([pred_dir '24066\34_training_pred.mat']);
p2 = u_load([pred_dir '24303\34_training_pred.mat']);
write_im_from_colormap(p1 - p2, [fig_dir 'retina\34_resampling_difference.png'], jet(256));

ret = u_load('C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\training\images\34_training.mat');
imwrite(ret, [fig_dir 'retina\34_DRIVE_ret.png']);
%%
%Orientation_maps
retroot = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\predictions\orientation\rf_regression\';
prediction_image = u_load([retroot '24308\02_test_pred.mat']);
imwrite(complex2rgb(prediction_image), [fig_dir 'retina\02_orientation_mag.png']);
%%
%Best/worst result in fibre test set
%Segmentation maps
rf_codes = cell(4,2);
rf_codes( 1,:) = {'88443',  'dt'};
rf_codes( 2,:) = {'88311', 'gabor'};
rf_codes( 3,:) = {'88444', 'gh2da'};
rf_codes( 4,:) = {'88445', 'mono'};
fibre_root = 'C:\isbe\asymmetry_project\data\fibre\test\';
image_dir = [fibre_root 'images\'];
pred_dir = [fibre_root 'predictions\detection\rf_classification\'];

f_list = dir([image_dir '\*.mat']);
p_list = dir([pred_dir rf_codes{1,1} '\*.mat']);

for i_im = 3
    ret = u_load([image_dir f_list(i_im).name]);
    imwrite(ret, [fig_dir 'fibre\' zerostr(i_im,2) '_fibre_ccm.png']);
    for i_rf = 1:4
        prediction_image = u_load([pred_dir rf_codes{i_rf,1} '\' p_list(i_im).name]);
        write_im_from_colormap(prediction_image, [fig_dir 'fibre\' zerostr(i_im,2) '_fibre_segmentation_' rf_codes{i_rf,2} '.png'], gray(256));
        write_im_from_colormap(1-prediction_image, [fig_dir 'fibre\' zerostr(i_im,2) '_fibre_segmentation_' rf_codes{i_rf,2} '_inv.png'], gray(256));
    end
end
%%
%--------------------------------------------------------------------------
% DRIVE ROC curves on test set
rf_codes = cell(5,2);
rf_codes( 1,:) = {'rf_classification\88126', 'Gaussian (G,H) + RF (A_z='};
rf_codes( 2,:) = {'rf_classification\88122', 'Gabor (Re,Im) + RF (A_z='};
rf_codes( 3,:) = {'rf_classification\9162',  'DT-CWT (Re,Im) + RF (A_z='};
rf_codes( 4,:) = {'rf_classification\12394', 'Monogenic + RF (A_z='};
rf_codes( 5,:) = {'staal', 'Staal 2004 (A_z='};
rf_codes( 6,:) = {'', '2nd manual segmentation'};

retroot = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\';
pred_dir = [retroot 'predictions\detection\'];
%
f1 = figure;
graph(f1);
set(gca,'box','on','fontsize', 8); 
axis equal; axis([0 1 0 1]); hold all;
title('ROC curves for vessel segmentation, DRIVE database');
xlabel('FPF');
ylabel('TPF');

for ii = 1:5
    load([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'])
    plot(roc_pts(:,1), roc_pts(:,2), '-', 'linewidth', 1.0);
    rf_codes{ii,2} = [rf_codes{ii,2} num2str(auc, 3) ')'];
end
plot(0.0275, 0.775, 'k*', 'markersize', 10);
legend(rf_codes(:,2), 'location', 'southeast');

exportfig([fig_dir 'retina\DRIVE_test_detection_roc.png']);
axis([0 0.5 0.5 1]); 
exportfig([fig_dir 'retina\DRIVE_test_detection_roc_zoom.png']);
%%
%--------------------------------------------------------------------------
% STARE ROC curves on test set
rf_codes = cell(4,2);
rf_codes( 1,:) = {'rf_classification\gh2da', 'Gaussian (G,H) + RF (A_z='};
rf_codes( 2,:) = {'rf_classification\gabor', 'Gabor (Re,Im) + RF (A_z='};
rf_codes( 3,:) = {'rf_classification\dt',  'DT-CWT (Re,Im) + RF (A_z='};
rf_codes( 4,:) = {'rf_classification\mono', 'Monogenic + RF (A_z='};
%rf_codes( 5,:) = {'staal', 'Staal 2004 (A_z='};
%rf_codes( 6,:) = {'', '2nd manual segmentation'};

retroot = 'C:\isbe\asymmetry_project\data\retinograms\STARE\training\';
pred_dir = [retroot 'predictions\detection\'];
%
f1 = figure;
graph(f1);
set(gca,'box','on','fontsize', 8); 
axis equal; axis([0 1 0 1]); hold all;
title('ROC curves for vessel segmentation, STARE database');
xlabel('FPF');
ylabel('TPF');

for ii = 1:4
    load([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'])
    plot(roc_pts(:,1), roc_pts(:,2), '-', 'linewidth', 1.0);
    rf_codes{ii,2} = [rf_codes{ii,2} num2str(auc, 3) ')'];
end
%plot(0.0275, 0.775, 'k*', 'markersize', 10);
legend(rf_codes(:,2), 'location', 'southeast');

exportfig([fig_dir 'retina\STARE_detection_roc.png']);
axis([0 0.5 0.5 1]); 
exportfig([fig_dir 'retina\STARE_detection_roc_zoom.png']);
%%
%--------------------------------------------------------------------------
% DRIVE CDF curves on test set
rf_codes = cell(4,2);
rf_codes( 1,:) = {'rf_regression\88171', 'Gaussian (G,H) + RF (mae='};
rf_codes( 2,:) = {'rf_regression\24308', 'Gabor (Re,Im) + RF (mae='};
rf_codes( 3,:) = {'rf_regression\42813',  'DT-CWT (Re,Im) + RF (mae='};
rf_codes( 4,:) = {'rf_regression\13295', 'Monogenic + RF (mae='};

retroot = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\';
pred_dir = [retroot 'predictions\orientation\'];
%
f1 = figure;
graph(f1);
set(gca,'box','on','fontsize', 8); 
axis equal; axis([0 90 0 100]); hold all;
title('CDF of absolute error in vessel orientation prediction, DRIVE database');
xlabel('Absolute error (degrees)');
ylabel('Cumulative %');

for ii = 1:4
    load([pred_dir rf_codes{ii,1} '\errors\error_stats.mat']);
    plot(error_stats.abs_percentiles, 1:100, '-', 'linewidth', 1.0);
    rf_codes{ii,2} = [rf_codes{ii,2} num2str(error_stats.abs_median, 3) ')'];
end
legend(rf_codes(:,2), 'location', 'southeast');

exportfig([fig_dir 'retina\DRIVE_test_orientation_cdf.png']);
axis([0 45 50 100]); 
exportfig([fig_dir 'retina\DRIVE_test_orientation_cdf_zoom.png']);
%%
%--------------------------------------------------------------------------
% STARE CDF curves on test set
rf_codes = cell(4,2);
rf_codes( 1,:) = {'rf_regression\88171', 'Gaussian (G,H) + RF (mae='};
rf_codes( 2,:) = {'rf_regression\24308', 'Gabor (Re,Im) + RF (mae='};
rf_codes( 3,:) = {'rf_regression\42813',  'DT-CWT (Re,Im) + RF (mae='};
rf_codes( 4,:) = {'rf_regression\13295', 'Monogenic + RF (mae='};

retroot = 'C:\isbe\asymmetry_project\data\retinograms\STARE\training\';
pred_dir = [retroot 'predictions\orientation\'];
%
f1 = figure;
graph(f1);
set(gca,'box','on','fontsize', 8); 
axis equal; axis([0 90 0 100]); hold all;
title('CDF of absolute error in vessel orientation prediction, STARE database');
xlabel('Absolute error (degrees)');
ylabel('Cumulative %');

for ii = 1:4
    load([pred_dir rf_codes{ii,1} '\errors\error_stats.mat']);
    plot(error_stats.abs_percentiles, 1:100, '-', 'linewidth', 1.0);
    rf_codes{ii,2} = [rf_codes{ii,2} num2str(error_stats.abs_median, 3) ')'];
end
legend(rf_codes(:,2), 'location', 'southeast');

exportfig([fig_dir 'retina\STARE_orientation_cdf.png']);
axis([0 45 50 100]); 
exportfig([fig_dir 'retina\STARE_orientation_cdf_zoom.png']);
