fig_root = 'K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\';
%%
clear, clc
data_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\comparing_representations\';
exp_dir = 'C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\';

dx = 64;
dy = 64;
cx = 32;
cy = 32;
warning('off', 'ASYM:unexpectedArgument');
%
%Make a load of data for bothing training and testing purposes
% Set up arguments for each decomposition type
d_args{1}.decomp_type = {'g1d'};
d_args{1}.sigma_range = [1 2 4];
            
d_args{2}.decomp_type = {'g2d'};
d_args{2}.sigma_range = [1 2 4];
            
num_decomps = 2;

%set common parameters for all decomp types and then compute vector sizes
D = zeros(num_decomps,1);
for i_decomp = 1:num_decomps
    d_args{i_decomp}.win_size = 3;
    d_args{i_decomp}.normalise = 0;
    d_args{i_decomp}.pca = [];
    
    D(i_decomp) = get_samples_per_channel(d_args{i_decomp});
end

%%
%--------------------------------------------------------------------------
% ***************** GENERATE DATA *****************************************
%--------------------------------------------------------------------------
%Sample properties of line
noise_level = 2;
pts_per_im = 10;
line_width = sample_uniform([1 8]);
line_contrast = sample_uniform([1 8]);
line_ori = sample_uniform([0 180]);
line_rad = pi * line_ori / 180;

%Generate line
[line, label, label_centre] =...
    create_ellipse_bar(line_width/2, line_contrast, line_ori, dy, dx, cx, cy);
test_image = ricernd(1 + line, 2);

label([1 end],:) = 0;
label(:, [1 end]) = 0;

%Select some random pixels in the signal image
shuffle = randperm(sum(label(:)));
idx = find(label);
r_idx = idx(shuffle(1:pts_per_im));
[r_rows r_cols] = ind2sub([dy dx], r_idx);

%For ecah decomposition compute the responses
g1d_responses = compute_filter_responses(line, d_args{1});
g2d_responses = compute_filter_responses(line, d_args{2});

figure; imgray(test_image);
%%
figure;
for i_level = 1:4
    subplot(2,4,i_level); imgray(g1d_responses.g1d(:,:,i_level,1));
    subplot(2,4,i_level+4); imgray(g1d_responses.g1d(:,:,i_level,2));
end
figure;
for i_level = 1:4
    subplot(3,4,i_level); imgray(g2d_responses.g2d(:,:,i_level,1));
    subplot(3,4,i_level+4); imgray(g2d_responses.g2d(:,:,i_level,2));
    subplot(3,4,i_level+8); imgray(g2d_responses.g2d(:,:,i_level,3));
end
%%
figure;
for i_level = 1:4
    g1dr = gaussian_1st_derivative_gradient(g1d_responses.g1d(:,:,i_level,:));
    g2dr = gaussian_2nd_derivative_line(g2d_responses.g2d(:,:,i_level,:));
    subplot(2,4,i_level); imgray(g1dr);
    subplot(2,4,i_level+4); imgray(g2dr);
end
%%
load('C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\training\rician_4\true_labels.mat')
g1df = u_load('C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\training\rician_4\responses_g1d.mat');
g2df = u_load('C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\training\rician_4\responses_g2d.mat');
g1db = u_load('C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\training\rician_4\bg_responses_g1d.mat');
g2db = u_load('C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\training\rician_4\bg_responses_g2d.mat');
%%
g1df = reshape(g1df, [], 9, 5, 2);
g2df = reshape(g2df, [], 9, 5, 3);
g1db = reshape(g1db, [], 9, 5, 2);
g2db = reshape(g2db, [], 9, 5, 3);
%%
for i_level = 1:4
    g1df_i = squeeze(g1df(:,5,i_level,:));
    g1db_i = squeeze(g2df(:,5,i_level,:));
    figure; 
    subplot(1,2,1); hold all;
    plot(g1df_i(:,1), g1df_i(:,2), 'rx');
    plot(g1db_i(:,1), g1db_i(:,2), 'gx');
    subplot(1,2,2); hold all;
    plot(g1df_i(centre_idx,1), g1df_i(centre_idx,2), 'rx');
    plot(g1db_i(centre_idx,1), g1db_i(centre_idx,2), 'gx');
end
%%
num_ims = 200;
pts_per_im = 10;

decomps = 1:num_decomps;

for data_type = {'training', 'test'}
    
    for noise_level = 2%[0 1 2 4];
        mkdir([exp_dir data_type{1} '/rician_' num2str(noise_level)]);

        %Pre-allocate space for the true line parameters
        true_oris = zeros(num_ims*pts_per_im,1);
        true_cons = zeros(num_ims*pts_per_im,1);
        true_widths = zeros(num_ims*pts_per_im,1);
        centre_idx = false(num_ims*pts_per_im,1);
        predicted_ori_1 = zeros(num_ims*pts_per_im,1);
        predicted_ori_2 = zeros(num_ims*pts_per_im,1);

        if strcmpi(data_type{1}, 'training')
            rng(1000, 'twister');
            display('setting twister to 1000');
        else
            rng(2000, 'twister');
            display('setting twister to 2000');
        end
        
        for ii = 1:num_ims

            display(['Testing image ' num2str(ii)]);

            %Sample properties of line
            line_width = sample_uniform([1 8]);
            line_contrast = sample_uniform([1 8]);
            line_ori = sample_uniform([0 180]);
            line_rad = pi * line_ori / 180;

            %Generate line
            [line, label, label_centre] =...
                create_ellipse_bar(line_width/2, line_contrast, line_ori, dy, dx, cx, cy);
%             [line, label, label_centre] =...
%                 create_sin_step(line_width/2, line_contrast, line_ori, dy, dx, cx, cy);
            test_image = ricernd(1 + line, noise_level);

            label([1 end],:) = 0;
            label(:, [1 end]) = 0;

            %Select some random pixels in the signal image
            shuffle = randperm(sum(label(:)));
            idx = find(label);
            r_idx = idx(shuffle(1:pts_per_im));
            [r_rows r_cols] = ind2sub([dy dx], r_idx);
            sample_idx = ((ii-1)*pts_per_im + 1):(ii*pts_per_im);
            
            if ii < 3
                figure; imgray(test_image);
                plot(r_cols, r_rows, 'rx');
                plot(r_cols(label_centre(r_idx)), r_rows(label_centre(r_idx)), 'gx');
            end

            %Save the line parameters
            true_oris(sample_idx,:) = line_rad;
            true_cons(sample_idx,:) = line_contrast;
            true_widths(sample_idx,:) = line_width;
            centre_idx(sample_idx,:) = label_centre(r_idx);

            %For ecah decomposition compute the responses
            g1d_responses = compute_filter_responses(test_image, d_args{1});
            g2d_responses = compute_filter_responses(test_image, d_args{2});
            
            [~, ori_1] = gaussian_1st_derivative_gradient(g1d_responses.g1d);
            [~, ori_2] = gaussian_2nd_derivative_line(g2d_responses.g2d);
            
            predicted_ori_1(sample_idx,:) = ori_1(r_idx);
            predicted_ori_2(sample_idx,:) = ori_2(r_idx);

        end
        
        complex_ori = complex(cos(2*true_oris), sin(2*true_oris));
        [~, err_stats_1] = ori_error(complex_ori(centre_idx), predicted_ori_1(centre_idx));
        [~, err_stats_2] = ori_error(complex_ori(centre_idx), predicted_ori_2(centre_idx));
        display(err_stats_1);
        display(err_stats_2);

        figure; hold all;
        title(['Noise level: ' num2str(noise_level)]);
        plot(err_stats_1.abs_percentiles, 1:100);
        plot(err_stats_2.abs_percentiles, 1:100);
        legend({'G''', 'G"'}, 'location', 'southeast');
        
    end
end
%%
g1d_predictions = u_load('C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\results\orientation\g1d\orig\1\results.mat');
g2d_predictions = u_load('C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\results\orientation\g2d\orig\1\results.mat');
load('C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\test\rician_4\true_labels.mat');
complex_ori = complex(cos(2*true_oris), sin(2*true_oris));

[~, err_stats_1] = ori_error(complex_ori(centre_idx), g1d_predictions(centre_idx));
[~, err_stats_2] = ori_error(complex_ori(centre_idx), g2d_predictions(centre_idx));
display(err_stats_1);
display(err_stats_2);

figure; hold all;
title(['Noise level: ' num2str(2)]);
plot(err_stats_1.abs_percentiles, 1:100);
plot(err_stats_2.abs_percentiles, 1:100);
legend({'G''', 'G"'}, 'location', 'southeast');
%%
% predictor_dt = u_load('C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\rfs\rician_1\1\detection\dt\orig\1\predictor.mat');
% predictor_dt.tree_dir = 'C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\rfs\rician_1\1\detection\dt\orig\1\trees\';
predictor_g2 = u_load('C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\rfs\rician_1\1\detection\g2d\orig\1\predictor.mat');
predictor_g2.tree_dir = 'C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\rfs\rician_1\1\detection\g2d\orig\1\trees\';

% d_args_dt.decomp_type = {'dt'};
% d_args_dt.levels = 1:4;
% d_args_dt.feature_shape = 'rect';
% d_args_dt.feature_type = 'conj';
% d_args_dt.do_max = 0;
% d_args_dt.rotate = 0;
% d_args_dt.use_nag = 0;
% d_args_dt.win_size = 1;
% d_args_dt.normalise = 0;
% d_args_dt.pca = [];
% d_args_dt.rgb_channel = 'rgb';

d_args_g2.decomp_type = {'g2d'};
d_args_g2.sigma_range = [1 2 4 8];
d_args_g2.feature_shape = 'rect';
d_args_g2.win_size = 1;
d_args_g2.normalise = 0;
d_args_g2.pca = [];
d_args_g2.rgb_channel = 'rgb';

noise_level = 4; 
pts_per_im = 10;
num_ims = 2;
dx = 64;
dy = 64;
cx = 32;
cy = 32;

class_labels = [false(pts_per_im*num_ims,1) true(pts_per_im*num_ims,1)];
% class_predictions_dt = zeros(pts_per_im*num_ims,2);
class_predictions_g2 = zeros(pts_per_im*num_ims,2);

for ii = 1:num_ims
    line_width = sample_uniform([1 8]);
    line_contrast = sample_uniform([1 8]);
    line_ori = sample_uniform([0 180]);
    line_rad = pi * line_ori / 180;

    %Generate line
    [line, label, label_centre] =...
        create_ellipse_bar(line_width/2, line_contrast, line_ori, dy, dx, cx, cy);
    test_image = ricernd(1 + line, noise_level);
    
    fov_mask = true(size(label));
    fov_mask([1:16 end-15:end],:) = 0;
    fov_mask(:, [1:16 end-15:end]) = 0;

    %Select some random pixels in the signal image
    shuffle = randperm(sum(label(:) & fov_mask(:)));
    idx = find(label & fov_mask);
    fg_idx = idx(shuffle(1:pts_per_im));

    %Select some random bg pixels in the signal image
    shuffle = randperm(sum(~label(:) & fov_mask(:)));
    idx = find(~label & fov_mask);
    bg_idx = idx(shuffle(1:pts_per_im));
    
    sample_idx = ((ii-1)*pts_per_im + 1):(ii*pts_per_im);   
    
%     line_prob_dt = predict_image(...
%         'image_in', test_image,... % the mandatory arguments
%         'decomposition_args', d_args_dt,...
%         'predictor', predictor_dt, ...
%         'prediction_type', 'rf_classification',...
%         'output_type', 'detection',...
%         'use_probs', 0,...
%         'mask', [],...
%         'tree_mask', [], ...
%         'num_trees', [], ...
%         'max_size', 128,...
%         'incremental_results', 0);

    line_prob_g2 = predict_image(...
        'image_in', test_image,... % the mandatory arguments
        'decomposition_args', d_args_g2,...
        'predictor', predictor_g2, ...
        'prediction_type', 'rf_classification',...
        'output_type', 'detection',...
        'use_probs', 1,...
        'mask', [],...
        'tree_mask', [], ...
        'num_trees', [], ...
        'max_size', 128,...
        'incremental_results', 0);
    
%     class_predictions_dt(sample_idx,1) = line_prob_dt(bg_idx);
%     class_predictions_dt(sample_idx,2) = line_prob_dt(fg_idx);
    class_predictions_g2(sample_idx,1) = line_prob_g2(bg_idx);
    class_predictions_g2(sample_idx,2) = line_prob_g2(fg_idx);

    figure; 
%     subplot(1,2,1); imgray(line_prob_dt);
    subplot(1,2,2); imgray(line_prob_g2);
end

operating_pts = (-1:101)/100;
% [roc_dt, auc_dt] = calculate_roc_curve(class_predictions_dt(:),class_labels(:),operating_pts);
[roc_g2, auc_g2] = calculate_roc_curve(class_predictions_g2(:),class_labels(:),operating_pts);
            
%%
for i_decomp = 6%[1 3 5 6 7 9 10]
    for repeat = 1:10
        comparing_representations_experiment_csf(i_decomp, repeat, ... % the user's input
            'noise_level',  0, ...
            'num_ims',      2000, ...
            'pts_per_im',   10,...
            'make_data',    1, ...
            'do_orientation', 0, ...
            'do_detection', 0);
    end
end
%%
for noise = [0 2]
    data_dir = [exp_dir 'test\rician_' num2str(noise) '\1\'];

    %Load and reformat data
    responses_g2d = u_load([data_dir 'responses_g2d.mat']);
    responses_g2d = reshape(responses_g2d, [], 9, 4, 3);
    responses_g2d = responses_g2d(:, 5, :, :);

    %Load ground truth data
    load([data_dir 'true_labels.mat'], 'true_oris', 'true_widths');
    complex_oris = complex(cos(2*true_oris), sin(2*true_oris));

    figure; hold all;
    for i_scale = 1:5
        if i_scale < 5
            [~, predicted_oris] = ...
                gaussian_2nd_derivative_line(responses_g2d(:,:,i_scale,:));
        else
            [~, predicted_oris, predicted_scales] = ...
                gaussian_2nd_derivative_line(responses_g2d(:,:,:,:));
        end

        predicted_oris = complex(cos(2*predicted_oris), sin(2*predicted_oris));
        [~, err_stats] = ori_error(complex_oris, predicted_oris);

        plot(err_stats.abs_percentiles, 1:100);
    end
    legend({'\sigma = 1', '\sigma = 2', '\sigma = 4', '\sigma = 8', '\sigma = all'});
    title(['Orientation prediction, noise level = ', num2str(noise)]);
    
    bins = linspace(1, 8, 50);
    counts = zeros(4, 50);
    for i_scale = 1:4
        sigma = 2^(i_scale-1);
        counts(i_scale,:) = hist(true_widths(round(predicted_scales) == round(sigma)), bins);
    end
    figure; plot(bsxfun(@rdivide, counts, sum(counts))'); 
end
%%
figure; plot(true_widths, predicted_scales, 'r.');
%%
exp_dir = 'C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\';

rad_angles = pi/2 + pi*(0:5)'/6;
complex_angles = complex(cos(2*rad_angles), sin(2*rad_angles));

num_bins = 50;
width_bins = linspace(1, 8, num_bins);
ori_bins = linspace(0, pi, num_bins);
for i_decomp = {'gabor'}%{'gabor', 'gabori', 'dt'}
    for noise = [0 1 2 3]
        data_dir = [exp_dir 'test\rician_' num2str(noise) '\1\'];

        %Load and reformat data
        try
            responses_gabor = u_load([data_dir 'responses_' i_decomp{1} '.mat']);
        catch
            continue;
        end
        if strcmpi(i_decomp{1}, 'dt')
            responses_gabor = reshape(responses_gabor, [], 9, 6, 4);
            responses_gabor = permute(responses_gabor, [1 2 4 3]);
        else
            responses_gabor = reshape(responses_gabor, [], 9, 4, 6);
        end
        responses_gabor = responses_gabor(:, 5, :, :);        
        %[responses_gabor] = reorder_bands(responses_gabor, 0);

        %Load ground truth data
        load([data_dir 'true_labels.mat'], 'true_oris', 'true_widths');
        complex_oris = complex(cos(2*true_oris), sin(2*true_oris));

        if noise == 0
            true_oris = true_oris / 2;
        end
        
        %figure; a1 = gca; hold all;
        for i_scale = 5
            if i_scale < 5
                [~, predicted_bands] = ...
                    max_response_line(responses_gabor(:,:,i_scale,:));
                %figure; hold all;          
            else
                [~, predicted_bands, predicted_scales] = ...
                    max_response_line(responses_gabor(:,:,2:4,:));
            end

            %predicted_oris = complex_angles(predicted_bands);
            %[~, err_stats] = ori_error(complex_oris, predicted_oris);

            %plot(a1, err_stats.abs_percentiles, 1:100);

            counts = zeros(6, num_bins);
            for i_band = 1:6
                counts(i_band,:) = hist(true_oris(predicted_bands == i_band), ori_bins);
                
%                 if i_scale < 5
%                     [ori_centres, smoothed_abs_response] =...
%                         kernel_smoother(true_oris, squeeze(abs(responses_gabor(:,:,i_scale,i_band))) );
%                     plot(ori_centres, smoothed_abs_response);
%                 end
            end
            figure; plot(ori_bins, bsxfun(@rdivide, counts, sum(counts))', 'linewidth', 3);
            set(gca, 'ylim', [0 1], 'fontsize', 18);
            xlabel('Orientation of line', 'fontsize', 18);
            ylabel('Percentage of pixels with maximum response', 'fontsize', 18);
            title({'Which oriented sub-band produces maximum response?'; ['Noise: ' num2str(noise) ', Decomp: ' i_decomp{1}]}, 'fontsize', 18);
%             exportfig(['K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\syn_lines_' i_decomp{1} '_ori_subbands_' num2str(noise) '.pdf']);
        end

        counts = zeros(4, num_bins);
        for i_scale = 1:4
            counts(i_scale,:) = hist(true_widths(predicted_scales == i_scale), width_bins);
        end
        figure; plot(width_bins, bsxfun(@rdivide, counts, sum(counts))', 'linewidth', 3);
        set(gca, 'ylim', [0 1], 'fontsize', 18);
        xlabel('Width of line', 'fontsize', 18);
        ylabel('Percentage of pixels with maximum response', 'fontsize', 18);
        title({'Which filter scale produces maximum response?'; ['Noise: ' num2str(noise) ', Decomp: ' i_decomp{1}]}, 'fontsize', 18);
        legend({'\sigma = 1', '\sigma = 2', '\sigma = 4', '\sigma = 8'}, 'fontsize', 18);
%         exportfig(['K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\syn_lines_' i_decomp{1} '_scales_' num2str(noise) '.pdf']);
    end
end
%%
exp_dir = 'C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\';

f1 = figure; a1 = gca; hold all;
for noise = [0 1 2 3]
    data_dir = [exp_dir 'test\rician_' num2str(noise) '\1\'];

    %Load and reformat data
    responses_g2d = u_load([data_dir 'responses_g2d.mat']);
    responses_g2d = reshape(responses_g2d, [], 9, 4, 3);
    responses_g2d = responses_g2d(:, 5, :, :);

    %Load ground truth data
    load([data_dir 'true_labels.mat'], 'true_oris', 'true_widths');
    complex_oris = complex(cos(2*true_oris), sin(2*true_oris));

    figure; a2 = gca; hold all;
    i_scale = 5;        
    [~, predicted_theta, predicted_scales] = ...
        gaussian_2nd_derivative_line(responses_g2d(:,:,:,:));

    predicted_ori = complex(cos(2*predicted_theta), sin(2*predicted_theta));
    [~, err_stats] = ori_error(complex_oris, predicted_ori);

    plot(a1, err_stats.abs_percentiles, 1:100, 'linewidth', 3);
    plot(a1, [0 err_stats.abs_percentiles(75)], [75 75], 'k--', 'linewidth', 1);
    plot(a1, [err_stats.abs_percentiles(75) err_stats.abs_percentiles(75)], [0 75], 'k--', 'linewidth', 1);
    plot(a2, err_stats.abs_percentiles, 1:100, 'linewidth', 3);    
    set(a2, 'xlim', [0 90]);
    xlabel(a2, 'Absolute error in orientation prediction', 'fontsize', 18);
    ylabel(a2, 'Percentage of pixels', 'fontsize', 18);
    title(a2, 'G" cumulative distribution of orientation prediction errors', 'fontsize', 18);
    legend(a2, {['Noise level = ' num2str(noise)]}, 'location', 'southeast', 'fontsize', 18);
    exportfig([fig_root 'syn_lines_g2d_noise_cdf_' num2str(noise) '.pdf']);
end
xlabel(a1, 'Absolute error in orientation prediction', 'fontsize', 18);
ylabel(a1, 'Percentage of pixels', 'fontsize', 18);
title(a1, 'G" cumulative distribution of orientation prediction errors', 'fontsize', 18);
%legend(a1, {'Noise level = 0', 'Noise level = 1', 'Noise level = 2', 'Noise level = 3'}, 'location', 'southeast', 'fontsize', 18);
figure(f1);
exportfig([fig_root 'syn_lines_g2d_noise_cdf.pdf']);
%%
f1 = figure; a1 = gca; hold all;
for noise = [0 1 2 3]
    %Load ground truth data
    load([exp_dir 'test\rician_' num2str(noise) '\1\true_labels.mat'], 'true_oris', 'true_widths');
    complex_oris = complex(cos(2*true_oris), sin(2*true_oris));

    predicted_ori = u_load([exp_dir 'results\rician_' num2str(noise) '\1\orientation\g2d\orig\3\results.mat']);
    [~, err_stats] = ori_error(complex_oris, predicted_ori);

    plot(a1, err_stats.abs_percentiles, 1:100, 'linewidth', 3);
    plot(a1, [0 err_stats.abs_percentiles(75)], [75 75], 'k--', 'linewidth', 1);
    plot(a1, [err_stats.abs_percentiles(75) err_stats.abs_percentiles(75)], [0 75], 'k--', 'linewidth', 1);
        
end
xlabel(a1, 'Absolute error in orientation prediction', 'fontsize', 18);
ylabel(a1, 'Percentage of pixels', 'fontsize', 18);
title(a1, 'G" cumulative distribution of orientation prediction errors', 'fontsize', 18);
legend(a1, {'Noise level = 0', 'Noise level = 1', 'Noise level = 2', 'Noise level = 3'}, 'location', 'southeast', 'fontsize', 18);
figure(f1);
exportfig([fig_root 'syn_lines_g2d_reg_noise_cdf.pdf']);
%%
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];

%load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
all_gt_widths = zeros(0,1);  
centre_idx = false(0,1);
for i_ret = 21:40
    display(['Analysing image ' num2str(i_ret)]);

    %load vessel mask
    ret = rgb2gray(u_load([retroot 'images\' zerostr(i_ret,2) '_training.mat']));
    v_mask = u_load([retroot 'vessel_masks\' zerostr(i_ret,2) '_training_v_mask.mat']);
    f_mask = u_load([retroot 'fov_masks\' zerostr(i_ret,2) '_training_f_mask.mat']);
    width_map = u_load([retroot 'width_maps\' zerostr(i_ret,2) '_training_width.mat']);
    v_mask = v_mask & f_mask;
    vc_mask = bwmorph(v_mask, 'skel', 'inf');
    all_gt_widths = [all_gt_widths; width_map(v_mask)];
    centre_idx = [centre_idx; vc_mask(v_mask)];
end
mkdir([retroot 'width_maps\gt\']);
save([retroot 'width_maps\gt\all_gt_widths.mat'], 'all_gt_widths', 'centre_idx');
    
%%
load('C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\results\rician_2\1\detection\dt\orig\3\all_votes.mat')
class_labels = [true(20000,1); false(20000,1)];
operating_pts = (-1:101)/100;
auc = zeros(200,1);
for i_trees = 1:200
    predicted_lines = sum(all_votes(:,2,1:i_trees),3) / i_trees;
    [~, auc(i_trees)] =...
        calculate_roc_curve(predicted_lines,class_labels,operating_pts);
end
figure; plot(1:200, auc);
%%
noise_level = 2;
i_pts = 20000;
for i_decomp = 1:10;
    for win_size = [1 3]


        new_results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/1'...
            '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/' num2str(i_pts) '/'];
        old_results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/1'...
            '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/'];
        mkdir(new_results_dir);
        copyfile([old_results_dir 'all_votes.mat'], [new_results_dir 'all_votes.mat']);
        copyfile([old_results_dir 'results.mat'], [new_results_dir 'results.mat']);
            
    end
end
%%
for jj = 10
    comparing_representations_experiment_csf(jj, 1, ...
        'noise_level',  unixenv('NOISE_LEVEL', 0), ...
        'num_ims',      5, ...
        'num_trees',    2, ...
        'pts_per_im',   10, ...
        'dx',           unixenv('DX',64),...
        'dy',           unixenv('DY',64),...
        'cx',           unixenv('CX',32),...
        'cy',           unixenv('CY',32),...
        'make_data',    unixenv('MAKE_DATA',1), ...
        'add_edge',    unixenv('ADD_EDGE',0), ...
        'add_circles',    unixenv('ADD_CIRCLES',0), ...
        'do_orientation', unixenv('DO_ORIENTATION',0), ...
        'do_detection', unixenv('DO_DETECTION',0), ...
        'do_width', unixenv('DO_WIDTH',0), ...
        'do_tests', unixenv('DO_TESTS',[1 2 3 4 5]), ...
        'win_sizes', unixenv('WIN_SIZES', [1]) ...
    );
end

%%
dx = 64; dy = 64;
noise_level = 0;
pts_per_im = 30;

for ii = 1:10
    line_width = sample_uniform([1 8]);
    line_contrast = sample_uniform([1 8]);
    line_ori = sample_uniform([0 180]);
    line_rad = pi * line_ori / 180;

    %Generate line
    [line, label, label_centre] =...
        create_ellipse_bar(line_width/2, line_contrast, line_ori, dy, dx, dx/2, dy/2);
    signal = 1 + line;

    shuffle = randperm(dx*dy);
    r_idx = shuffle(1:pts_per_im);
    [r_rows r_cols] = ind2sub([dy dx], r_idx);
    radii = sample_uniform([1 2], [pts_per_im 1]);
    circ_contrast = sample_uniform([1 8], [pts_per_im 1]);
    for i_circ = 1:pts_per_im
        signal = signal + create_gauss_blob(...
            radii(i_circ), circ_contrast(i_circ), dy, dx, r_cols(i_circ), r_rows(i_circ));
    end
    test_image = ricernd(signal, noise_level);
    figure; imgray(test_image);
end
%%
warning('off', 'ASYM:unexpectedArgument');
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];

load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
load([retroot '/width_maps/gt/all_gt_widths.mat'], 'all_gt_widths');
num_angles = 6;
        
rad_angles = pi/2 + pi*(0:(num_angles-1))'/num_angles;
complex_angles = complex(cos(2*rad_angles), sin(2*rad_angles));

d_args.decomp_type = {'gabor'};
d_args.num_angles = 6;
d_args.sigma_range = [1 2 4 8 16];	
d_args.do_max = 0;
d_args.rotate = 0;
d_args.feature_type = 'complex';
d_args.win_size = 1;
d_args.normalise = 0;
d_args.pca = [];
D = get_samples_per_channel(d_args);
v_responses = zeros(length(gt_orientations), D);

sample_idx = 0;
for i_ret = 21:40
    display(['Analysing image ' num2str(i_ret)]);

    %load vessel mask
    ret = rgb2gray(u_load([retroot 'images\' zerostr(i_ret,2) '_training.mat']));
    v_mask = u_load([retroot 'vessel_masks\' zerostr(i_ret,2) '_training_v_mask.mat']);
    f_mask = u_load([retroot 'fov_masks\' zerostr(i_ret,2) '_training_f_mask.mat']);
    gt_ori_map = u_load([retroot 'orientations\' zerostr(i_ret,2) '_ori1.mat']);
    v_mask = v_mask & f_mask;
    sample_idx = sample_idx(end)+(1:sum(v_mask(:)));
    im_responses = compute_filter_responses(ret, d_args);
    [r_rows r_cols] = find(v_mask);
    v_responses(sample_idx, :) = sample_image_features(im_responses, r_rows(:), r_cols(:), d_args);

end
save('C:\isbe\asymmetry_project\experiments\DRIVE\gabor_orientation_prediction\v_responses.mat', 'v_responses');
%%
%%
warning('off', 'ASYM:unexpectedArgument');
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];

load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations');       

d_args.decomp_type = {'dt'};
d_args.levels = 1:5;
d_args.feature_shape = 'rect';
d_args.feature_type = 'complex';
d_args.do_max = 0;
d_args.rotate = 0;
d_args.use_nag = 0;
d_args.win_size = 1;
d_args.normalise = 0;
d_args.pca = [];
d_args.rgb_channel = 'rgb';
D = get_samples_per_channel(d_args);
v_responses = zeros(length(gt_orientations), D);
clear gt_orientations;

sample_idx = 0;
for i_ret = 21:40
    display(['Analysing image ' num2str(i_ret)]);

    %load vessel mask
    ret = rgb2gray(u_load([retroot 'images\' zerostr(i_ret,2) '_training.mat']));
    v_mask = u_load([retroot 'vessel_masks\' zerostr(i_ret,2) '_training_v_mask.mat']);
    f_mask = u_load([retroot 'fov_masks\' zerostr(i_ret,2) '_training_f_mask.mat']);
    gt_ori_map = u_load([retroot 'orientations\' zerostr(i_ret,2) '_ori1.mat']);
    v_mask = v_mask & f_mask;
    sample_idx = sample_idx(end)+(1:sum(v_mask(:)));
    im_responses = compute_filter_responses(ret, d_args);
    [r_rows r_cols] = find(v_mask);
    v_responses(sample_idx, :) = sample_image_features(im_responses, r_rows(:), r_cols(:), d_args);

end
save('C:\isbe\asymmetry_project\experiments\DRIVE\dt_orientation_prediction\v_responses.mat', 'v_responses');
%%
warning('off', 'ASYM:unexpectedArgument');
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];

load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
load([retroot '/width_maps/gt/all_gt_widths.mat'], 'all_gt_widths');

d_args.decomp_type = {'g2d'};
d_args.sigma_range = [1 2 4 8 16];
d_args.win_size = 1;
d_args.normalise = 0;
d_args.pca = [];
D = get_samples_per_channel(d_args);
v_responses = zeros(length(gt_orientations), D);

sample_idx = 0;
for i_ret = 21:40
    display(['Analysing image ' num2str(i_ret)]);

    %load vessel mask
    ret = rgb2gray(u_load([retroot 'images\' zerostr(i_ret,2) '_training.mat']));
    v_mask = u_load([retroot 'vessel_masks\' zerostr(i_ret,2) '_training_v_mask.mat']);
    f_mask = u_load([retroot 'fov_masks\' zerostr(i_ret,2) '_training_f_mask.mat']);
    gt_ori_map = u_load([retroot 'orientations\' zerostr(i_ret,2) '_ori1.mat']);
    v_mask = v_mask & f_mask;
    sample_idx = sample_idx(end)+(1:sum(v_mask(:)));
    im_responses = compute_filter_responses(ret, d_args);
    [r_rows r_cols] = find(v_mask);
    v_responses(sample_idx, :) = sample_image_features(im_responses, r_rows(:), r_cols(:), d_args);

end
mkdir('C:\isbe\asymmetry_project\experiments\DRIVE\g2d_orientation_prediction');
save('C:\isbe\asymmetry_project\experiments\DRIVE\g2d_orientation_prediction\v_responses.mat', 'v_responses');
%%
warning('off', 'ASYM:unexpectedArgument');
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];

load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
load([retroot '/width_maps/gt/all_gt_widths.mat'], 'all_gt_widths');
num_angles = 6;
        
rad_angles = pi/2 + pi*(0:(num_angles-1))'/num_angles;
complex_angles = complex(cos(2*rad_angles), sin(2*rad_angles));

d_args.decomp_type = {'dt'};
d_args.levels = 1:5;	
d_args.do_max = 0;
d_args.rotate = 0;
d_args.feature_shape = 'rect';
d_args.feature_type = 'complex';
d_args.win_size = 1;
d_args.normalise = 0;
d_args.pca = [];
d_args.use_nag = 0;
D = get_samples_per_channel(d_args);
v_responses = zeros(length(gt_orientations), D);

sample_idx = 0;
for i_ret = 21:40
    display(['Analysing image ' num2str(i_ret)]);

    %load vessel mask
    ret = rgb2gray(u_load([retroot 'images\' zerostr(i_ret,2) '_training.mat']));
    v_mask = u_load([retroot 'vessel_masks\' zerostr(i_ret,2) '_training_v_mask.mat']);
    f_mask = u_load([retroot 'fov_masks\' zerostr(i_ret,2) '_training_f_mask.mat']);
    gt_ori_map = u_load([retroot 'orientations\' zerostr(i_ret,2) '_ori1.mat']);
    v_mask = v_mask & f_mask;
    sample_idx = sample_idx(end)+(1:sum(v_mask(:)));
    im_responses = compute_filter_responses(ret, d_args);
    [r_rows r_cols] = find(v_mask);
    v_responses(sample_idx, :) = sample_image_features(im_responses, r_rows(:), r_cols(:), d_args);

end
save('C:\isbe\asymmetry_project\experiments\DRIVE\dt_orientation_prediction\v_responses.mat', 'v_responses');
%%
%Load and reformat data
clear; pack;
fig_root = 'K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\';
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];
load('C:\isbe\asymmetry_project\experiments\DRIVE\gabor_orientation_prediction\v_responses.mat', 'v_responses');
load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations');
sample_idx = randperm(length(gt_orientations));
sample_idx = sample_idx(1:2e4);

v_responses = reshape(v_responses(sample_idx,:), [], 1, 5, 6);
gt_orientations = gt_orientations(sample_idx,:);

true_oris = angle(gt_orientations) / 2;
clear gt_orientations;
pack;

num_bins = 50;
ori_bins = linspace(-pi/2, pi/2, num_bins);

%figure; a1 = gca; hold all;
for i_scale = 6%1:5
    if i_scale < 6
        [~, predicted_bands] = ...
            max_response_line(v_responses(:,:,i_scale,:));
        scale_str = num2str(i_scale);
    else
        [~, predicted_bands, predicted_scales] = ...
            max_response_line(v_responses(:,:,1:end,:));
        scale_str = 'all';
    end

    counts = zeros(6, num_bins);
    for i_band = 1:6
        counts(i_band,:) = hist(true_oris(predicted_bands == i_band), ori_bins);
    end
    figure; plot(ori_bins, bsxfun(@rdivide, counts, sum(counts))', 'linewidth', 3);
    set(gca, 'ylim', [0 1], 'xlim', [-pi/2 pi/2], 'fontsize', 18);
    xlabel('Orientation of vessel', 'fontsize', 18);
    ylabel('Percentage of pixels with maximum response', 'fontsize', 18);
    title({'Which oriented sub-band produces maximum response?';}, 'fontsize', 18);
    exportfig([fig_root 'ret_vessels_gabor_ori_subbands.pdf']);
end
clear v_responses;
%%
load([retroot '/width_maps/gt/all_gt_widths.mat'], 'all_gt_widths');
all_gt_widths = all_gt_widths(sample_idx,:);
counts = zeros(5, 16);
for i_scale = 1:5
    counts(i_scale,:) = hist(all_gt_widths(predicted_scales == i_scale), 1:16);
end
figure; plot(1:16, bsxfun(@rdivide, counts, sum(counts))', 'linewidth', 3); 
xlabel('Width of vessel', 'fontsize', 18);
set(gca, 'ylim', [0 1], 'fontsize', 18);
ylabel('Percentage of pixels with maximum response', 'fontsize', 18);
title({'Which filter scale produces maximum response?'}, 'fontsize', 18);
legend({'\sigma = 1', '\sigma = 2', '\sigma = 4', '\sigma = 8', '\sigma = 16'}, 'fontsize', 18);
exportfig([fig_root 'ret_vessels_gabor_scales.pdf']);
%%
%Load and reformat data
clear; pack;
fig_root = 'K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\';
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];
load('C:\isbe\asymmetry_project\experiments\DRIVE\gabor_orientation_prediction\v_responses.mat', 'v_responses');
load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations');
load([retroot '/width_maps/gt/all_gt_widths.mat'], 'all_gt_widths');
sample_idx = all_gt_widths < 4;

v_responses = reshape(v_responses(sample_idx,:), [], 1, 5, 6);
gt_orientations = gt_orientations(sample_idx,:);

true_oris = angle(gt_orientations) / 2;
clear gt_orientations;
pack;

num_bins = 50;
ori_bins = linspace(-pi/2, pi/2, num_bins);

%figure; a1 = gca; hold all;
for i_scale = 6%1:5
    if i_scale < 6
        [~, predicted_bands] = ...
            max_response_line(v_responses(:,:,i_scale,:));
        scale_str = num2str(i_scale);
    else
        [~, predicted_bands, predicted_scales] = ...
            max_response_line(v_responses(:,:,1:end,:));
        scale_str = 'all';
    end

    counts = zeros(6, num_bins);
    for i_band = 1:6
        counts(i_band,:) = hist(true_oris(predicted_bands == i_band), ori_bins);
    end
    figure; plot(ori_bins, bsxfun(@rdivide, counts, sum(counts))', 'linewidth', 3);
    set(gca, 'ylim', [0 1], 'xlim', [-pi/2 pi/2], 'fontsize', 18);
    xlabel('Orientation of vessel', 'fontsize', 18);
    ylabel('Percentage of pixels with maximum response', 'fontsize', 18);
    title({'Which oriented sub-band produces maximum response?';}, 'fontsize', 18);
    exportfig([fig_root 'ret_thin_vessels_gabor_ori_subbands.pdf']);
end
clear v_responses;
%
%%
%Load and reformat data
clear; pack;
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];
load('C:\isbe\asymmetry_project\experiments\DRIVE\dt_orientation_prediction\v_responses.mat', 'v_responses');
load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations');

sample_idx = randperm(length(gt_orientations));
sample_idx = sample_idx(1:2e4);

v_responses = reshape(v_responses(sample_idx,:), [], 1, 6, 5);
v_responses = permute(v_responses, [1 2 4 3]);
gt_orientations = gt_orientations(sample_idx,:);
true_oris = angle(gt_orientations) / 2;
clear gt_orientations;

pack;

num_bins = 50;
ori_bins = linspace(-pi/2, pi/2, num_bins);

%figure; a1 = gca; hold all;
for i_scale = 6%1:5
    if i_scale < 6
        [~, predicted_bands] = ...
            max_response_line(v_responses(:,:,:,i_scale));
        scale_str = num2str(i_scale);
    else
        [~, predicted_bands, predicted_scales] = ...
            max_response_line(v_responses(:,:,:,:));
        scale_str = 'all';
    end

    counts = zeros(6, num_bins);
    for i_band = 1:6
        counts(i_band,:) = hist(true_oris(predicted_bands == i_band), ori_bins);
    end
    figure; plot(ori_bins, bsxfun(@rdivide, counts, sum(counts))');
    xlabel('Orientation of vessel');
    ylabel('Percentage of pixels with maximum response');
    title({'Which oriented sub-band produces maximum response?'; ['Scale: ' scale_str]});
end
clear v_responses;
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
title({'Which filter scale produces maximum response?'});
legend({'\sigma = 1', '\sigma = 2', '\sigma = 4', '\sigma = 8', '\sigma = 16'});
%%
%Gaussian prediction of REtinograms - all vessels
clear; pack;
fig_root = 'K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\';
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];
load('C:\isbe\asymmetry_project\experiments\DRIVE\g2d_orientation_prediction\v_responses.mat', 'v_responses');
load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
load([retroot '/width_maps/gt/all_gt_widths.mat'], 'all_gt_widths');

v_responses = reshape(v_responses, [], 1, 5, 3);
true_oris = angle(gt_orientations) / 2;
%
figure; hold all; a1 = gca;
figure; hold all; a2 = gca;
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
    if i_scale == 6
        plot(a2, err_stats.abs_percentiles, 1:100, 'linewidth', 3);
    end
end
xlabel(a1, 'Absolute error in orientation preditcion', 'fontsize', 18);
ylabel(a1, 'Percentage of pixels', 'fontsize', 18);
title(a1, {'G" cumulative distribution of orientation error'}, 'fontsize', 18);
xlabel(a2, 'Absolute error in orientation preditcion', 'fontsize', 18);
ylabel(a2, 'Percentage of pixels', 'fontsize', 18);
title(a2, {'G" cumulative distribution of orientation error'}, 'fontsize', 18);
legend(a1, {'All scales', '\sigma = 1', '\sigma = 2', '\sigma = 4', '\sigma = 8', '\sigma = 16'}, 'location', 'southeast');
axes(a1); exportfig([fig_root 'ret_vessels_g2d_orientation_prediction_separate.pdf']);
legend(a2, {'All scales'}, 'location', 'southeast');
axes(a2); exportfig([fig_root 'ret_vessels_g2d_orientation_prediction_all.pdf']);
%%
%Gabor prediction RF
retroot = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\training\';
pred_dir = [retroot 'predictions\orientation\rf_regression\'];
label_dir = [retroot 'orientations\'];
fov_mask_dir = [retroot 'fov_masks\'];
fg_mask_dir = [retroot 'vessel_masks\'];

rf_codes( 6,:) = {'24069'};
rf_codes(1,:) = {'30469'};
rf_codes(2,:) = {'30470'};
rf_codes(3,:) = {'30471'};
rf_codes(4,:) = {'30472'};
rf_codes(5,:) = {'30473'};
figure; hold all; a1 = gca;
figure; hold all; a2 = gca;
for i_scale = [6 1:5]
    [~, ~, err_stats] =...
        compute_image_orientation_errors([pred_dir rf_codes{i_scale,1} '\'], fg_mask_dir,...
        'gt_orientations', gt_orientations, 'fov_mask_dir', fov_mask_dir);
    plot(a1, err_stats.abs_percentiles, 1:100, 'linewidth', 3);
    if i_scale == 6
        plot(a2, err_stats.abs_percentiles, 1:100, 'linewidth', 3);
    end
end

xlabel(a1, 'Absolute error in orientation preditcion', 'fontsize', 18);
ylabel(a1, 'Percentage of pixels', 'fontsize', 18);
title(a1, {'G" cumulative distribution of orientation error'}, 'fontsize', 18);
xlabel(a2, 'Absolute error in orientation preditcion', 'fontsize', 18);
ylabel(a2, 'Percentage of pixels', 'fontsize', 18);
title(a2, {'G" cumulative distribution of orientation error'}, 'fontsize', 18);
legend(a1, {'All scales', '\sigma = 1', '\sigma = 2', '\sigma = 4', '\sigma = 8', '\sigma = 16'}, 'location', 'southeast');
axes(a1); exportfig([fig_root 'ret_vessels_gabor_RF_orientation_prediction_separate.pdf']);
legend(a2, {'All scales'}, 'location', 'southeast');
axes(a2); exportfig([fig_root 'ret_vessels_gabor_RF_orientation_prediction_all.pdf']);
%%
%Gaussian prediction of REtinograms - thin vessels
clear; pack;
fig_root = 'K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\';
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];
load('C:\isbe\asymmetry_project\experiments\DRIVE\g2d_orientation_prediction\v_responses.mat', 'v_responses');
load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
load([retroot '/width_maps/gt/all_gt_widths.mat'], 'all_gt_widths');

sample_idx = all_gt_widths < 2;

v_responses = reshape(v_responses(sample_idx,:), [], 1, 5, 3);
gt_orientations = gt_orientations(sample_idx,:);
true_oris = angle(gt_orientations) / 2;
%
figure; hold all; a1 = gca;
figure; hold all; a2 = gca;
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
    if i_scale == 6
        plot(a2, err_stats.abs_percentiles, 1:100, 'linewidth', 3);
    end
end
xlabel(a1, 'Absolute rror in orientation preidtcion', 'fontsize', 18);
ylabel(a1, 'Percentage of pixels', 'fontsize', 18);
title(a1, {'G" cumulative distribution of orientation error'}, 'fontsize', 18);
xlabel(a2, 'Absolute rror in orientation preidtcion', 'fontsize', 18);
ylabel(a2, 'Percentage of pixels', 'fontsize', 18);
title(a2, {'G" cumulative distribution of orientation error'}, 'fontsize', 18);
legend(a1, {'\sigma = 1', '\sigma = 2', '\sigma = 4', '\sigma = 8', '\sigma = 16', 'All scales'}, 'location', 'southeast');
axes(a1); exportfig([fig_root 'ret_thin_vessels_g2d_orientation_prediction_separate.pdf']);
legend(a2, {'All scales'}, 'location', 'southeast');
axes(a2); exportfig([fig_root 'ret_thin_vessels_g2d_orientation_prediction_all.pdf']);
%%
[orientation_errs, err_stats] = ori_error(predicted_ori, true_oris);
idx = ~isnan(orientation_errs);
[centres,smoothed_y] = kernel_smoother(all_gt_widths(idx), 180*abs(orientation_errs(idx))/pi);


figure; plot(centres, smoothed_y);
%%
load([exp_dir 'test/rician_' num2str(noise_level) '/1/true_labels.mat'], '*idx');                       
noise_level = 2;
for i_decomp = 1:10;
    for win_size = [1 3]
        figure; hold all;
        title(['Decomp: ' d_args{i_decomp}.decomp_type{1} ', w = ' num2str(win_size)]);
        leg_text = cell(0,1);
        for i_pts = [10000:2000:20000]% 25000:5000:40000]
            auc = zeros(200,1);
            results_dir = [exp_dir '/results/rician_' num2str(noise_level) '/1'...
                '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/' num2str(i_pts) '/'];
            
            try
                load([results_dir 'all_votes.mat'], 'all_votes'); 
            catch
                continue;
            end
            for i_trees = 5:5:100

                predicted_lines = sum(all_votes(:,2,1:i_trees),3) / i_trees;
                [~, auc(i_trees)] =...
                    calculate_roc_curve(predicted_lines([centre_idx; edge_idx]),class_labels([centre_idx; edge_idx]),operating_pts);
            end
            plot(5:5:100, auc(5:5:100));
            leg_text{end+1} = ['N =' num2str(i_pts)];
        end
        legend(leg_text, 'location', 'southeast');
        set(gca, 'ylim', [0.9 0.99]);
    end
end
%%
exp_dir2 = [asymmetryroot 'experiments/synthetic_lines/line_on_edge/'];
exp_dir1 = [asymmetryroot 'experiments/synthetic_lines/comparing_representations/']; 

load([exp_dir 'test/rician_' num2str(noise_level) '/1/true_labels.mat'], '*idx');                       
noise_level = 2;
for i_decomp = 1:10;
    for win_size = [1 3]
        i_pts = 20000;
        results_dir1 = [exp_dir1 '/results/rician_' num2str(noise_level) '/1'...
                '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/' num2str(i_pts) '/'];
        load([results_dir1 'all_votes.mat'], 'all_votes');
        predicted_lines1 = sum(all_votes(:,2,:),3) / 200;
        
        results_dir2 = [exp_dir2 '/results/rician_' num2str(noise_level) '/1'...
                '/detection/' d_args{i_decomp}.decomp_type{1} '/orig/' num2str(win_size) '/' num2str(i_pts) '/'];
        load([results_dir2 'all_votes.mat'], 'all_votes');
        predicted_lines2 = sum(all_votes(:,2,:),3) / 100;
        
        [~, auc_orig] = calculate_roc_curve(...
            predicted_lines1,...
            class_labels, operating_pts);
        [~, auc_all] = calculate_roc_curve(...
            predicted_lines2,...
            class_labels, operating_pts); 
        [~, auc_roi] = calculate_roc_curve(...
            predicted_lines2([true(20000,1); edge_idx]),...
            class_labels([true(20000,1); edge_idx]), operating_pts); 
        display([results_dir2 ':']);
        display(['Az (orig) = ' num2str(auc_orig,4) ', Az (edge) = ' num2str(auc_all,4) ', Az (roi) = ' num2str(auc_roi,4)]); 
    end
end
%%
%Sample properties of line
line_width = 4;
line_contrast = 4;
line_ori = 27;

for noise_level = 0:3

    %Generate line
    [line, label, label_centre] =...
        create_ellipse_bar(line_width/2, line_contrast, line_ori, 64, 64, 32, 32);
    
    test_image = ricernd(1 + line, noise_level);
    figure; imgray(test_image);
    write_im_from_colormap(test_image, [fig_root 'syn_line_rician_' num2str(noise_level) '.png']);
end
%%
for ii = 1:4
    line_width = sample_uniform([1 8]);
    line_contrast = sample_uniform([1 8]);
    line_ori = sample_uniform([0 360]);
    [line, label, label_centre] =...
        create_ellipse_bar(line_width/2, line_contrast, line_ori, 64, 64, 32, 32);
    write_im_from_colormap(line, [fig_root 'syn_line_rician_i' num2str(ii) '.png'], gray(256), [0 8]);
end
%%
[filters] = gabor_filters(6, 16);
for i_filter = 1:6
    %figure; imgray(real(filters(:,:,i_filter)));
    write_im_from_colormap(real(filters(33:end-32,33:end-32,i_filter)),...
        [fig_root 'gabor_filter_' num2str(i_filter) '.png']);
end
%%
retroot = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\';
fov_mask_dir = [retroot 'fov_masks\'];
fg_mask_dir = [retroot 'vessel_masks\'];

pred_dir = [retroot 'predictions\detection\rf_classification\'];
label_dir = [retroot 'vessel_masks\'];
    
tp_counts = zeros(20, 102);
fp_counts = zeros(20, 102);
t_counts = zeros(20, 1);
f_counts = zeros(20, 1);
aucci = zeros(4,2);

f1 = figure('windowstyle', 'normal');
graph(f1);
set(gca,'box','on'); 
axis equal; axis([0 1 0 1]); hold all;
%title('ROC curves for vessel segmentation, DRIVE database');
xlabel('FPF');
ylabel('TPF');
legend_label = cell(4,1);
colors = lines(4);
colors = colors([4 2 1 3],:);
%
for data_type = 1:5
    for ii = 1:20

        switch data_type

            case 1
                label = 'Forest Gabor';
                vessel_prob = load_uint8([pred_dir '24303\' zerostr(ii, 2) '_test_pred.mat']);

            case 2
                label = 'Forest DT-CWT';
                vessel_prob = load_uint8([pred_dir '9162\' zerostr(ii, 2) '_test_pred.mat']);
            case 3
                label = 'Staal';
                vessel_prob = double(rgb2gray(imread(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\predictions\staal\' zerostr(ii-1, 2) '.bmp'])))/255;
            case 4
                label = 'Niemeijer';
                vessel_prob = double(rgb2gray(imread(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\predictions\niemeijer\' zerostr(ii, 2) '_PC_soft.bmp'])))/255;
                vessel_prob(1,:) = [];

        end

        v_mask = u_load([fg_mask_dir zerostr(ii,2) '_test_v_mask.mat']);
        f_mask = u_load([fov_mask_dir zerostr(ii,2) '_test_f_mask.mat']);

        %Compute ROC counts for image
        [~, ~, tp_count fp_count] = calculate_roc_image(vessel_prob, v_mask,(-1:100)/100, f_mask, 'dilate', 0);
        t_counts(ii) = sum(v_mask(f_mask));
        f_counts(ii) = sum(~v_mask(f_mask));

        %Increment total counts
        tp_counts(ii,:) = tp_count;
        fp_counts(ii,:) = fp_count;


    end

    aucci(data_type,:) = bootci(2000,@compute_auc, fp_counts, tp_counts, f_counts, t_counts);
        
    %Compute ROC points for complete set of data
    roc_pts = [sum(fp_counts)' / sum(f_counts) sum(tp_counts)' / sum(t_counts)];

    %Compute AUC for ROC curve
    auc = sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* roc_pts(1:end-1,2)) + ...
            0.5*sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* (roc_pts(2:end,2)-roc_pts(1:end-1,2)) );
        
    figure(f1);
    plot(roc_pts(:,1), roc_pts(:,2), '-', 'color', colors(data_type,:), 'linewidth', 1); axis([0 1 0 1]);
    legend_label{data_type,1} = [ label ' AUC: ' num2str(auc)];%label;%
end
% legend(legend_label, 'location', 'southeast');
% exportfig('K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\segmentation_roc.png', 'png');
%exportfig('K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\segmentation_roc.pdf');
%%
pred_ga1 = u_load('C:\isbe\asymmetry_project\experiments\DRIVE\comparing_representations\results\1\orientation\gabor\orig\3\results.mat');
pred_ga2 = u_load('C:\isbe\asymmetry_project\experiments\DRIVE\comparing_representations\results\2\orientation\gabor\orig\3\results.mat');

pred_dt1 = u_load('C:\isbe\asymmetry_project\experiments\DRIVE\comparing_representations\results\1\orientation\dt\orig\3\results.mat');
pred_dt2 = u_load('C:\isbe\asymmetry_project\experiments\DRIVE\comparing_representations\results\2\orientation\dt\orig\3\results.mat');

gt1 = load('C:\isbe\asymmetry_project\experiments\DRIVE\comparing_representations\test\1\true_labels.mat', 'true_oris_test');
gt2 = load('C:\isbe\asymmetry_project\experiments\DRIVE\comparing_representations\test\2\true_labels.mat', 'true_oris_test');

ori_ga1 = ori_error(pred_ga1, gt1.true_oris_test);
ori_ga2 = ori_error(pred_ga2, gt2.true_oris_test);

ori_dt1 = ori_error(pred_dt1, gt1.true_oris_test);
ori_dt2 = ori_error(pred_dt2, gt2.true_oris_test);

figure; 
subplot(2,3,1); plot(abs(ori_dt1), abs(ori_dt2), 'r.'); axis equal; title('DT_1 v DT_2');
subplot(2,3,2); plot(abs(ori_dt1), abs(ori_ga1), 'r.'); axis equal; title('DT_1 v Gabor_1');
subplot(2,3,3); plot(abs(ori_dt1), abs(ori_ga2), 'r.'); axis equal; title('DT_1 v Gabor_2');
subplot(2,3,4); plot(abs(ori_ga1), abs(ori_ga2), 'r.'); axis equal; title('DT_1 v DT_2');
subplot(2,3,5); plot(abs(ori_ga1), abs(ori_dt2), 'r.'); axis equal; title('DT_1 v DT_2');
subplot(2,3,6); plot(abs(ori_dt2), abs(ori_ga2), 'r.'); axis equal; title('DT_1 v DT_2');
%%
retroot = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\training\';
fov_mask_dir = [retroot 'fov_masks\'];
fg_mask_dir = [retroot 'vessel_masks\'];
pred_dir = [retroot 'predictions\orientation\rf_regression\'];
label_dir = [retroot 'orientations\'];

load([label_dir 'gt\all_gt_orientations.mat'], 'gt_orientations');
load([retroot 'width_maps\gt\all_gt_widths.mat'], 'all_gt_widths');
[uni_widths,~, width_idx] = unique(all_gt_widths);

rf_codes(1,:) = {'30469', 'gabor_1', '3', 'orig'};
rf_codes(2,:) = {'30470', 'gabor_2', '3', 'orig'};
rf_codes(3,:) = {'30471', 'gabor_4', '3', 'orig'};
rf_codes(4,:) = {'30472', 'gabor_8', '3', 'orig'};
rf_codes(5,:) = {'30473', 'gabor_16', '3', 'orig'};
rf_codes(6,:) = {'24069', 'gabor', '3', 'orig'};

figure; hold all; a1 = gca;
for ii = [6 1:5]

    [orientation_errors] =...
        compute_image_orientation_errors([pred_dir rf_codes{ii,1} '\'], fg_mask_dir,...
        'gt_orientations', gt_orientations, 'fov_mask_dir', fov_mask_dir);
    valid_idx = ~isnan(orientation_errors);
    [width_centres, smoothed_abs_errs] =...
        kernel_smoother(all_gt_widths(valid_idx), 180*abs(orientation_errors(valid_idx))/pi );
    plot(a1, width_centres, smoothed_abs_errs, 'linewidth', 2);
    
    %figure; boxplot(abs(orientation_errors), width_idx);
end
set(gca, 'xlim', [1 10], 'ylim', [0 22.5])
legend({'all scales', '\sigma=1', '\sigma=2', '\sigma=4', '\sigma=8', '\sigma=16'});
ylabel(a1, 'Absolute error (degrees)', 'fontsize', 18);
xlabel(a1, 'Vessel width (pixels)', 'fontsize', 18);
title(a1, {'Estimate of Mean Abolsute Error as a function of vessel width'}, 'fontsize', 18);
exportfig('K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\g2d_rf_ori_errors_v_width.pdf');
%
clear; pack;
fig_root = 'K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\';
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
set(gca, 'xlim', [1 10], 'ylim', [0 45])
legend({'all scales', '\sigma=1', '\sigma=2', '\sigma=4', '\sigma=8', '\sigma=16'});
ylabel(a1, 'Absolute error (degrees)', 'fontsize', 18);
xlabel(a1, 'Vessel width (pixels)', 'fontsize', 18);
title(a1, {'Estimate of Mean Abolsute Error as a function of vessel width'}, 'fontsize', 18);
exportfig('K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\g2d_ori_errors_v_width.pdf');
%%
retroot = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\training\';
fov_mask_dir = [retroot 'fov_masks\'];
fg_mask_dir = [retroot 'vessel_masks\'];
pred_dir = [retroot 'predictions\orientation\rf_regression\'];
label_dir = [retroot 'orientations\'];

load([label_dir 'gt\all_gt_orientations.mat'], 'gt_orientations');
load([retroot 'width_maps\gt\all_gt_widths.mat'], 'all_gt_widths');
[uni_widths,~, width_idx] = unique(all_gt_widths);
num_widths = length(uni_widths);
rf_codes(1,:) = {'30469', 'gabor_1', '3', 'orig'};
rf_codes(2,:) = {'30470', 'gabor_2', '3', 'orig'};
rf_codes(3,:) = {'30471', 'gabor_4', '3', 'orig'};
rf_codes(4,:) = {'30472', 'gabor_8', '3', 'orig'};
rf_codes(5,:) = {'30473', 'gabor_16', '3', 'orig'};
rf_codes(6,:) = {'24069', 'gabor', '3', 'orig'};

figure; hold all; a1 = gca;
for ii = [6 1:5]

    [orientation_errors] =...
        compute_image_orientation_errors([pred_dir rf_codes{ii,1} '\'], fg_mask_dir,...
        'gt_orientations', gt_orientations, 'fov_mask_dir', fov_mask_dir);
    valid_idx = ~isnan(orientation_errors);
    width_meds = zeros(num_widths, 1);
    for i_width = 1:num_widths
        width_meds(i_width) = ...
            naNmedian(180*abs(orientation_errors(width_idx==i_width))/pi);
    end
    plot(a1, uni_widths, width_meds, 'linewidth', 2);
    
end
set(gca, 'xlim', [1 10], 'ylim', [0 22.5])
legend({'all scales', '\sigma=1', '\sigma=2', '\sigma=4', '\sigma=8', '\sigma=16'});
ylabel(a1, 'Absolute error (degrees)', 'fontsize', 18);
xlabel(a1, 'Vessel width (pixels)', 'fontsize', 18);
title(a1, {'Median Abolsute Error as a function of vessel width'}, 'fontsize', 18);
%%
exportfig('K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\g2d_rf_ori_meds_v_width.pdf');
%
clear; pack;
fig_root = 'K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\';
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];
load('C:\isbe\asymmetry_project\experiments\DRIVE\g2d_orientation_prediction\v_responses.mat', 'v_responses');
load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
load([retroot '/width_maps/gt/all_gt_widths.mat'], 'all_gt_widths');
[uni_widths,~, width_idx] = unique(all_gt_widths);
num_widths = length(uni_widths);

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
    width_meds = zeros(num_widths, 1);
    for i_width = 1:num_widths
        width_meds(i_width) = ...
            naNmedian(180*abs(orientation_errors(width_idx==i_width))/pi);
    end
    plot(a1, uni_widths, width_meds, 'linewidth', 2);
end
set(gca, 'xlim', [1 10], 'ylim', [0 45])
legend({'all scales', '\sigma=1', '\sigma=2', '\sigma=4', '\sigma=8', '\sigma=16'});
ylabel(a1, 'Absolute error (degrees)', 'fontsize', 18);
xlabel(a1, 'Vessel width (pixels)', 'fontsize', 18);
title(a1, {'Median Abolsute Error as a function of vessel width'}, 'fontsize', 18);
exportfig('K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\g2d_ori_meds_v_width.pdf');





