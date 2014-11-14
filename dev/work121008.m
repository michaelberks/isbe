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
    for noise = 0%[0 1 2 3]
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
        for i_scale = [3 5]
            if i_scale < 5
                [~, predicted_bands] = ...
                    max_response_line(responses_gabor(:,:,i_scale,:));
                %figure; hold all;          
            else
                [~, predicted_bands, predicted_scales] = ...
                    max_response_line(responses_gabor(:,:,:,:));
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
            exportfig(['K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\syn_lines_' i_decomp{1} '_ori_subbands_' num2str(noise) '.pdf']);
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
        exportfig(['K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\syn_lines_' i_decomp{1} '_scales_' num2str(noise) '.pdf']);
    end
end
%%
exp_dir = 'C:\isbe\asymmetry_project\experiments\synthetic_lines\comparing_representations\';

num_bins = 50;
width_bins = linspace(1, 8, num_bins);

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

    %figure; a1 = gca; hold all;
    for i_scale = 5
        if i_scale < 5
            [~, predicted_theta] = ...
                gaussian_2nd_derivative_line(responses_gabor(:,:,i_scale,:));
        else
            [~, predicted_theta, predicted_scales] = ...
                gaussian_2nd_derivative_line(responses_g2d(:,:,:,:));
        end

        predicted_ori = complex(cos(2*predicted_theta), sin(2*predicted_theta));
        [~, err_stats] = ori_error(complex_oris, predicted_ori);

        plot(a1, err_stats.abs_percentiles, 1:100);
        
    end
    
    counts = zeros(4, num_bins);
    for i_scale = 1:4
        counts(i_scale,:) = hist(true_widths(predicted_scales == round(2^(i_scale-1))), width_bins);
    end
    figure; plot(width_bins, bsxfun(@rdivide, counts, sum(counts))'); 
    xlabel('Width of line');
    ylabel('Percentage of pixels with maximum response');
    title({'Which filter scale produces maximum response?'; ['Noise: ' num2str(noise)]});
    legend({'\sigma = 1', '\sigma = 2', '\sigma = 4', '\sigma = 8'});
    exportfig(['K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\syn_lines_g2d_scales_' num2str(noise) '.pdf']);
end
xlabel(a1, 'Absolute orientation error');
title(a1, {'CDF of orientation prediction errors, G" filters'; ['Noise: ' num2str(noise)]});
legend(a1, {'Noise level = 0', 'Noise level = 1', 'Noise level = 2', 'Noise level = 3'}, 'location', 'southeast');
figure(f1);
exportfig(['K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\syn_lines_g2d_noise_cdf.pdf']);
%%
f1 = figure; a1 = gca; hold all;
for noise = [0 1 2 3]
    %Load ground truth data
    load([exp_dir 'test\rician_' num2str(noise) '\1\true_labels.mat'], 'true_oris', 'true_widths');
    complex_oris = complex(cos(2*true_oris), sin(2*true_oris));

    predicted_ori = u_load([exp_dir 'results\rician_' num2str(noise) '\1\orientation\g2d\orig\3\results.mat']);
    [~, err_stats] = ori_error(complex_oris, predicted_ori);

    plot(a1, err_stats.abs_percentiles, 1:100);
        
end
xlabel(a1, 'Absolute orientation error');
title(a1, {'CDF of orientation prediction errors, G" filters'; ['Noise: ' num2str(noise)]});
legend(a1, {'Noise level = 0', 'Noise level = 1', 'Noise level = 2', 'Noise level = 3'}, 'location', 'southeast');
figure(f1);
exportfig(['K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\syn_lines_g2d_reg_noise_cdf.pdf']);
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
        'make_data',    0, ...
        'add_edge',    0, ...
        'add_circles',    1, ...
        'do_orientation', 0, ...
        'do_detection', 1, ...
        'do_width', 0, ...
        'do_tests', 0, ...
        'win_sizes', 1 ...
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
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];
load('C:\isbe\asymmetry_project\experiments\DRIVE\dt_orientation_prediction\v_responses.mat', 'v_responses');
load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations');

v_responses = reshape(v_responses, [], 1, 6, 5);
true_oris = angle(gt_orientations) / 2;
clear gt_orientations;
pack;

num_bins = 50;
ori_bins = linspace(-pi/2, pi/2, num_bins);

%figure; a1 = gca; hold all;
for i_scale = 6%1:5
    if i_scale < 6
        [~, predicted_bands] = ...
            max_response_line(permute(v_responses(:,:,:,i_scale), [1 2 4 3]));
        scale_str = num2str(i_scale);
    else
        [~, predicted_bands, predicted_scales] = ...
            max_response_line(permute(v_responses(:,:,:,:), [1 2 4 3]));
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
[~, predicted_bands, predicted_scales] = ...
    max_response_line(v_responses(:,:,:,:));
for i_scale = 1:5
    counts = zeros(6, num_bins);
    for i_band = 1:6
        counts(i_band,:) = hist(true_oris((predicted_bands == i_band) & (predicted_scales == i_scale)), ori_bins);
    end
    figure; plot(ori_bins, bsxfun(@rdivide, counts, sum(counts))'); 
    xlabel('Orientation of vessel');
    ylabel('Percentage of pixels with maximum response');
    title({'Which oriented sub-band produces maximum response?'; ['Scale: ' num2str(i_scale)]});
end
%%
clear; pack;
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];
load('C:\isbe\asymmetry_project\experiments\DRIVE\g2d_orientation_prediction\v_responses.mat', 'v_responses');
load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
load([retroot '/width_maps/gt/all_gt_widths.mat'], 'all_gt_widths');

v_responses = reshape(v_responses, [], 1, 5, 3);
true_oris = angle(gt_orientations) / 2;
%%
for i_scale = 6%1:6
    if i_scale < 6
        [~, predicted_theta] = gaussian_2nd_derivative_line(v_responses(:,:,i_scale,:));
    else
        [~, predicted_theta] = gaussian_2nd_derivative_line(v_responses(:,:,1:4,:));
    end
    predicted_ori = complex(cos(2*predicted_theta), sin(2*predicted_theta));
    [~, err_stats] = ori_error(predicted_ori, true_oris);
    display(['Scale: ' num2str(i_scale)]);
    display(err_stats);
end
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
comparing_DRIVE_experiment_csf(1, 1,  ... % non-strict mode
    'image_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/images/'], ...
    'vessel_mask_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/vessel_masks/'], ...
    'fov_mask_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/fov_masks/'], ...
    'ori_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/orientations/'], ...
    'width_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/width_maps/'], ...
    'selected_images', [],...
    'num_pts',      unixenv('NUM_SAMPLES', 100000), ...
    'num_trees',    unixenv('NUM_TREES', 2), ...
    'make_data',    unixenv('MAKE_DATA',1), ...
    'do_orientation', unixenv('DO_ORIENTATION',0), ...
    'do_detection', unixenv('DO_DETECTION',0), ...
    'do_width', unixenv('DO_WIDTH',0), ...
    'do_tests', unixenv('DO_TESTS',1:5), ...
    'win_sizes', unixenv('WIN_SIZES',1) ...
);
%%
DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 qsub -V -t 1-100 -l short matlab_code/trunk/hydra/cuc/generic_batch_fun.sh
%%
warning('off', 'ASYM:unexpectedArgument');
pts_per_im = 10;
dx = 64;
dy = 64;
cx = 32;
cy = 32;

d_args{1}.decomp_type = {'dt'};
d_args{1}.levels = 1:4;
d_args{1}.feature_shape = 'rect';
d_args{1}.feature_type = 'complex';
d_args{1}.do_max = 0;
d_args{1}.rotate = 0;
d_args{1}.use_nag = 0;

d_args{2}.decomp_type = {'linop'};
d_args{2}.num_levels = 4;
d_args{2}.num_angles = 6;
d_args{2}.do_max = 0;
d_args{2}.rotate = 0;       

d_args{3}.decomp_type = {'gabor'};
d_args{3}.num_angles = 6;
d_args{3}.sigma_range = [1 2 4 8];	
d_args{3}.do_max = 0;
d_args{3}.rotate = 0;
d_args{3}.feature_type = 'complex';

d_args{4}.decomp_type = {'mono'};
d_args{4}.num_levels = 4;
d_args{4}.min_wavelength = 4;
d_args{4}.onf = 0.65;

d_args{5}.decomp_type = {'g1d'};
d_args{5}.sigma_range = [1 2 4 8];
            
d_args{6}.decomp_type = {'g2d'};
d_args{6}.sigma_range = [1 2 4 8];
            
d_args{7}.decomp_type = {'g2di'};
d_args{7}.sigma_range = [1 4];
            
d_args{8}.decomp_type = {'h2d'};
d_args{8}.sigma_range = [1 2 4 8];

d_args{9}.decomp_type = {'g2da'};
d_args{9}.sigma_range = [1 2 4 8];
d_args{9}.num_angles = 6;
d_args{9}.do_max = 0;
d_args{9}.rotate = 0;

d_args{10}.decomp_type = {'gabori'};
d_args{10}.sigma_range = [1 4];
d_args{10}.num_angles = 6;
d_args{10}.do_max = 0;
d_args{10}.feature_type = 'complex';

num_decomps = 10;

%set common parameters for all decomp types and then compute vector sizes
    

%--------------------------------------------------------------------------
% ***************** GENERATE DATA *****************************************
%--------------------------------------------------------------------------
sampling_time = zeros(100,2,10);

for i_decomp = [1 3 4 5 6 7 10]
    d_args{i_decomp}.win_size = 3;
    d_args{i_decomp}.normalise = 0;
    d_args{i_decomp}.pca = [];

    D = zeros(num_decomps,1);
    D(i_decomp) = get_samples_per_channel(d_args{i_decomp});

    for ii = 1:100

        display(['Testing image ' num2str(ii)]);

        %Sample properties of line
        line_width = sample_uniform([1 8]);
        line_contrast = sample_uniform([1 8]);
        line_ori = sample_uniform([0 360]);
        line_rad = pi * line_ori / 180;

        %Generate line
        [line, label, label_centre] =...
            create_ellipse_bar(line_width/2, line_contrast, line_ori, dy, dx, cx, cy);
        
        tic;
        %For ecah decomposition compute the responses
        im_responses = compute_filter_responses(line, d_args{i_decomp});
        sampling_time(ii,1,i_decomp) = toc;
        
        fov_mask = true(size(label));
        fov_mask([1:16 end-15:end],:) = 0;
        fov_mask(:, [1:16 end-15:end]) = 0;

        %Select some random pixels in the signal image
        shuffle = randperm(sum(label(:) & fov_mask(:)));
        idx = find(label & fov_mask);
        r_idx = idx(shuffle(1:pts_per_im));
        [r_rows r_cols] = ind2sub([dy dx], r_idx);
        sample_idx = ((ii-1)*pts_per_im + 1):(ii*pts_per_im);

        tic;
        sample_image_features(im_responses, r_rows(:), r_cols(:), d_args{i_decomp});
        toc;
        sampling_time(ii,2,i_decomp) = toc;
    end
end
%%
for i_decomp = [1 3 4 5 6 7 10]
    display([d_args{i_decomp}.decomp_type ':, ' num2str(1000*mean(sampling_time(:,:,i_decomp)),3) num2str(1000*mean(sum(sampling_time(:,:,i_decomp),2)),3)]);
end






