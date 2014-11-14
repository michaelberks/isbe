retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];

load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
num_angles = [4 5 6 8 12 15 18 24 30];
n_tests = length(num_angles);
mean_abs_err = zeros(n_tests,1);
median_abs_err = zeros(n_tests,1);
results_dir = ['C:\isbe\asymmetry_project\experiments\DRIVE\gabor_orientation_prediction\'...
    'sigma_all\'];
        
for i_angle = 1:n_tests
    
    rad_angles = pi/2 + pi*(0:(num_angles(i_angle)-1))'/num_angles(i_angle);
    complex_angles = complex(cos(2*rad_angles), sin(2*rad_angles));
    predicted_orientations = [];
    
    for i_ret = 21:40
        display(['Analysing image ' num2str(i_ret)]);
        
        %load vessel mask
        ret = rgb2gray(u_load([retroot 'images\' zerostr(i_ret,2) '_training.mat']));
        v_mask = u_load([retroot 'vessel_masks\' zerostr(i_ret,2) '_training_v_mask.mat']);
        f_mask = u_load([retroot 'fov_masks\' zerostr(i_ret,2) '_training_f_mask.mat']);
        gt_ori_map = u_load([retroot 'orientations\' zerostr(i_ret,2) '_ori1.mat']);
        v_mask = v_mask & f_mask;
        
        for i_sigma = [1 2 4 8]
            
            %Load responses
%             responses = u_load([retroot 'filter_responses\gabor\'...
%                         zerostr(i_ret,2) '_responses_' zerostr(i_sigma, 2) '.mat']);
%             
%             %Select only the relevant angles
%             responses = responses(:,:,:,angles);
            [responses] = compute_gabor_responses(ret, i_sigma, num_angles(i_angle));
            
            %Compute orientation as arg max of absolute response
            [max_responses_i, max_angles_i] = max(abs(responses), [], 4);
            clear responses;                
            
            if i_sigma == 1
                %For sigma = 1, create the max_responses and angles
                max_responses = max_responses_i;%(v_mask)
                max_angles = max_angles_i;%(v_mask)
            else
                %For further scales, swap in the responses and angles if
                %they're greater than the existing responses
                %max_responses_i = max_responses_i(v_mask);
                %max_angles_i = max_angles_i(v_mask);
                
                swap_idx = max_responses_i > max_responses;
                max_responses(swap_idx) = max_responses_i(swap_idx);
                max_angles(swap_idx) = max_angles_i(swap_idx);
            end
        end
        predicted_ori_map = complex_angles(max_angles);
        
        if i_ret == 1           
            figure; 
            subplot(1,2,1); imgray(complex2rgb(predicted_ori_map .* abs(max_responses)));
            subplot(1,2,2); imgray(complex2rgb(gt_ori_map));
        end
        
        %Convert the max angles to an angle in complex form and add to the
        %predictions list
        predicted_orientations = [predicted_orientations; predicted_ori_map(v_mask)]; %#ok
        
    end
    
    %Now get the error stats
    [orientation_errors, error_stats] = ori_error(predicted_orientations, gt_orientations);
    mean_abs_err(i_angle) = error_stats.abs_mean;
    median_abs_err(i_angle) = error_stats.abs_median;    
    save([results_dir 'ori_errors_' zerostr(num_angles(i_angle),2) '.mat'], 'orientation_errors');
end

figure; hold all;
plot(num_angles, mean_abs_err, '-.x');
plot(num_angles, median_abs_err, '-.x');
legend({'Mean absolute errors', 'Median absolute errors'});
title('Orientation errors for increasing angluar resolution of Gabor filter - DRIVE data');
xlabel('Number of filters');
%%
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];

load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
num_angles = [4 5 6 8 12 15 18 24 30 36 45 60];
n_tests = length(num_angles);
mean_abs_err = zeros(n_tests,1);
median_abs_err = zeros(n_tests,1);

for i_sigma = [1 2 4 8]
    
    results_dir = ['C:\isbe\asymmetry_project\experiments\DRIVE\gabor_orientation_prediction\'...
        'sigma' num2str(i_sigma) '\'];
    mkdir(results_dir);
    
    for i_angle = 1:n_tests

        rad_angles = pi/2 + pi*(0:(num_angles(i_angle)-1))'/num_angles(i_angle);
        complex_angles = complex(cos(2*rad_angles), sin(2*rad_angles));
        predicted_orientations = [];

        for i_ret = 21:40
            display(['Analysing image ' num2str(i_ret)]);

            %load vessel mask
            ret = rgb2gray(u_load([retroot 'images\' zerostr(i_ret,2) '_training.mat']));
            v_mask = u_load([retroot 'vessel_masks\' zerostr(i_ret,2) '_training_v_mask.mat']);
            f_mask = u_load([retroot 'fov_masks\' zerostr(i_ret,2) '_training_f_mask.mat']);
            gt_ori_map = u_load([retroot 'orientations\' zerostr(i_ret,2) '_ori1.mat']);
            v_mask = v_mask & f_mask;

            %Compute responses
            [responses] = compute_gabor_responses(ret, i_sigma, num_angles(i_angle));
            
            %Compute orientation as arg max of absolute response
            [~, max_angles] = max(abs(responses), [], 4);
            clear responses;               
            predicted_ori_map = complex_angles(max_angles);
        
            %Convert the max angles to an angle in complex form and add to the
            %predictions list
            predicted_orientations = [predicted_orientations; predicted_ori_map(v_mask)]; %#ok
        end
        
        %Now get the error stats
        [orientation_errors, error_stats] = ori_error(predicted_orientations, gt_orientations);
        mean_abs_err(i_angle) = error_stats.abs_mean;
        median_abs_err(i_angle) = error_stats.abs_median;  
        save([results_dir 'ori_errors_' zerostr(num_angles(i_angle),2) '.mat'], 'orientation_errors');
        
    end
    
    figure; hold all;
    plot(num_angles, mean_abs_err, '-.x');
    plot(num_angles, median_abs_err, '-.x');
    legend({'Mean absolute errors', 'Median absolute errors'});
    title('Orientation errors for increasing angluar resolution of Gabor filter - DRIVE data');
    xlabel('Number of filters');
end
%%
num_angles = [4 5 6 8 12 15 18 24 30 36 45 60];
n_tests = length(num_angles);
leg_text = cell(0,1);
figure; a1 = gca; hold all; title('Mean absolute errors for varying angular resolution of Gabor filters');
figure; a2 = gca; hold all; title('Median absolute errors for varying angular resolution of Gabor filters');
figure; a3 = gca; hold all; title('CDF of absolute errors for varying angular resolution of Gabor filters'); 
for i_sigma = [1 2 4 8 9]
    
    if i_sigma == 9
        results_dir = ['C:\isbe\asymmetry_project\experiments\DRIVE\gabor_orientation_prediction\'...
            'sigma_all\'];
        leg_text{end+1} = '\sigma = [1 2 4 8]';
    else
        results_dir = ['C:\isbe\asymmetry_project\experiments\DRIVE\gabor_orientation_prediction\'...
            'sigma' num2str(i_sigma) '\'];
        leg_text{end+1} = ['\sigma = ' num2str(i_sigma)];
    end
    
    mean_abs_err = zeros(n_tests,1);
    median_abs_err = zeros(n_tests,1);

    for i_angle = 1:n_tests        
        load([results_dir 'ori_errors_' zerostr(num_angles(i_angle),2) '.mat'], 'orientation_errors');
        orientation_errors(isnan(orientation_errors)) = [];
        mean_abs_err(i_angle) = mean(abs(orientation_errors))*180/pi;
        median_abs_err(i_angle) = median(abs(orientation_errors))*180/pi;
        
        if num_angles(i_angle) == 30 
            sorted_errors = sort(abs(orientation_errors))*180/pi;
            cdf_pts = round(linspace(1, length(orientation_errors), 1001));
            %plot(a3, sorted_errors(cdf_pts), 0:100, '-x');
            plot(a3, sorted_errors(cdf_pts), 0:1000, '-x');
        end
        
    end
    plot(a1, num_angles, mean_abs_err, '-.x');
    plot(a2, num_angles, median_abs_err, '-.x');
    
end
legend(a1, leg_text);
legend(a2, leg_text);
legend(a3, leg_text);

    
%%
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];

decomp_types = cell(4, 2);
decomp_types( 1,:) = {'g2d', 1};
decomp_types( 2,:) = {'g1d', 0};
decomp_types( 3,:) = {'gabor', 1};
decomp_types( 4,:) = {'mono', 0};
decomp_types( 5,:) = {'h2d', 0};
decomp_types( 6,:) = {'g', 0};


for i_decomp = 1:length(decomp_types)
    decomp_type = decomp_types{i_decomp,1};
    if decomp_types{i_decomp,2}
        mkdir([retroot 'filter_responses\' decomp_type]);
    else
        continue;
    end

    for i_ret = 21:40
        %load retinogram and merge RGB channels
        ret = u_load([retroot '\images\' zerostr(i_ret,2) '_training.mat']);
        ret = rgb2gray(ret);

        switch decomp_type

            case {'gabor'}
                sigma_range = [1 2 4 8 16];
                for i_sigma = sigma_range
                    
                    load([retroot 'filter_responses\' decomp_type '/'...
                        zerostr(i_ret,2) '_responses_' zerostr(i_sigma, 2) '.mat'], 'responses');
                    save([retroot 'filter_responses\' decomp_type '/'...
                        zerostr(i_ret,2) '_responses_' zerostr(i_sigma, 2) '.mat'], 'responses');
                end

            case {'g2d'}
                sigma_range = [1 2 4 8 16];
                for i_sigma = sigma_range
                    load([retroot 'filter_responses\' decomp_type '/'...
                        zerostr(i_ret,2) '_responses_' zerostr(i_sigma, 2) '.mat'], 'responses');
                    save([retroot 'filter_responses\' decomp_type '/'...
                        zerostr(i_ret,2) '_responses_' zerostr(i_sigma, 2) '.mat'], 'responses');
                end

        end
    end
end
%%
[rows cols] = size(ret);



g2di_responses = compute_gaussian_2nd_derivatives_d(ret, 1, 4);
g2d_responses = compute_gaussian_2nd_derivatives(ret, [1 2 4 8]);

for i_scale = 1:4
    figure;
    for i_band = 1:3
        subplot(2,3,i_band); imgray(g2d_responses(:,:,i_scale,i_band));
        subplot(2,3,3+i_band); imgray(g2di_responses{i_scale}(:,:,i_band));
    end
end
%%
xx = repmat(1:cols,rows,1);
yy = repmat((1:rows)',1,cols);
%%
xx = 1:cols;
yy = (1:rows)';

for i_scale = 1:4
    figure;
    for i_band = 1:3
        scaling = 2^(i_scale-1);
        offset = (i_scale-1) / (2^i_scale);
        band = g2d_responses(:,:,i_scale,i_band);
        band_i = g2di_responses{i_scale}(:,:,i_band);
        band_i = interp2(band_i, offset+xx/scaling, offset+yy/scaling, 'cubic');
        %band_i = imresize(band_i, [rows cols], 'cubic');
        
        subplot(2,3,i_band); imgray(band);
        subplot(2,3,3+i_band); imgray(abs(band-band_i)); colorbar;
    end
end
%%
[line_response, orientation, scale, normal_response] = gaussian_2nd_derivative_line(ret, [1 2 4 8]);
[line_response_r, orientation_r, scale_r, normal_response_r] = gaussian_2nd_derivative_line(g2d_responses);
[line_response_i, orientation_i, scale_i, normal_response_i] = gaussian_2nd_derivative_line(g2di_responses);
%%
figure; 
subplot(2,3,1); imgray(line_response);
subplot(2,3,2); imgray(line_response_r);
subplot(2,3,3); imgray(line_response_i);
%subplot(2,3,4); imgray(line_response);
subplot(2,3,5); imgray(abs(line_response-line_response_r));
subplot(2,3,6); imgray(abs(line_response-line_response_i));
%%
discard_idx = false(0,1);
for i_ret = 21:40
    display(['Analysing image ' num2str(i_ret)]);

    %load vessel mask
    v_mask = u_load([retroot 'vessel_masks\' zerostr(i_ret,2) '_training_v_mask.mat']);
    f_mask = u_load([retroot 'fov_masks\' zerostr(i_ret,2) '_training_f_mask.mat']);
    discard_idx = [discard_idx; ~f_mask(v_mask)];
end
%%
predict_image_set(...
    'model_id', '14976',...
    'image_dir',        'images',...
    'prediction_dir',   'predictions',...
    'model_name',       'predictor',...
    'use_sampled_maps', 1,...
    'model_root',		[asymmetryroot,'data/models/vessel/detection/rf_classification/'], ...
    'image_root',       [asymmetryroot,'data/retinograms/DRIVE/test'],...
    'task_id',			1, ...
    'num_jobs',			20, ...
    'mask_dir',			'fov_masks');
%%
deg_angles = 90 + 180*(0:11)/12;
complex_angles = complex(cosd(2*rad_angles), sind(2*rad_angles));
for theta = 0:15:165
    line = create_ellipse_bar(2, 1, theta, 32, 32, 16, 16);
    [gabor_responses] = compute_gabor_responses(line, 2, 12);
    [max_responses_i, max_angles_i] = max((real(gabor_responses)), [], 4);

    ori_map = complex_angles(max_angles_i);
    figure; 
    subplot(1,3,1); imgray(complex2rgb(ori_map .* abs(max_responses_i)));
    %subplot(1,3,2); imgray(complex2rgb(ori_map .* real(max_responses_i)));
    %subplot(1,3,3); imgray(complex2rgb(ori_map .* imag(max_responses_i)));
    
    max_band = max_angles_i(16,16);
    display(['Theta = ' num2str(theta) ', max band = ' num2str(max_band) ', predicted angle = ' num2str(deg_angles(max_band))]);
end
%%
figure;
a1 = subplot(1,2,1); hold all;
a2 = subplot(1,2,2); hold all;
for line_width = 1:8
    line = create_ellipse_bar(line_width/2, 1, theta, 32, 32, 16, 16);
    
    scale_responses_gabor = zeros(8,1);
    scale_responses_gauss = zeros(8,1);
    for sigma = 1:8
        [gabor_responses] = compute_gabor_responses(line, sigma, 2);
        [gauss_responses] = gaussian_2nd_derivative_line(line, sigma);
        scale_responses_gabor(sigma) = max(abs(gabor_responses(16,16,:,:)), [], 4);
        scale_responses_gauss(sigma) = abs(gauss_responses(16,16));
    end
    plot(a1, 1:8, scale_responses_gabor, '-x');
    plot(a2, 1:8, scale_responses_gauss, '-x');
end
legend(a1, num2str((1:8)'));
legend(a2, num2str((1:8)'));
        
%%
ret = u_load('C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\images\21_training.mat');
ret = ret(:,:,2);
rad_angles = pi/2 + pi*(0:59)'/60;
complex_angles = complex(cos(2*rad_angles), sin(2*rad_angles));
%    
for num_angles = [6 12 20]
    i_step = 60 / num_angles;
    angles = 1:i_step:60;
    
    for sigma = [1 2 4]
    
        [gabor_responses] = compute_gabor_responses(ret, sigma, num_angles);
        [max_responses_i, max_angles_i] = max(abs(gabor_responses), [], 4);
        
        if sigma == 1
            %For sigma = 1, create the max_responses and angles
            max_responses = max_responses_i;
            max_angles = max_angles_i;
        else
            %For further scales, swap in the responses and angles if
            %they're greater than the existing responses
            swap_idx = max_responses_i > max_responses;
            max_responses(swap_idx) = max_responses_i(swap_idx);
            max_angles(swap_idx) = max_angles_i(swap_idx);
        end
    end
    
    ori_map = complex_angles(angles(max_angles));
    figure; imgray(complex2rgb(ori_map .* abs(max_responses)));
end
%%
figure;
a1 = subplot(1,2,1); hold all;
a2 = subplot(1,2,2); hold all;
for line_width = 1:8
    line = create_sin_bar(line_width/2, 1, 0, 32, 32, 0, 16, 16);
    
    scale_responses_gabor = zeros(8,1);
    scale_responses_gauss = zeros(8,1);
    for sigma = 1:8
        [gabor_responses] = compute_gabor_responses(line, sigma, 2);
        [gauss_responses] = gaussian_2nd_derivative_line(line, sigma);
        scale_responses_gabor(sigma) = max(abs(gabor_responses(16,16,:,:)), [], 4);
        scale_responses_gauss(sigma) = abs(gauss_responses(16,16));
    end
    plot(a1, 1:8, scale_responses_gabor, '-x');
    plot(a2, 1:8, scale_responses_gauss, '-x');
end
legend(a1, num2str((1:8)'));
legend(a2, num2str((1:8)'));
%%
% Gabor filter responses for synthetic lines at different scales - try
% swapping over the normalising constant in gabor_filters
for f_type = 1:2
    if f_type == 1
        feat = 'line';
    else
        feat = 'edge';
    end
    
    f1 = figure;
    a1 = subplot(2,2,1); hold all; 
    title(['Abs Gabor response to ' feat]); xlabel('\sigma of filter');
    a2 = subplot(2,2,2); hold all; 
    title(['G" response to ' feat]); xlabel('\sigma of filter');
    a3 = subplot(2,2,3); hold all; 
    title(['Real Gabor response to ' feat]); xlabel('\sigma of filter');
    a4 = subplot(2,2,4); hold all; 
    title(['Imag Gabor response to ' feat]); xlabel('\sigma of filter');
    
    for line_width = 1:8
        if f_type == 1
            f_image = create_sin_bar(line_width, 1, 0, 32, 32,0, 16, 16) + 3;
        else
            f_image = create_sin_step(line_width, 1, 0, 32, 32, 16, 16) + 3;
        end

        scale_responses_gabor = zeros(8,1);
        scale_responses_gauss = zeros(8,1);
        for sigma = 1:8
            [gabor_responses] = compute_gabor_responses(f_image, sigma, 2);
            [gauss_responses] = gaussian_2nd_derivative_line(f_image, sigma);
            scale_responses_gabor(sigma) = max(gabor_responses(16,16,:,:), [], 4);
            scale_responses_gauss(sigma) = abs(gauss_responses(16,16));
            
            if line_width == sigma
                figure;
                %Take maximum response across orientations
                [max_responses] = max(gabor_responses(:,:,:,:), [], 4);

                %Look at abs, real and imag parts of responses
                subplot(1,3,1); imgray(abs(max_responses)); colorbar('southoutside'); 
                title(['Absolute response, \sigma = ' num2str(sigma)]);
                subplot(1,3,2); imgray(real(max_responses)); colorbar('southoutside');
                title(['Real response, \sigma = ' num2str(sigma)]);
                subplot(1,3,3); imgray(imag(max_responses)); colorbar('southoutside');
                title(['Imaginary response, \sigma = ' num2str(sigma)]);
            end
        end
        figure(f1);
        plot(a1, 1:8, abs(scale_responses_gabor), '-x');
        plot(a2, 1:8, scale_responses_gauss, '-x');
        plot(a3, 1:8, real(scale_responses_gabor), '-x');
        plot(a4, 1:8, imag(scale_responses_gabor), '-x');
            
    end
    legend(a1, strcat('Feature width =  ', num2str((1:8)')));
    legend(a2, strcat('Feature width =  ', num2str((1:8)')));
    legend(a2, strcat('Feature width =  ', num2str((1:8)')));
    legend(a3, strcat('Feature width =  ', num2str((1:8)')));
end

%%
%Gabor responses for the retinograms
%Load retinogram
ret = u_load([asymmetryroot 'data\retinograms\DRIVE\training\images\21_training.mat']);
ret = ret(:,:,2);
for sigma = [1 2 4 8]
    %Compute Gabor responses
    [gabor_responses] = compute_gabor_responses(ret, sigma, 6);
    
    %Take maximum response across orientations
    [max_responses] = max(gabor_responses(:,:,:,:), [], 4);

    %Look at abs, real and imag parts of responses
    figure;
    a(1) = subplot(1,3,1); imgray(abs(max_responses)); colorbar('southoutside'); 
    title(['Absolute response, \sigma = ' num2str(sigma)]);
    a(2) = subplot(1,3,2); imgray(real(max_responses)); colorbar('southoutside');
    title(['Real response, \sigma = ' num2str(sigma)]);
    a(3) = subplot(1,3,3); imgray(imag(max_responses)); colorbar('southoutside');
    title(['Imaginary response, \sigma = ' num2str(sigma)]);
    linkaxes(a);
    
    %You might need to zoom to see the responses at vessel centres - the
    %axes are linked
end
%%
ret = u_load([asymmetryroot 'data\retinograms\DRIVE\training\images\21_training.mat']);
ret = ret(:,:,2);
[gabor_responses] = compute_gabor_responses(ret, [1 2 4], 3);
[gabor_responses_d] = compute_gabor_responses_d(ret, 1, 3, 3);

for level = 1:3
    figure;
    for band = 1:3
        a(band) = subplot(2,3,band); imgray(abs(gabor_responses(:,:,level,band)));
        a(band+3) = subplot(2,3,band+3); imgray(abs(gabor_responses_d{level}(:,:,band)));    
    end
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
%Load and reformat data
clear; pack;
retroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];
load('C:\isbe\asymmetry_project\experiments\DRIVE\gabor_orientation_prediction\v_responses.mat', 'v_responses');
load([retroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
load([retroot '/width_maps/gt/all_gt_widths.mat'], 'all_gt_widths');

v_responses = reshape(v_responses, [], 1, 5, 6);
true_oris = angle(gt_orientations) / 2;
num_bins = 50;
width_bins = linspace(min(all_gt_widths), max(all_gt_widths), num_bins);
ori_bins = linspace(-pi/2, pi/2, num_bins);

%figure; a1 = gca; hold all;
for i_scale = 1:6
    if i_scale < 6
        [~, predicted_bands] = ...
            max_response_line(v_responses(:,:,i_scale,:));
        scale_str = num2str(i_scale);
    else
        [~, predicted_bands, predicted_scales] = ...
            max_response_line(v_responses(:,:,:,:));
        scale_str = 'all';
    end

    %predicted_oris = complex_angles(predicted_bands);
    %[~, err_stats] = ori_error(complex_oris, predicted_oris);

    %plot(a1, err_stats.abs_percentiles, 1:100);

    counts = zeros(6, num_bins);
    for i_band = 1:6
        counts(i_band,:) = hist(true_oris(predicted_bands == i_band), ori_bins);
    end
    figure; plot(ori_bins, bsxfun(@rdivide, counts, sum(counts))');
    xlabel('Orientation of vessel');
    ylabel('Percentage of pixels with maximum response');
    title({'Which oriented sub-band produces maximum response?'; ['Scale: ' scale_str]});
end
%%
counts = zeros(5, 16);
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
my_root = [asymmetryroot 'experiments/synthetic_lines/comparing_representations/'];
isbe_root = 'Y:\asymmetry_project\experiments\synthetic_lines\comparing_representations\';
for i_decomp = 1:10
    for repeat = 1:11
        for noise_level = 1:4
            for data_type = {'training', 'test'}
                my_dir = [my_root data_type{1} '/rician_' num2str(noise_level) '/' num2str(repeat)];
                isbe_dir = [isbe_root data_type{1} '/rician_' num2str(noise_level) '/' num2str(repeat)];
                %mkdir(my_dir);
                if ~exist([my_dir '/true_labels.mat'], 'file') && exist([isbe_dir '/true_labels.mat'], 'file')
                    copyfile([isbe_dir '/true_labels.mat'], [my_dir '/true_labels.mat']);
                end
                if ~exist([my_dir '/responses_g2d.mat'], 'file') && exist([isbe_dir '/responses_g2d.mat'], 'file')
                    copyfile([isbe_dir '/responses_g2d.mat'], [my_dir '/responses_g2d.mat']);
                end
                if ~exist([my_dir '/responses_g2da.mat'], 'file') && exist([isbe_dir '/responses_g2da.mat'], 'file')
                    copyfile([isbe_dir '/responses_g2da.mat'], [my_dir '/responses_g2da.mat']);
                end
                if ~exist([my_dir '/responses_gabor.mat'], 'file') && exist([isbe_dir '/responses_gabor.mat'], 'file')
                    copyfile([isbe_dir '/responses_gabor.mat'], [my_dir '/responses_gabor.mat']);
                end
            end
        end
    end
end
    
    