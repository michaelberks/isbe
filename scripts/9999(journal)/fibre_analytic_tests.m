%**************************************************************************
%***** Script to produce results for CVPR submission Nov 2011 using *******
%****************** the DRIVE retinography dataset ************************
%**************************************************************************
%**************************************************************************
warning('off', 'load_uint8:missing_variables');

clear;

%% 1. Produce response, orientation and scale maps for the analytic methods

fibreroot = [asymmetryroot 'data\fibre\test\'];
im_list = dir([fibreroot 'images\*.mat']);

pred_dir = [fibreroot 'predictions\orientation\analytic\'];
do_mono = 0;
do_g1d = 0;
do_g2d = 0;
do_gabor = 1;

if do_mono
    mkdir([pred_dir 'mono\orientations']);
    mkdir([pred_dir 'mono\responses']);
    mkdir([pred_dir 'mono\scales']);
end
if do_g1d
    mkdir([pred_dir 'g1d\orientations']);
    mkdir([pred_dir 'g1d\responses']);
    mkdir([pred_dir 'g1d\scales']);
end
if do_g2d
    mkdir([pred_dir 'g2d\orientations']);
    mkdir([pred_dir 'g2d\responses']);
    mkdir([pred_dir 'g2d\scales']);
end
if do_gabor
    mkdir([pred_dir 'gabor\orientations']);
    mkdir([pred_dir 'gabor\responses']);
    mkdir([pred_dir 'gabor\scales']);
end

for ii = 1:length(im_list)
    %load retinogram and merge RGB channels
    fibre_im = u_load([fibreroot 'images\' im_list(ii).name]);
    if size(fibre_im,3) == 3
        fibre_im = rgb2gray(fibre_im);
    end
    
    if do_mono
        display(['Computing monogenic response for image ' num2str(ii)]);
        %Compute repsonses for mono and save
        [response_map d ori_map scale_map] = monogenic_multiscale(fibre_im, 2, 2, 2, 0.65);
        ori_map = abs(response_map).*exp(2i*ori_map);
        save_uint8([pred_dir 'mono\orientations\' zerostr(ii,3) '_' im_name '_ori.mat'], ori_map);
        save_uint8([pred_dir 'mono\responses\' zerostr(ii,3) '_' im_name '_response.mat'], response_map);
        save_uint8([pred_dir 'mono\scales\' zerostr(ii,3) '_' im_name '_scale.mat'], scale_map);
    end
    if do_g1d
        display(['Computing G1d response for image ' num2str(ii)]);
        %Compute repsonses for mono and save
        [response_map ori_map scale_map] = gaussian_1st_derivative_gradient(fibre_im, [2 4]);
        ori_map = abs(response_map).*exp(2i*ori_map);
        save_uint8([pred_dir 'g1d\orientations\' zerostr(ii,3) '_' im_name '_ori.mat'], ori_map);
        save_uint8([pred_dir 'g1d\responses\' zerostr(ii,3) '_' im_name '_response.mat'], response_map);
        save_uint8([pred_dir 'g1d\scales\' zerostr(ii,3) '_' im_name '_scale.mat'], scale_map);
    end
    if do_g2d
        display(['Computing G2d response for image ' num2str(ii)]);
        %Compute repsonses for mono and save
        [response_map ori_map scale_map] = gaussian_2nd_derivative_line(fibre_im, [2 4]);
        ori_map = abs(response_map).*exp(2i*ori_map);
        save_uint8([pred_dir 'g2d\orientations\' zerostr(ii,3) '_' im_name '_ori.mat'], ori_map);
        save_uint8([pred_dir 'g2d\responses\' zerostr(ii,3) '_' im_name '_response.mat'], response_map);
        save_uint8([pred_dir 'g2d\scales\' zerostr(ii,3) '_' im_name '_scale.mat'], scale_map);
    end
    if do_gabor
        display(['Computing Gabor response for image ' num2str(ii)]);
        %Compute repsonses for mono and save
        num_angles = 12;
        [responses] = compute_gabor_responses(fibre_im, 4, num_angles);
        [response_map, max_band, scale_map] = max_response_line(responses);
        clear responses;
        
        rad_angles = pi/2 + pi*(0:(num_angles-1))'/num_angles;
        complex_angles = complex(cos(2*rad_angles), sin(2*rad_angles));
        ori_map = abs(response_map).*complex_angles(max_band);
        
        save_uint8([pred_dir 'gabor\orientations\' zerostr(ii,3) '_' im_name '_ori.mat'], ori_map);
        save_uint8([pred_dir 'gabor\responses\' zerostr(ii,3) '_' im_name '_response.mat'], response_map);
        save_uint8([pred_dir 'gabor\scales\' zerostr(ii,3) '_' im_name '_scale.mat'], scale_map);
    end
end
%%
fibreroot = [asymmetryroot 'data\fibre\test\'];
fov_mask_dir = [fibreroot 'fov_masks\'];
fg_mask_dir = [fibreroot 'fibre_masks\'];
label_dir = [fibreroot 'orientation_maps\'];

warning('off', 'ori_error:nans');
pred_dir = [fibreroot 'predictions\orientation\analytic\'];
label_dir = [fibreroot 'orientations\'];

vessel_centres = logical(vessel_centres);

analytic_codes = cell(4, 1);
analytic_codes( 1,:) = {'g1d'};
analytic_codes( 2,:) = {'g2d'};
analytic_codes( 3,:) = {'mono'}; %Old version of filters '11480'
analytic_codes( 4,:) = {'gabor'};

error_medians = zeros(4,2);
error_means = zeros(4,2);

for ii = 1:4

    if 0%exist([pred_dir analytic_codes{ii,1} '\errors\error_stats.mat'], 'file')
        load([pred_dir analytic_codes{ii,1} '\errors\error_stats.mat']);
    else
        [orientation_errors, ~, error_stats] =...
            compute_image_orientation_errors([pred_dir analytic_codes{ii,1} '\orientations\'], fg_mask_dir,...
            'label_dir', label_dir, 'fov_mask_dir', fov_mask_dir,...
            'save_dir', [pred_dir analytic_codes{ii,1} '\errors\']);
    end
    centre_stats = ori_error_stats(orientation_errors(vessel_centres));
    
    error_medians(ii,1) = error_stats.abs_median;
    error_medians(ii,2) = centre_stats.abs_median;
    
    error_means(ii,1) = error_stats.abs_mean;
    error_means(ii,2) = centre_stats.abs_mean;
    
    display(['Errors for ' analytic_codes{ii,1} ...
        ', abs mean: ' num2str(error_stats.abs_mean,3) ', abs_median: ' num2str(error_stats.abs_median,3)]); 
       
end
%%
% 
% fibreroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];
% fov_mask_dir = [fibreroot 'fov_masks\'];
% fg_mask_dir = [fibreroot 'vessel_masks\'];
% 
% decomp_types = cell(7, 2);
% decomp_types( 1,:) = {'gabor', 0};
% decomp_types( 2,:) = {'g', 0};
% decomp_types( 3,:) = {'g1d', 0};
% decomp_types( 4,:) = {'g2d', 0};
% decomp_types( 5,:) = {'g2di', 1};
% decomp_types( 6,:) = {'h2d', 0};
% decomp_types( 7,:) = {'mono', 0};
% 
% 
% for i_decomp = 1:length(decomp_types)
%     decomp_type = decomp_types{i_decomp,1};
%     if decomp_types{i_decomp,2}
%         mkdir([fibreroot 'filter_responses\' decomp_type]);
%     else
%         continue;
%     end
% 
%     for i_ret = 21:40
%         %load retinogram and merge RGB channels
%         ret = u_load([fibreroot '\images\' zerostr(i_ret,2) '_training.mat']);
%         ret = rgb2gray(ret);
% 
%         switch decomp_type
% 
%             case {'gabor'}
%                 sigma_range = [1 2 4 8 16];
%                 num_angles = 60;
%                 for i_sigma = sigma_range
%                     
%                     responses = compute_gabor_responses(ret, i_sigma, num_angles);
%                     save([fibreroot 'filter_responses\' decomp_type '/'...
%                         zerostr(i_ret,2) '_responses_' zerostr(i_sigma, 2) '.mat'], 'responses');
%                 end
%             case {'g'}
%                 sigma_range = [1 2 4 8 16];
%                 for i_sigma = sigma_range                  
%                     responses = compute_gaussian_responses(ret, i_sigma);
%                     save([fibreroot 'filter_responses\' decomp_type '/'...
%                         zerostr(i_ret,2) '_responses_' zerostr(i_sigma, 2) '.mat'], 'responses');
%                 end
% 
%             case {'g1d'}
%                 sigma_range = [1 2 4 8 16];
%                 for i_sigma = sigma_range 
%                     responses = compute_gaussian_1st_derivatives(ret, i_sigma);
%                     save([fibreroot 'filter_responses\' decomp_type '/'...
%                         zerostr(i_ret,2) '_responses_' zerostr(i_sigma, 2) '.mat'], 'responses');
%                 end
% 
%             case {'g2d'}
%                 sigma_range = [1 2 4 8 16];
%                 for i_sigma = sigma_range
%                     responses = compute_gaussian_2nd_derivatives(ret, i_sigma);
%                     save([fibreroot 'filter_responses\' decomp_type '/'...
%                         zerostr(i_ret,2) '_responses_' zerostr(i_sigma, 2) '.mat'], 'responses');
%                 end
%                 
%             case {'g2di'}
%                 responses = compute_gaussian_2nd_derivatives_d(ret, 1, 5);
%                 save([fibreroot 'filter_responses\' decomp_type '/'...
%                     zerostr(i_ret,2) '_responses.mat'], 'responses');
% 
%             case {'h2d'}
%                 sigma_range = [1 2 4 8 16];
%                 for i_sigma = sigma_range
%                     responses = compute_hilbert_2nd_derivatives(ret, sigma_range);
%                     save([fibreroot 'filter_responses\' decomp_type '/'...
%                         zerostr(i_ret,2) '_responses_' zerostr(i_sigma, 2) '.mat'], 'responses');
%                 end
% 
%             case {'mono'}
%                 min_wavelength = 4;
%                 [local_amp local_phase local_ori] = ...
%                     monogenic(ret, 5, min_wavelength, 2, 0.65, 1);
%                 responses = cat(4, local_amp, local_phase, local_ori);
%                 save([fibreroot 'filter_responses\' decomp_type '/'...
%                     zerostr(i_ret,2) '_responses_' zerostr(min_wavelength, 2) '.mat'], 'responses');
%         end
%     end
% end
% 
% %%
% decomp_types = cell(2, 2);
% decomp_types( 1,:) = {'g2d', 1};
% decomp_types( 2,:) = {'g2di', 1};
% 
% for i_decomp = 1:size(decomp_types,1)
%     decomp_type = decomp_types{i_decomp,1};
%     if decomp_types{i_decomp,2}
%         pred_dir = [fibreroot 'predictions\orientation\analytic\' decomp_type];
%         mkdir([pred_dir '\response']);
%         mkdir([pred_dir '\orientation']);
%         mkdir([pred_dir '\scale']);
%         mkdir([pred_dir '\normal_response']);
%     else
%         continue;
%     end
% 
%     for i_ret = 21:40
%         %load retinogram and merge RGB channels
%         ret = u_load([fibreroot '\images\' zerostr(i_ret,2) '_training.mat']);
%         ret = rgb2gray(ret);
% 
%         switch decomp_type
% 
%             case {'gabor'}
%                 %TO DO
% 
%             case {'g2d'}
%                 sigma_range = [1 2 4 8];
%                 raw_responses = compute_gaussian_2nd_derivatives(ret, sigma_range);
%                 [response, orientation, scale, normal_response] =...
%                     gaussian_2nd_derivative_line(raw_responses);
%                 
%                 save([pred_dir '\response\' zerostr(i_ret,2) '_training.mat'], 'response');
%                 save([pred_dir '\orientation\' zerostr(i_ret,2) '_training.mat'], 'orientation');
%                 save([pred_dir '\scale\' zerostr(i_ret,2) '_training.mat'], 'scale');
%                 save([pred_dir '\normal_response\' zerostr(i_ret,2) '_training.mat'], 'normal_response');
%                 
%             case {'g2di'}
%                 raw_responses = compute_gaussian_2nd_derivatives_d(ret, 1, 4);
%                 [response, orientation, scale, normal_response] =...
%                     gaussian_2nd_derivative_line(raw_responses);
%                 
%                 save([pred_dir '\response\' zerostr(i_ret,2) '_training.mat'], 'response');
%                 save([pred_dir '\orientation\' zerostr(i_ret,2) '_training.mat'], 'orientation');
%                 save([pred_dir '\scale\' zerostr(i_ret,2) '_training.mat'], 'scale');
%                 save([pred_dir '\normal_response\' zerostr(i_ret,2) '_training.mat'], 'normal_response');
%         end
% 
%     end
% end
% %%
% load([fibreroot 'orientations\gt\all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
% 
% decomp_types = cell(2, 2);
% decomp_types( 1,:) = {'g2d', 1};
% decomp_types( 2,:) = {'g2di', 1};
% 
% for i_decomp = 1:size(decomp_types,1)
%     decomp_type = decomp_types{i_decomp,1};
%     if decomp_types{i_decomp,2}
%         pred_dir = [fibreroot 'predictions\orientation\analytic\' decomp_type '\orientation'];
%         display(['Errors for ' decomp_type]); 
%         [ori_errors, ~, error_stats] =...
%             compute_image_orientation_errors(pred_dir, fg_mask_dir,...
%             'gt_orientations', gt_orientations, 'fov_mask_dir', fov_mask_dir);
%         display(error_stats);
% 
%         centre_stats = ori_error_stats(ori_errors(vessel_centres));
%         display(centre_stats);
%     end
% end
% %%
% %--------------------------------------------------------------------------
% % Testing Gabor filters...
% fibreroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];
% 
% load([fibreroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
% num_angles = [4 5 6 8 12 15 18 24 30 36 45 60];
% n_tests = length(num_angles);
% mean_abs_err = zeros(n_tests,1);
% median_abs_err = zeros(n_tests,1);
% results_dir = ['C:\isbe\asymmetry_project\experiments\DRIVE\gabor_orientation_prediction\'...
%     'sigma_all\'];
%         
% for i_angle = 1:n_tests
%     
%     rad_angles = pi/2 + pi*(0:(num_angles(i_angle)-1))'/num_angles(i_angle);
%     complex_angles = complex(cos(2*rad_angles), sin(2*rad_angles));
%     predicted_orientations = [];
%     
%     for i_ret = 21:40
%         display(['Analysing image ' num2str(i_ret)]);
%         
%         %load vessel mask
%         ret = rgb2gray(u_load([fibreroot 'images\' zerostr(i_ret,2) '_training.mat']));
%         v_mask = u_load([fibreroot 'vessel_masks\' zerostr(i_ret,2) '_training_v_mask.mat']);
%         f_mask = u_load([fibreroot 'fov_masks\' zerostr(i_ret,2) '_training_f_mask.mat']);
%         v_mask = v_mask & f_mask;
%         
%         for i_sigma = [1 2 4 8]
%             
%             %Compute responses
%             [responses] = compute_gabor_responses(ret, i_sigma, num_angles(i_angle));
%             
%             %Compute orientation as arg max of absolute response
%             [max_responses_i, max_angles_i] = max(abs(responses), [], 4);
%             clear responses;                
%             
%             if i_sigma == 1
%                 %For sigma = 1, create the max_responses and angles
%                 max_responses = max_responses_i;%(v_mask)
%                 max_angles = max_angles_i;%(v_mask)
%             else
%                 %For further scales, swap in the responses and angles if
%                 %they're greater than the existing responses               
%                 swap_idx = max_responses_i > max_responses;
%                 max_responses(swap_idx) = max_responses_i(swap_idx);
%                 max_angles(swap_idx) = max_angles_i(swap_idx);
%             end
%         end
%         predicted_ori_map = complex_angles(max_angles);
%         
%         if i_ret == 1           
%             figure; 
%             subplot(1,2,1); imgray(complex2rgb(predicted_ori_map .* abs(max_responses)));
%             subplot(1,2,2); imgray(complex2rgb(gt_ori_map));
%         end
%         
%         %Convert the max angles to an angle in complex form and add to the
%         %predictions list
%         predicted_orientations = [predicted_orientations; predicted_ori_map(v_mask)]; %#ok
%         
%     end
%     
%     %Now get the error stats
%     [orientation_errors, error_stats] = ori_error(predicted_orientations, gt_orientations);
%     mean_abs_err(i_angle) = error_stats.abs_mean;
%     median_abs_err(i_angle) = error_stats.abs_median;    
%     save([results_dir 'ori_errors_' zerostr(num_angles(i_angle),2) '.mat'], 'orientation_errors');
% end
% 
% figure; hold all;
% plot(num_angles, mean_abs_err, '-.x');
% plot(num_angles, median_abs_err, '-.x');
% legend({'Mean absolute errors', 'Median absolute errors'});
% title('Orientation errors for increasing angluar resolution of Gabor filter - DRIVE data');
% xlabel('Number of filters');
% %%
% fibreroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];
% 
% load([fibreroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
% num_angles = [4 5 6 8 12 15 18 24 30 36 45 60];
% n_tests = length(num_angles);
% mean_abs_err = zeros(n_tests,1);
% median_abs_err = zeros(n_tests,1);
% 
% for i_sigma = [1 2 4 8]
%     
%     results_dir = ['C:\isbe\asymmetry_project\experiments\DRIVE\gabor_orientation_prediction\'...
%         'sigma' num2str(i_sigma) '\'];
%     mkdir(results_dir);
%     
%     for i_angle = 1:n_tests
% 
%         rad_angles = pi/2 + pi*(0:(num_angles(i_angle)-1))'/num_angles(i_angle);
%         complex_angles = complex(cos(2*rad_angles), sin(2*rad_angles));
%         predicted_orientations = [];
% 
%         for i_ret = 21:40
%             display(['Analysing image ' num2str(i_ret)]);
% 
%             %load vessel mask
%             ret = rgb2gray(u_load([fibreroot 'images\' zerostr(i_ret,2) '_training.mat']));
%             v_mask = u_load([fibreroot 'vessel_masks\' zerostr(i_ret,2) '_training_v_mask.mat']);
%             f_mask = u_load([fibreroot 'fov_masks\' zerostr(i_ret,2) '_training_f_mask.mat']);
%             v_mask = v_mask & f_mask;
% 
%             %Compute responses
%             [responses] = compute_gabor_responses(ret, i_sigma, num_angles(i_angle));
%             
%             %Compute orientation as arg max of absolute response
%             [~, max_angles] = max(abs(responses), [], 4);
%             clear responses;               
%             predicted_ori_map = complex_angles(max_angles);
%         
%             %Convert the max angles to an angle in complex form and add to the
%             %predictions list
%             predicted_orientations = [predicted_orientations; predicted_ori_map(v_mask)]; %#ok
%         end
%         
%         %Now get the error stats
%         [orientation_errors, error_stats] = ori_error(predicted_orientations, gt_orientations);
%         mean_abs_err(i_angle) = error_stats.abs_mean;
%         median_abs_err(i_angle) = error_stats.abs_median;  
%         save([results_dir 'ori_errors_' zerostr(num_angles(i_angle),2) '.mat'], 'orientation_errors');
%         
%     end
%     
%     figure; hold all;
%     plot(num_angles, mean_abs_err, '-.x');
%     plot(num_angles, median_abs_err, '-.x');
%     legend({'Mean absolute errors', 'Median absolute errors'});
%     title('Orientation errors for increasing angluar resolution of Gabor filter - DRIVE data');
%     xlabel('Number of filters');
% end
% %%
% num_angles = [4 5 6 8 12 15 18 24 30 36 45 60];
% n_tests = length(num_angles);
% leg_text = cell(0,1);
% figure; a1 = gca; hold all; title('Mean absolute errors for varying angular resolution of Gabor filters');
% figure; a2 = gca; hold all; title('Median absolute errors for varying angular resolution of Gabor filters');
% figure; a3 = gca; hold all; title('CDF of absolute errors for varying angular resolution of Gabor filters'); 
% for i_sigma = [1 2 4 8 9]
%     
%     if i_sigma == 9
%         results_dir = ['C:\isbe\asymmetry_project\experiments\DRIVE\gabor_orientation_prediction\'...
%             'sigma_all\'];
%         leg_text{end+1} = '\sigma = [1 2 4 8]';
%     else
%         results_dir = ['C:\isbe\asymmetry_project\experiments\DRIVE\gabor_orientation_prediction\'...
%             'sigma' num2str(i_sigma) '\'];
%         leg_text{end+1} = ['\sigma = ' num2str(i_sigma)];
%     end
%     
%     mean_abs_err = zeros(n_tests,1);
%     median_abs_err = zeros(n_tests,1);
% 
%     for i_angle = 1:n_tests        
%         load([results_dir 'ori_errors_' zerostr(num_angles(i_angle),2) '.mat'], 'orientation_errors');
%         orientation_errors(isnan(orientation_errors)) = [];
%         mean_abs_err(i_angle) = mean(abs(orientation_errors))*180/pi;
%         median_abs_err(i_angle) = median(abs(orientation_errors))*180/pi;
%         
%         if num_angles(i_angle) == 30 
%             sorted_errors = sort(abs(orientation_errors))*180/pi;
%             cdf_pts = round(linspace(1, length(orientation_errors), 1001));
%             %plot(a3, sorted_errors(cdf_pts), 0:100, '-x');
%             plot(a3, sorted_errors(cdf_pts), 0:1000, '-x');
%         end
%         
%     end
%     plot(a1, num_angles, mean_abs_err, '-.x');
%     plot(a2, num_angles, median_abs_err, '-.x');
%     
% end
% legend(a1, leg_text);
% legend(a2, leg_text);
% legend(a3, leg_text);
% %%
% % Testing Gaussian 2nd derivative filters...
% fibreroot = [asymmetryroot,'data\retinograms\DRIVE_clean\training\'];
% load([fibreroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
% 
% results_dir = ['C:\isbe\asymmetry_project\experiments\DRIVE\g2d_orientation_prediction\'...
%         'sigma' num2str(i_sigma) '\'];
% 
% figure; hold all; a1 = gca;
% for i_sigma = 1:8
%     
%     predicted_orientations = [];
%     for i_ret = 21:40
%         %load retinogram and merge RGB channels
%         ret = u_load([fibreroot '\images\' zerostr(i_ret,2) '_training.mat']);
%         ret = rgb2gray(ret);
%         v_mask = u_load([fibreroot 'vessel_masks\' zerostr(i_ret,2) '_training_v_mask.mat']);
%         f_mask = u_load([fibreroot 'fov_masks\' zerostr(i_ret,2) '_training_f_mask.mat']);
%         v_mask = v_mask & f_mask;
% 
%         raw_responses = compute_gaussian_2nd_derivatives(ret, i_sigma);
%         [response, orientation] =...
%             gaussian_2nd_derivative_line(raw_responses);
%         
%         ori_map = complex(cos(2*orientation), sin(2*orientation));
%         predicted_orientations = [predicted_orientations; ori_map(v_mask)]; %#ok
%     end
%     
%     %Now get the error stats
%     [orientation_errors] = ori_error(predicted_orientations, gt_orientations);   
%     sorted_errors = sort(abs(orientation_errors))*180/pi;
%     cdf_pts = round(linspace(1, length(orientation_errors), 1001));
%     plot(a1, sorted_errors(cdf_pts), (0:1000)/1000, '-x');
%    
%         
%     results_dir = [pred_dir '\sigma_' num2str(i_sigma)];
%     mkdir(results_dir);
%     save([results_dir 'ori_errors_' zerostr(num_angles(i_angle),2) '.mat'], 'orientation_errors');
% end
% legend(a1, strcat('\sigma =  ', num2str((1:8)')));
% 
% %%
% fibreroot = [asymmetryroot 'data\retinograms\DRIVE_clean\test\predictions\orientation\'];
% 
% mkdir([fibreroot 'analytic\mono\orientations']);
% mkdir([fibreroot 'analytic\mono\responses']);
% mkdir([fibreroot 'analytic\mono\scales']);
% 
% for ii = 1:20
%     %load retinogram and merge RGB channels
%     ret = u_load([asymmetryroot('shared'),'data\retinograms\DRIVE\test\images_extended\' zerostr(ii,2) '_test.mat']);
%     ret = rgb2gray(ret);
%     
%     %Compute repsonses for mono and save
%     [response_map d ori_map scale_map] = monogenic_multiscale(ret, 4, 2, 2, 0.65);
%     save_uint8([fibreroot 'mono\analytic\orientations\' zerostr(ii,3) '_test_ori.mat'], ori_map);
%     save_uint8([fibreroot 'mono\analytic\responses\' zerostr(ii,3) '_test_response.mat'], response_map);
%     save_uint8([fibreroot 'mono\analytic\scales\' zerostr(ii,3) '_test_scale.mat'], scale_map);
% end
% %%
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% fibreroot = [asymmetryroot,'data\retinograms\STARE\training\'];
% 
% mkdir([fibreroot '/orientations/gt/']);
% gt_orientations = [];        
% vessel_centres = [];
% 
% for i_ret = 1:20
%     v_mask = u_load([fibreroot 'vessel_masks\' zerostr(i_ret,2) '_training_v_mask.mat']);
%     f_mask = u_load([fibreroot 'fov_masks\' zerostr(i_ret,2) '_training_f_mask.mat']);
%     v_mask = v_mask & f_mask;
%     vc_mask = bwmorph(v_mask, 'thin', inf);
%     
%     ori_map = u_load([fibreroot 'orientations\' zerostr(i_ret,2) '_training_ori.mat']);
%     
%     vessel_centres = [vessel_centres; vc_mask(v_mask)]; %#ok
%     gt_orientations = [gt_orientations; ori_map(v_mask)]; %#ok
% end
% save([fibreroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
% 
% %%
% 
% load([fibreroot '/orientations/gt/all_gt_orientations.mat'], 'gt_orientations', 'vessel_centres');
% num_angles = [4 5 6 8 12];
% n_tests = length(num_angles);
% mean_abs_err = zeros(n_tests,1);
% median_abs_err = zeros(n_tests,1);
% 
% for i_sigma = [1 2 4 8]
%     
%     results_dir = ['C:\isbe\asymmetry_project\experiments\STARE\gabor_orientation_prediction\'...
%         'sigma' num2str(i_sigma) '\'];
%     mkdir(results_dir);
%     
%     for i_angle = 1:n_tests
% 
%         rad_angles = pi/2 + pi*(0:(num_angles(i_angle)-1))'/num_angles(i_angle);
%         complex_angles = complex(cos(2*rad_angles), sin(2*rad_angles));
%         predicted_orientations = [];
% 
%         for i_ret = 1:20
%             display(['Analysing image ' num2str(i_ret)]);
% 
%             %load vessel mask
%             ret = rgb2gray(u_load([fibreroot 'images\' zerostr(i_ret,2) '_training.mat']));
%             v_mask = u_load([fibreroot 'vessel_masks\' zerostr(i_ret,2) '_training_v_mask.mat']);
%             f_mask = u_load([fibreroot 'fov_masks\' zerostr(i_ret,2) '_training_f_mask.mat']);
%             v_mask = v_mask & f_mask;
% 
%             %Compute responses
%             [responses] = compute_gabor_responses(ret, i_sigma, num_angles(i_angle));
%             
%             %Compute orientation as arg max of absolute response
%             [~, max_angles] = max(abs(responses), [], 4);
%             clear responses;               
%             predicted_ori_map = complex_angles(max_angles);
%         
%             %Convert the max angles to an angle in complex form and add to the
%             %predictions list
%             predicted_orientations = [predicted_orientations; predicted_ori_map(v_mask)]; %#ok
%         end
%         
%         %Now get the error stats
%         [orientation_errors, error_stats] = ori_error(predicted_orientations, gt_orientations);
%         mean_abs_err(i_angle) = error_stats.abs_mean;
%         median_abs_err(i_angle) = error_stats.abs_median;  
%         save([results_dir 'ori_errors_' zerostr(num_angles(i_angle),2) '.mat'], 'orientation_errors');
%         
%     end
%     
%     figure; hold all;
%     plot(num_angles, mean_abs_err, '-.x');
%     plot(num_angles, median_abs_err, '-.x');
%     legend({'Mean absolute errors', 'Median absolute errors'});
%     title('Orientation errors for increasing angluar resolution of Gabor filter - DRIVE data');
%     xlabel('Number of filters');
% end
% %%
% num_angles = [4 5 6 8 12];
% n_tests = length(num_angles);
% leg_text = cell(0,1);
% figure; a1 = gca; hold all; title('Mean absolute errors for varying angular resolution of Gabor filters');
% figure; a2 = gca; hold all; title('Median absolute errors for varying angular resolution of Gabor filters');
% figure; a3 = gca; hold all; title('CDF of absolute errors for varying angular resolution of Gabor filters'); 
% for i_sigma = [1 2 4 8 9]
%     
%     if i_sigma == 9
%         results_dir = ['C:\isbe\asymmetry_project\experiments\STARE\gabor_orientation_prediction\'...
%             'sigma_all\'];
%         leg_text{end+1} = '\sigma = [1 2 4 8]';
%     else
%         results_dir = ['C:\isbe\asymmetry_project\experiments\STARE\gabor_orientation_prediction\'...
%             'sigma' num2str(i_sigma) '\'];
%         leg_text{end+1} = ['\sigma = ' num2str(i_sigma)];
%     end
%     
%     mean_abs_err = zeros(n_tests,1);
%     median_abs_err = zeros(n_tests,1);
% 
%     for i_angle = 1:n_tests        
%         load([results_dir 'ori_errors_' zerostr(num_angles(i_angle),2) '.mat'], 'orientation_errors');
%         orientation_errors(isnan(orientation_errors)) = [];
%         mean_abs_err(i_angle) = mean(abs(orientation_errors))*180/pi;
%         median_abs_err(i_angle) = median(abs(orientation_errors))*180/pi;
%         
%         if num_angles(i_angle) == 30 
%             sorted_errors = sort(abs(orientation_errors))*180/pi;
%             cdf_pts = round(linspace(1, length(orientation_errors), 1001));
%             %plot(a3, sorted_errors(cdf_pts), 0:100, '-x');
%             plot(a3, sorted_errors(cdf_pts), 0:1000, '-x');
%         end
%         
%     end
%     plot(a1, num_angles, mean_abs_err, '-.x');
%     plot(a2, num_angles, median_abs_err, '-.x');
%     
% end
% legend(a1, leg_text);
% legend(a2, leg_text);
% legend(a3, leg_text);