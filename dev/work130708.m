warning('off', 'ori_error:nans');
pred_dir = 'C:\isbe\nailfold\data\rsa_study\set2\predictions\orientation\rf_regression\';
label_dir = 'C:\isbe\nailfold\data\rsa_study\set2\orientations\';
fov_mask_dir = 'C:\isbe\nailfold\data\rsa_study\set2\fov_masks\';
fg_mask_dir = 'C:\isbe\nailfold\data\rsa_study\set2\vessel_masks\';
fgc_mask_dir = 'C:\isbe\nailfold\data\rsa_study\set2\vessel_centre_masks\';
%%
rf_codes = cell(4, 2);
rf_codes( 1,:) = {'219121', 'OrigIms'};
rf_codes( 2,:) = {'219122', 'NormIms'};
rf_codes( 3,:) = {'219133', 'OrigImsAllVess'};
rf_codes( 4,:) = {'222210', 'RandBG'};
rf_codes( 5,:) = {'222835', 'RandBG?'};

error_medians = zeros(4,2);
error_means = zeros(4,2);
%%
for ii = 5%1:5

    if 0%exist([pred_dir rf_codes{ii,1} '\errors\error_stats.mat'], 'file')
        load([pred_dir rf_codes{ii,1} '\errors\error_stats.mat']);
    else
        [orientation_errors, ~, error_stats] =...
            compute_image_orientation_errors([pred_dir rf_codes{ii,1} '\'], fgc_mask_dir,...
            'label_dir', label_dir, 'fov_mask_dir', fov_mask_dir,...
            'save_dir', [pred_dir rf_codes{ii,1} '\errors\']);
    end
    %centre_stats = ori_error_stats(orientation_errors(vessel_centres));
    
    error_medians(ii,1) = error_stats.abs_median;
    %error_medians(ii,2) = centre_stats.abs_median;
    
    error_means(ii,1) = error_stats.abs_mean;
    %error_means(ii,2) = centre_stats.abs_mean;
    
    display(['Errors for ' rf_codes{ii,2} ...
        ', abs mean: ' num2str(error_stats.abs_mean,3) ', abs_median: ' num2str(error_stats.abs_median,3)]); 
       
end
%%
pred_dir = 'C:\isbe\nailfold\data\rsa_study\set2\predictions\detection\rf_classification\';
label_dir = 'C:\isbe\nailfold\data\rsa_study\set2\vessel_centre_masks\';
rf_codes = cell(2, 2);
rf_codes( 1,:) = {'219128', 'NormIms'};
rf_codes( 2,:) = {'219131', 'OrigIms'};
rf_codes( 3,:) = {'182321', 'OrigIms'};
rf_codes( 4,:) = {'182321', 'OrigIms'};

auc_allc = [];
roc_pts_allc = [];
%
for ii = 1:3
    if 0%exist([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'], 'file')
        load([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'])
    else
        [roc_pts, auc, auc_ind] = compute_roc_image_set([pred_dir rf_codes{ii,1} '\'], label_dir, 'fov_mask_dir', fov_mask_dir);
        %create_folder([pred_dir rf_codes{ii,1} '\roc\']);
        %save([pred_dir rf_codes{ii,1} '\roc\roc_stats.mat'], 'roc_pts', 'auc', 'auc_ind');
    end
    roc_pts_allc(:,:,ii) = roc_pts; %#ok
    auc_allc(ii,:) = [auc auc_ind']; %#ok
end
%%
%Define constants
vessel_prob_smoothing_sigma = 2;
strong_vessel_thresh = 0.25;
g = gaussian_filters_1d(vessel_prob_smoothing_sigma);
g = g / sum(g);
do_plot = 1;
rsa_dir = 'rsa_study';
max_length = 1000;

i_sz = {'normal', 'enlarged'};
vessel_size = i_sz{1};


%Define directories
prob_dir = [nailfoldroot 'data/' rsa_dir '/set2/predictions/detection/rf_classification/219131/'];
ori_dir = [nailfoldroot 'data/' rsa_dir '/set2/predictions/orientation/rf_regression/222210/'];
ori_truth_dir = [nailfoldroot 'data/' rsa_dir '/set2/orientations/' vessel_size];

apex_dir = [nailfoldroot 'data/' rsa_dir '/apexes/' vessel_size '/'];
contour_dir = [nailfoldroot 'data/' rsa_dir '/vessel_contours/' vessel_size '/'];

%Get list of vessel names
v_list = dir([prob_dir, vessel_size '*.mat']);
%%    
for i_ve = 89%81:102;%:length(v_list) 
        
    %load vessel
    apex_name = v_list(i_ve).name(length(vessel_size)+(1:8));
    apex_struc = load([apex_dir apex_name '.mat']);
    vessel_struc = u_load([apex_dir apex_name '_vessel.mat']);
    contour_struc = load([contour_dir apex_name '_vessel_contour.mat']);

    %Load in the prob and ori maps
    prob_patch = u_load([prob_dir v_list(i_ve).name]);
    prob_patch = conv2(g', g, prob_patch, 'same');
    ori_patch = u_load([ori_dir v_list(i_ve).name]);
    ori_patch = conv2(g', g, ori_patch, 'same');
    ori_truth = u_load([ori_truth_dir apex_name '_vessel_ori.mat']);



    display(['Processing ' apex_name ' from ' vessel_size]);     

    vessel_nms = mb_non_maximal_supp(prob_patch, angle(ori_patch)/2);
    strong_vessels = vessel_nms > strong_vessel_thresh;
    if any(strong_vessels(:))
        [rstrong cstrong] = find(strong_vessels);
        vessel_centre = bwselect(vessel_nms > 0, cstrong, rstrong, 8);
    else
        vessel_centre = strong_vessels;
    end
    [vcy vcx] = find(vessel_centre);

    %Correct apex coordinates frame
    apex_xy(:,1) = apex_struc.apex_xy(:,1) +...
        apex_struc.apex_properties.sc - vessel_struc.vessel_properties.sc;
    apex_xy(:,2) = apex_struc.apex_xy(:,2) +...
        apex_struc.apex_properties.sr - vessel_struc.vessel_properties.sr;

    %Resample the vessel pts to be equalled spaced, 1 pixel apart
    v_pts = spline_contour(contour_struc.vessel_centre, [], 1);
    dists = sum(bsxfun(@minus, v_pts, mean(apex_xy)).^2, 2);
    [~, apex_idx] = min(dists);

    [shape_contexts_v, path_contexts_v, path_map_v] = shape_context_prob_track_mult(...
                    ori_patch.*prob_patch, prob_patch, ori_patch,...
                    round(v_pts(apex_idx,2)), round(v_pts(apex_idx,1)),...
                    'num_streams', 1e4, 'stopping_prob', 0.2);
                
    [shape_contexts_r, path_contexts_r, path_map_r] = shape_context_prob_track_mult(...
                    ori_patch.*prob_patch, prob_patch, ori_patch,...
                    round(v_pts(1,2)), round(v_pts(1,1)),...
                    'num_streams', 1e4, 'stopping_prob', 0.2);

%     [~,~, path_map_a] = shape_context_prob_track_mult(...
%                     ori_patch.*prob_patch, prob_patch, ori_patch,...
%                     round(mean(apex_xy(:,2)))+5, round(mean(apex_xy(:,1))),...
%                     'num_streams', 1e4, 'stopping_prob', 0.2);

%     [~,~, path_map_t] = shape_context_prob_track_mult(...
%                     ori_patch.*prob_patch, prob_patch, ori_patch,...
%                     54, 52,...
%                     'num_streams', 1e4, 'stopping_prob', 0.2);

    v_im = double(vessel_struc.vessel_patch);
    v_im = (v_im - min(v_im(:))) / (max(v_im(:))-min(v_im(:)));
    v_im = cat(3, v_im, v_im, v_im)*.5;
    v_im(:,:,1) = v_im(:,:,1) + (.5*path_map_v/max(path_map_v(:)));
 
    figure; 
    subplot(2,3,1); imgray(display_shape_context(shape_contexts_v, 10, 'p'));
    subplot(2,3,2); imgray(display_shape_context(shape_contexts_v, 10, 'm'));
    subplot(2,3,4); imgray(display_shape_context(shape_contexts_r, 10, 'p'));
    subplot(2,3,5); imgray(display_shape_context(shape_contexts_r, 10, 'm'));
    
    subplot(2,3,[3 6]); imgray(v_im);
    plot(v_pts(:,1), v_pts(:,2), 'g--', 'linewidth', 1);
    
%     figure; 
%     
%     subplot(1,3,1); imgray(path_map_v);
%     plot(v_pts(:,1), v_pts(:,2), 'k', 'linewidth', 1.5);
    
    %figure; imgray(path_map_a);
    %plot(v_pts(:,1), v_pts(:,2), 'g');
    
    %figure; imgray(path_map_t);
    %plot(v_pts(:,1), v_pts(:,2), 'g');

%     theta_patch = angle(ori_patch) / 2;
%     figure; imgray(prob_patch);
%     quiver(cos(theta_patch).*abs(ori_patch), -sin(theta_patch).*abs(ori_patch), 'r');
%     quiver(-cos(theta_patch).*abs(ori_patch), sin(theta_patch).*abs(ori_patch), 'r');
%     plot(v_pts(:,1), v_pts(:,2), 'g');
% 
%     vessel_curv = abs(complex_gaussian_curvature(ori_patch, 2));
%     figure; imgray(vessel_curv);
%     plot(v_pts(:,1), v_pts(:,2), 'g');

%     subplot(1,3,2); imgray(complex2rgb(ori_patch));
%     
%     subplot(1,3,3); imgray(abs(ori_error(ori_truth, ori_patch))); colormap jet; %colorbar;
%     plot(v_pts(:,1), v_pts(:,2), 'w', 'linewidth', 1.5);
end
%%
for i_ve = 21:40
    ori_patch = u_load([ori_dir v_list(i_ve).name]);
    figure; imgray(complex2rgb(ori_patch));
end
%%
job_args.sampling_args.sampling_method ='generate_training_data';
job_args.sampling_args.image_root ='C:/isbe/nailfold/data/rsa_study/set1';
job_args.sampling_args.image_type = 'real';
job_args.sampling_args.num_samples = 5000;
job_args.sampling_args.max_n_images = [];
job_args.sampling_args.shift_images = 0;
job_args.sampling_args.shrink_fov = 0;
job_args.sampling_args.win_size = 1;
job_args.sampling_args.replace_sample = 0;
job_args.sampling_args.image_dir = 'C:/isbe/nailfold/data/rsa_study/set1/images';
job_args.sampling_args.fov_mask_dir = 'C:/isbe/nailfold/data/rsa_study/set1/fov_masks';
job_args.sampling_args.fg_mask_dir = 'C:/isbe/nailfold/data/rsa_study/set1/vessel_masks';
job_args.sampling_args.probability_dir = [];
job_args.sampling_args.output_type = 'orientation';
job_args.sampling_args.bg_ratio = 0.5;
job_args.sampling_args.ori_dir = 'C:/isbe/nailfold/data/rsa_study/set1/orientations';
job_args.sampling_args.width_dir = [];
job_args.sampling_args.sampled_data_dir = [];
    
job_args.decomposition_args.decomp_type = {'g2da'; 'h2da'};
job_args.decomposition_args.win_size = 1;
job_args.decomposition_args.rgb_channel = 'rgb';
job_args.decomposition_args.normalise = 0;
job_args.decomposition_args.sigma_range = 2;
job_args.decomposition_args.num_angles = 3;
job_args.decomposition_args.do_max = 0;
job_args.decomposition_args.rotate = 0;
job_args.decomposition_args.pca = [];

warning('off', 'ASYM:unexpectedArgument');
[training_data training_labels] = generate_training_data(job_args);
%%
theta_patch = angle(ori_patch) / 2;
x1 = round(mean(apex_xy(:,1)));
y1 = round(mean(apex_xy(:,2)));
vx1 = cos(theta_patch(y1,x1));
vy1 = -sin(theta_patch(y1,x1));

% [particles_x particles_y] = jetstream(prob_patch, theta_patch, [], [x1 y1], -[vx1 vy1], ...
%     'd', 1,...
%     'nu', 0.01,...
%     'M', 100,...
%     'sigma_theta', 0.1,...
%     'sigma_psi', 1,...
%     'lambda', [],...
%     'N', 200, ...
%     'plot', 1);
% [particles_x particles_y] = jetstream(prob_patch, theta_patch, [], [x1 y1], [vx1 vy1], ...
%     'd', 1,...
%     'nu', 0.01,...
%     'M', 100,...
%     'sigma_theta', 0.1,...
%     'sigma_psi', 1,...
%     'lambda', [],...
%     'N', 200, ...
%     'plot', 1);
[particles_x particles_y] = jetstream_rf(prob_patch, ori_patch, [], [x1 y1], [vx1 vy1], ...
    'd', 1,...
    'M', 100,...
    'step_length', 1,...
    'N', 200, ...
    'double_angle', 1,...
    'plot', 1);
[particles_x particles_y] = jetstream_rf(prob_patch, ori_patch, [], [x1 y1], -[vx1 vy1], ...
    'd', 1,...
    'M', 100,...
    'step_length', 1,...
    'N', 200, ...
    'double_angle', 1,...
    'plot', 1);