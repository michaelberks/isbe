for i_pt = 1:size(v_pts,1)-2
            
    if angle(ori_diff(i_pt)) > .5
        plot(gca, ...
            v_pts(i_pt+1,1)+[-1 1]*cos(angle(predicted_ori(i_pt))/2),... %cos(marked_ori_change(i_pt)/2) + 
            v_pts(i_pt+1,2)-[-1 1]*sin(angle(predicted_ori(i_pt))/2), 'r');
    elseif angle(ori_diff(i_pt)) > -.5
        plot(gca, ...
            v_pts(i_pt+1,1)+[-1 1]*cos(angle(predicted_ori(i_pt))/2),... %cos(marked_ori_change(i_pt)/2) + 
            v_pts(i_pt+1,2)-[-1 1]*sin(angle(predicted_ori(i_pt))/2), 'y');
    else
        plot(gca, ...
            v_pts(i_pt+1,1)+[-1 1]*cos(angle(predicted_ori(i_pt))/2),... %cos(marked_ori_change(i_pt)/2) + 
            v_pts(i_pt+1,2)-[-1 1]*sin(angle(predicted_ori(i_pt))/2), 'g');
    end
end

DATA_ROOT="scratch/nailfold/" MODEL_ROOT="models/vessel/orientation/rf_regression" MODEL_PATH="'222835'" NUM_JOBS=32 IMAGE_ROOT="data/set12" USE_SAMPLED_MAPS=1 MASK_DIR="fov_masks" qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
DATA_ROOT="scratch/nailfold/" MODEL_ROOT="models/vessel/detection/rf_classification" MODEL_PATH="'222836'" NUM_JOBS=32 IMAGE_ROOT="data/set12" USE_SAMPLED_MAPS=1 MASK_DIR="fov_masks" qsub -V -t 1-32 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh

%%
apex_radius = 5;
max_apex_guess = 500;
vessel_prob_smoothing_sigma = 2;
curvature_smoothing_sigma = 2;
strong_vessel_thresh = 0.25;
curv_max = 0.5;

g = gaussian_filters_1d(vessel_prob_smoothing_sigma);
g = g / sum(g);
neg_curv_vessels = find(v_stats.vessel_curv.apex_max < 0);

for i_ve = 1:20%length(neg_curv_vessels);
    curr_ves = neg_curv_vessels(i_ve);
    
    i_im = v_stats.vessel_properties.im_num(curr_ves);
    im_name = pred_list(i_im).name(1:end-9);
    vessel_im = u_load([image_dir im_name '.mat']);
    vessel_im = sample_window(vessel_im, 401,...
        round(v_stats.vessel_properties.apex_pos(curr_ves,2)),...
        round(v_stats.vessel_properties.apex_pos(curr_ves,1)), 150);
    
    vessel_ori = u_load([ori_dir pred_list(i_im).name]);
    vessel_ori = sample_window(vessel_ori, 401,...
        round(v_stats.vessel_properties.apex_pos(curr_ves,2)),...
        round(v_stats.vessel_properties.apex_pos(curr_ves,1)));
    
    vessel_ori = conv2(g', g, vessel_ori, 'same');
    vessel_curv = -(complex_gaussian_curvature(vessel_ori ./ (abs(vessel_ori)+1e-6), curvature_smoothing_sigma));
    
    figure; 
    subplot(1,3,1); imgray(vessel_im);
    plot(v_stats.vessel_properties.apex_pos(curr_ves,1), v_stats.vessel_properties.apex_pos(curr_ves,2), 'rx');
    subplot(1,3,2); imgray(complex2rgb(vessel_ori));
    plot(v_stats.vessel_properties.apex_pos(curr_ves,1), v_stats.vessel_properties.apex_pos(curr_ves,2), 'rx');
    subplot(1,3,3); imgray(vessel_curv);
    plot(v_stats.vessel_properties.apex_pos(curr_ves,1), v_stats.vessel_properties.apex_pos(curr_ves,2), 'rx');
end
%%
rsa_dir = 'rsa_study';
vessel_prob_smoothing_sigma = 2;
g = gaussian_filters_1d(vessel_prob_smoothing_sigma);
g = g / sum(g);

%Define directories
prob_dir = [nailfoldroot 'data/' rsa_dir '/set12/predictions/detection/rf_classification/222836/'];
ori_dir = [nailfoldroot 'data/' rsa_dir '/set12/predictions/orientation/rf_regression/222835/'];

ngrid_rows = 5;
ngrid_cols = 20;
grid_row_spacing = 5;
grid_col_spacing = 5;
grid_yoffset = 20;
grid_x = linspace(-grid_col_spacing*ngrid_cols/2, grid_col_spacing*ngrid_cols/2, ngrid_cols);
grid_x = repmat(grid_x, ngrid_rows, 1);

grid_y = grid_yoffset + linspace(-grid_row_spacing*ngrid_rows/2, grid_row_spacing*ngrid_rows/2, ngrid_rows);
grid_y = repmat(grid_y', 1, ngrid_cols);


%Loop through each vessel
apex_angles = zeros(0,1);
angle_thresh = 0.9;
for i_sz = {'normal', 'enlarged'}
    
    vessel_size = i_sz{1};
    apex_dir = [nailfoldroot 'data/' rsa_dir '/apexes/' vessel_size '/'];
    contour_dir = [nailfoldroot 'data/' rsa_dir '/vessel_contours/' vessel_size '/'];

    %Get list of vessel names
    v_list = dir([ori_dir, vessel_size '*.mat']);
    
    for i_ve = 25:36%1:length(v_list)%1:20%
        plot_num = rem(i_ve - 1, 12) + 1;
        if plot_num == 1
            figure;
        end

        %load vessel
        apex_name = v_list(i_ve).name(length(vessel_size)+(1:8));
        contour_struc = load([contour_dir apex_name '_vessel_contour.mat']);
        
        display(['Processing ' apex_name ' from ' vessel_size]);

        %Resample the vessel pts to be equalled spaced, 1 pixel apart
        v_pts = spline_contour(contour_struc.vessel_centre, [], 1);

        apex_struc = load([apex_dir apex_name '.mat']);
        vessel_struc = u_load([apex_dir apex_name '_vessel.mat']);

        %Correct apex coordinates frame
        apex_xy(:,1) = apex_struc.apex_xy(:,1) +...
            apex_struc.apex_properties.sc - vessel_struc.vessel_properties.sc;
        apex_xy(:,2) = apex_struc.apex_xy(:,2) +...
            apex_struc.apex_properties.sr - vessel_struc.vessel_properties.sr;

        %Find the vessel pt nearets the apex        
        dists = sum(bsxfun(@minus, v_pts, mean(apex_xy)).^2, 2);
        [~, apex_idx] = min(dists);
        
        apex_vec = v_pts(apex_idx-1,:) - v_pts(apex_idx+1,:);
        apex_angle = atan2(-apex_vec(2), apex_vec(1));
        if (apex_angle > -pi/2) && (apex_angle < pi/2);
            v_pts = flipud(v_pts);
            apex_idx = size(v_pts,1) - apex_idx + 1;
            apex_vec = v_pts(apex_idx-1,:) - v_pts(apex_idx+1,:);
            apex_angle = atan2(-apex_vec(2), apex_vec(1));
        end
        apex_angles(end+1) = apex_angle;
        
        xy_steps = diff(v_pts);
        marked_ori = exp(2i*atan(-xy_steps(:,2) ./ xy_steps(:,1)));
        
        left_tri = find(abs(ori_error(marked_ori(1:apex_idx-1), marked_ori(apex_idx))) > (angle_thresh*pi/2), 1, 'last');
        right_tri = find(abs(ori_error(marked_ori(apex_idx+1:end), marked_ori(apex_idx))) > (angle_thresh*pi/2), 1, 'first');
        
        if ~isempty(right_tri)
            right_tri = right_tri + apex_idx;
        end
        
        apex_theta = angle(marked_ori(apex_idx))/2;
        rot_mat = [cos(apex_theta) -sin(apex_theta); sin(apex_theta) cos(apex_theta)];

        grid_xy_i = bsxfun(@plus, [grid_x(:) grid_y(:)]*rot_mat, v_pts(apex_idx,:));

        prob_patch = u_load([prob_dir v_list(i_ve).name]);
        prob_patch = conv2(g', g, prob_patch, 'same');
        ori_patch = u_load([ori_dir v_list(i_ve).name]);
        ori_patch = conv2(g', g, ori_patch, 'same');
        
        sampled_prob_patch = interp2(prob_patch, grid_xy_i(:,1), grid_xy_i(:,2),...
            'nearest');
        sampled_prob_patch = reshape(sampled_prob_patch, ngrid_rows, ngrid_cols);
        
        sampled_ori_patch = interp2(ori_patch, grid_xy_i(:,1), grid_xy_i(:,2),...
            'nearest');
        sampled_ori_patch = reshape(sampled_ori_patch, ngrid_rows, ngrid_cols);
        sampled_ori_patch = sampled_ori_patch * conj(marked_ori(apex_idx));
        
        poss_points_ori = ...
            (angle(sampled_ori_patch)/2 > angle_thresh*pi/2);
            
        poss_points_prob = ...
            (sampled_prob_patch > 0.5);
            
        [poss_rp poss_cp] = find(poss_points_prob);
        [poss_ro poss_co] = find(poss_points_ori);
        
        subplot(3, 4, plot_num);
        
%         imgray(complex2rgb(sampled_ori_patch));
%         imgray(sampled_prob_patch);
%         plot(poss_cp, poss_rp, 'r.');
%         plot(poss_co, poss_ro, 'g.');

%         imgray(prob_patch);
%         plot(v_pts(:,1), v_pts(:,2));
%         plot(v_pts(apex_idx,1), v_pts(apex_idx,2), 'rx');
% 
%         plot(grid_xy_i(:,1), grid_xy_i(:,2), 'rx')
%         
%         axis equal ij; hold on;
        imgray(complex2rgb(ori_patch));
        plot(v_pts(:, 1), v_pts(:, 2));
        plot(v_pts(apex_idx, 1), v_pts(apex_idx, 2), 'mo');
        
        if ~isempty(left_tri) && ~isempty(right_tri)
            plot(...
                v_pts([apex_idx; left_tri; right_tri; apex_idx], 1),...
                v_pts([apex_idx; left_tri; right_tri; apex_idx], 2),...
                'r', 'linewidth', 2);
        end                   
        
    end
            


end
apex_degrees = 180*(apex_angles)/pi;
apex_degrees(apex_degrees <0) = apex_degrees(apex_degrees<0) + 360;

%%
btf = [62 109 119 151 176];
vessel_size = 'normal';
apex_dir = [nailfoldroot 'data/' rsa_dir '/apexes/' vessel_size '/'];
contour_dir = [nailfoldroot 'data/' rsa_dir '/vessel_contours/' vessel_size '/'];

%Get list of vessel names
v_list = dir([ori_dir, vessel_size '*.mat']);

for i_ve = btf
    
    %load vessel
    apex_name = v_list(i_ve).name(length(vessel_size)+(1:8));
    contour_struc = load([contour_dir apex_name '_vessel_contour.mat']);

    display(['Processing ' apex_name ' from ' vessel_size]);

    %Resample the vessel pts to be equalled spaced, 1 pixel apart
    v_pts = spline_contour(contour_struc.vessel_centre, [], 1);

    apex_struc = load([apex_dir apex_name '.mat']);
    vessel_struc = u_load([apex_dir apex_name '_vessel.mat']);

    %Correct apex coordinates frame
    apex_xy(:,1) = apex_struc.apex_xy(:,1) +...
        apex_struc.apex_properties.sc - vessel_struc.vessel_properties.sc;
    apex_xy(:,2) = apex_struc.apex_xy(:,2) +...
        apex_struc.apex_properties.sr - vessel_struc.vessel_properties.sr;

    %Find the vessel pt nearets the apex        
    dists = sum(bsxfun(@minus, v_pts, mean(apex_xy)).^2, 2);
    [~, apex_idx] = min(dists);
        
    figure; axis equal ij; hold on;
    plot(v_pts(:, 1), v_pts(:, 2));
    plot(v_pts(1, 1), v_pts(1, 2), 'rx');
    plot(v_pts(apex_idx, 1), v_pts(apex_idx, 2), 'mo');
    plot(v_pts(end, 1), v_pts(end, 2), 'gx');
end
%%
image_dir = 'C:\isbe\nailfold\data\rsa_study\test\images\';
vessel_im = u_load([image_dir '10147c.mat']);
vessel_patch = vessel_im(200:500, 300:500);
figure; imgray(vessel_im); caxis([-20 10]);
figure; imgray(vessel_patch);
%%
vessel_patch = u_load('C:\isbe\nailfold\data\rsa_study\set1\images\normalapex0008_vessel.mat');
vessel_patch = u_load('C:\isbe\nailfold\data\rsa_study\set2\images\normalapex1871_vessel.mat');
rf = u_load('C:\isbe\nailfold\models\vessel\orientation\rf_regression\222835\predictor.mat');
rf.tree_root = 'C:\isbe\nailfold\models\vessel\orientation\rf_regression/';
load('C:\isbe\nailfold\models\vessel\orientation\rf_regression\222835\job_args.mat');

x = 69;
y = 154;

interest_mask = false(size(vessel_patch));
interest_mask(y,x) = 1;
[patch_ori] = predict_image(...
    'image_in', vessel_patch,...
    'decomposition_args', job_args.decomposition_args,...
    'predictor', rf, ...
    'prediction_type', 'rf_regression',...
    'output_type', 'orientation',...
    'use_probs', 0,...
    'mask', interest_mask,...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 128,...
    'incremental_results', 0);
figure; imgray(complex2rgb(patch_ori));
%%
vessel_prob_smoothing_sigma = 2;

g = gaussian_filters_1d(vessel_prob_smoothing_sigma);
g = g / sum(g);
rsa_dir = 'rsa_study';

model_dir = [nailfoldroot 'data/' rsa_dir '/models/apex_templates/'];
ori_dir = [nailfoldroot 'data/' rsa_dir '/set12/predictions/orientation/rf_regression/222835/'];
prob_dir = [nailfoldroot 'data/' rsa_dir '/set12/predictions/detection/rf_classification/222836/'];
vessel_size = 'normal';
load([model_dir 'apex_templates.mat']);

%Get list of vessel names
v_list = dir([ori_dir, vessel_size '*.mat']);

prior_left = bsxfun(@rdivide, curv_stats.l.marked_curv_by_step.hist, sum(curv_stats.l.marked_curv_by_step.hist,2));
prior_right = bsxfun(@rdivide, curv_stats.r.marked_curv_by_step.hist, sum(curv_stats.r.marked_curv_by_step.hist,2));

ori_hist_r_smoothed = imfilter(curv_stats.r.pred_ori_by_step.hist, conv(g,g), 'circular');
ori_hist_r_smoothed = imfilter(ori_hist_r_smoothed, conv(g,g)', 'replicate');

ori_hist_l_smoothed = imfilter(curv_stats.l.pred_ori_by_step.hist, conv(g,g), 'circular');
ori_hist_l_smoothed = imfilter(ori_hist_l_smoothed, conv(g,g)', 'replicate');

apex_prior = cat(3, ori_hist_l_smoothed, ori_hist_r_smoothed);
apex_prior_scores_all = zeros(0,1);
smoo_prior_scores_all = zeros(0,1);
rand_prior_scores_all = zeros(0,1);
apex_corr_scores_all = zeros(0,1);
apex_corr_dists_all = zeros(0,1);
rand_idxs = zeros(0,1);

num_streams = 1e3;
%
for i_sz = {'normal', 'enlarged'}
    
    vessel_size = i_sz{1};
    apex_dir = [nailfoldroot 'data/' rsa_dir '/apexes/' vessel_size '/'];
    contour_dir = [nailfoldroot 'data/' rsa_dir '/vessel_contours/' vessel_size '/'];

    %Get list of vessel names
    v_list = dir([ori_dir, vessel_size '*.mat']);
    %
    for i_ve = 1:length(v_list)%1:20%
    
        %load vessel
        apex_name = v_list(i_ve).name(length(vessel_size)+(1:8));
        contour_struc = load([contour_dir apex_name '_vessel_contour.mat']);
        apex_struc = load([apex_dir apex_name '.mat']);
        vessel_struc = u_load([apex_dir apex_name '_vessel.mat']);
        
        display(['Processing ' apex_name ' from ' vessel_size]);

        prob_patch = u_load([prob_dir v_list(i_ve).name]);
        prob_patch = conv2(g', g, prob_patch, 'same');
        ori_patch = u_load([ori_dir v_list(i_ve).name]);
        %ori_patch = conv2(g', g, ori_patch, 'same');        

        %Resample the vessel pts to be equalled spaced, 1 pixel apart
        v_pts = spline_contour(contour_struc.vessel_centre, [], 1);


        %Correct apex coordinates frame
        apex_xy(:,1) = apex_struc.apex_xy(:,1) +...
            apex_struc.apex_properties.sc - vessel_struc.vessel_properties.sc;
        apex_xy(:,2) = apex_struc.apex_xy(:,2) +...
            apex_struc.apex_properties.sr - vessel_struc.vessel_properties.sr;

        %Find the vessel pt nearets the apex        
        dists = sum(bsxfun(@minus, v_pts, mean(apex_xy)).^2, 2);
        [~, apex_idx] = min(dists);


    %     [~,~, path_map] = shape_context_prob_track_mult(...
    %                     ori_patch.*prob_patch, prob_patch, ori_patch,...
    %                     round(v_pts(apex_idx,2)), round(v_pts(apex_idx,1)),...
    %                     'num_streams', 1e4, 'stopping_prob', 0.1);
    %     
    %     [~,~, path_map_prior] = shape_context_prob_track_mult(...
    %                 ori_patch.*prob_patch, prob_patch, ori_patch,...
    %                 round(v_pts(apex_idx,2)), round(v_pts(apex_idx,1)),...
    %                 'num_streams', 1e4, 'stopping_prob', 0.1,...
    %                 'prior_left', prior_left,...
    %                 'prior_right', prior_right);            
    
        rand_idx = ceil(size(v_pts,1)*rand);
        rand_idxs(end+1) = rand_idx; %#ok
        
        if i_ve > 20 && i_ve <= 30
            [apex_prior_scores, apex_path_map] = apex_prior_prob_track_mult(...
                        prob_patch, ori_patch,...
                        round(v_pts(apex_idx,2)), round(v_pts(apex_idx,1)),...
                        'num_streams', num_streams, 'stopping_prob', 0.1,...
                        'apex_prior', apex_prior);

    %         [smoo_prior_scores, smoo_path_map] = apex_prior_prob_track_mult(...
    %                     prob_patch, conv2(g', g, ori_patch, 'same'),...
    %                     round(v_pts(apex_idx,2)), round(v_pts(apex_idx,1)),...
    %                     'num_streams', num_streams, 'stopping_prob', 0.1,...
    %                     'apex_prior', apex_prior);

            [rand_prior_scores, rand_path_map] = apex_prior_prob_track_mult(...
                        prob_patch, ori_patch,...
                        round(v_pts(rand_idx,2)), round(v_pts(rand_idx,1)),...
                        'num_streams', num_streams, 'stopping_prob', 0.1,...
                        'apex_prior', apex_prior);
        else

            [apex_prior_scores] = apex_prior_prob_track_mult(...
                        prob_patch, ori_patch,...
                        round(v_pts(apex_idx,2)), round(v_pts(apex_idx,1)),...
                        'num_streams', num_streams, 'stopping_prob', 0.1,...
                        'apex_prior', apex_prior);

            [rand_prior_scores] = apex_prior_prob_track_mult(...
                        prob_patch, ori_patch,...
                        round(v_pts(rand_idx,2)), round(v_pts(rand_idx,1)),...
                        'num_streams', num_streams, 'stopping_prob', 0.1,...
                        'apex_prior', apex_prior);
        end
        
        [maxima_pos, maxima_vals] = template_match_apexes(vessel_struc.vessel_patch, ...
            {'g2d', apex_template.g2}, 'threshold', -1);
        
        dists = sum(bsxfun(@minus, maxima_pos, mean(apex_xy)).^2, 2);
        [min_dist, min_idx] = min(dists);
        apex_corr_scores_all(end+1,:) = maxima_vals(min_idx); %#ok
        apex_corr_dists_all(end+1,:) = min_dist; %#ok
        
        
        apex_prior_scores_all(end+1,:) = sum(sum(apex_prior_scores(:,:))) / num_streams; %#ok
        %smoo_prior_scores_all(end+1,:) = sum(sum(smoo_prior_scores(10:40,:))) / num_streams; %#ok
        rand_prior_scores_all(end+1,:) = sum(sum(rand_prior_scores(:,:))) / num_streams; %#ok

        if i_ve > 20 && i_ve <= 30
            figure; 
            subplot(2,3,1); imgray(apex_path_map);
            plot(v_pts(:,1), v_pts(:,2));
            plot(v_pts(apex_idx,1), v_pts(apex_idx,2), 'rx');
            plot(v_pts(rand_idx,1), v_pts(rand_idx,2), 'go');
            
            subplot(2,3,2); imgray(complex2rgb(ori_patch));
            plot(v_pts(:,1), v_pts(:,2));
            plot(v_pts(apex_idx,1), v_pts(apex_idx,2), 'rx');
            
            subplot(2,3,3); plot(apex_prior_scores);
            
            title(num2str(sum(sum(apex_prior_scores(10:40,:)))));
            
            subplot(2,3,4); imgray(rand_path_map);
            plot(v_pts(:,1), v_pts(:,2));
            plot(v_pts(apex_idx,1), v_pts(apex_idx,2), 'rx');
            plot(v_pts(rand_idx,1), v_pts(rand_idx,2), 'go');
            
            subplot(2,3,5); imgray(prob_patch);
            plot(v_pts(:,1), v_pts(:,2));
            plot(v_pts(apex_idx,1), v_pts(apex_idx,2), 'rx');
            
            subplot(2,3,6); plot(rand_prior_scores);
            title(num2str(sum(sum(rand_prior_scores(10:40,:)))));
            
        end
    end
    %
end
%%
num_streams = 100;
marker_colors = jet(256);
v_count = 0;
curvature_smoothing_sigma = 0.5;
strong_vessel_thresh = 0.25;
stop_count = 0;
stop_at = 100;
for i_sz = {'normal', 'enlarged'}
    
    vessel_size = i_sz{1};
    apex_dir = [nailfoldroot 'data/' rsa_dir '/apexes/' vessel_size '/'];
    contour_dir = [nailfoldroot 'data/' rsa_dir '/vessel_contours/' vessel_size '/'];

    %Get list of vessel names
    v_list = dir([ori_dir, vessel_size '*.mat']);
    %
    for i_ve = 1:length(v_list)%1:20%
    
        v_count = v_count + 1;
        
        if apex_corr_scores_all(v_count) > 0.4 %apex_prior_scores_all(v_count) > 100
            continue;
        end
        stop_count = stop_count + 1;
        if stop_count > stop_at
            break;
        end
        
        %load vessel
        apex_name = v_list(i_ve).name(length(vessel_size)+(1:8));
        contour_struc = load([contour_dir apex_name '_vessel_contour.mat']);
        apex_struc = load([apex_dir apex_name '.mat']);
        vessel_struc = u_load([apex_dir apex_name '_vessel.mat']);
        
        prob_patch = u_load([prob_dir v_list(i_ve).name]);
        prob_patch = conv2(g', g, prob_patch, 'same');
        ori_patch = u_load([ori_dir v_list(i_ve).name]);    

        %Resample the vessel pts to be equalled spaced, 1 pixel apart
        v_pts = spline_contour(contour_struc.vessel_centre, [], 1);

        %Correct apex coordinates frame
        apex_xy(:,1) = apex_struc.apex_xy(:,1) +...
            apex_struc.apex_properties.sc - vessel_struc.vessel_properties.sc;
        apex_xy(:,2) = apex_struc.apex_xy(:,2) +...
            apex_struc.apex_properties.sr - vessel_struc.vessel_properties.sr;

        %Find the vessel pt nearets the apex        
        dists = sum(bsxfun(@minus, v_pts, mean(apex_xy)).^2, 2);
        [~, apex_idx] = min(dists);
        
        %Compute curvature and correlation
        curv_patch = -(complex_gaussian_curvature(ori_patch ./ (abs(ori_patch)+1e-6), 0.5));
        [~, ~, corr_map] = template_match_apexes(vessel_struc.vessel_patch, ...
            {'g2d', apex_template.g2}, 'threshold', -1);
        
        
        
        %Find vessel centres
        vessel_nms = mb_non_maximal_supp(prob_patch, angle(ori_patch)/2);
        strong_vessels = vessel_nms > strong_vessel_thresh;
        if any(strong_vessels(:))
            [rstrong cstrong] = find(strong_vessels);
            vessel_centre = bwselect(vessel_nms > 0, cstrong, rstrong, 8);
        else
            vessel_centre = strong_vessels;
        end
        
        %Update the centres to include only high curvatur epoints
        %vessel_centre = vessel_centre & (curv_patch > 0.05) & (prob_patch > 0.8);
        
        %Now take only the highest curvature from each connected section
%         centre_labels = bwlabel(vessel_centre, 8);
%         for i_lab = 1:max(centre_labels(:))
%             label_mask = centre_labels == i_lab;
%             label_curv = curv_patch(label_mask);
%             label_bw = false(size(label_curv));
%             [~, max_idx] = max(label_curv);
%             label_bw(max_idx) = 1;
%             vessel_centre(label_mask) = label_bw;
%         end
        [vessel_centre_y vessel_centre_x] = find(vessel_centre);
        
        %Scale curv between 1 and 256 for plotting
        vessel_curv = curv_patch(vessel_centre);
        vessel_curv_idx = (vessel_curv - min(vessel_curv)) / (max(vessel_curv)-min(vessel_curv));
        vessel_curv_idx = round(255*vessel_curv_idx)+1;
        
        [apex_prior_scores, apex_path_map] = apex_prior_prob_track_mult(...
                    prob_patch, ori_patch,...
                    vessel_centre_y, vessel_centre_x,...
                    'num_streams', num_streams, 'stopping_prob', 0.1,...
                    'apex_prior', apex_prior);
                
                
        apex_prior_scores_sum = squeeze(sum(sum(apex_prior_scores(:,:,:),2)));
        [max_prior_score, max_prior_idx] = max(apex_prior_scores_sum);
        
        corr_scores = interp2(corr_map, vessel_centre_x, vessel_centre_y, 'nearest');
        [max_corr_score, max_corr_idx] = max(corr_scores);
        
        figure;

        subplot(1,3,1); imgray(complex2rgb(ori_patch));
        plot(v_pts(:,1), v_pts(:,2));
        plot(v_pts(apex_idx,1), v_pts(apex_idx,2), 'rx');       
        plot(vessel_centre_x, vessel_centre_y, 'w.');
        plot(vessel_centre_x(max_prior_idx), vessel_centre_y(max_prior_idx), 'go');
        plot(vessel_centre_x(max_corr_idx), vessel_centre_y(max_corr_idx), 'mo');
        title(num2str(apex_prior_scores_all(v_count)));
        
        subplot(1,3,2); imgray(prob_patch);
        plot(v_pts(:,1), v_pts(:,2));
        plot(v_pts(apex_idx,1), v_pts(apex_idx,2), 'rx');
        
        for i_pt = 1:length(vessel_centre_x)
            pt_color = marker_colors(vessel_curv_idx(i_pt),:);
            plot(vessel_centre_x(i_pt), vessel_centre_y(i_pt), '.', 'markeredgecolor', pt_color);
        end
        plot(vessel_centre_x(max_prior_idx), vessel_centre_y(max_prior_idx), 'go');
        plot(vessel_centre_x(max_corr_idx), vessel_centre_y(max_corr_idx), 'mo');        
        title(num2str(max_prior_score));
        
        subplot(1,3,3); imgray(corr_map);
        plot(vessel_centre_x(max_prior_idx), vessel_centre_y(max_prior_idx), 'go');
        plot(vessel_centre_x(max_corr_idx), vessel_centre_y(max_corr_idx), 'mo');
        title(num2str(max_corr_score));
        
    end
    %
end
%%
v_count = 0;
for i_sz = {'normal', 'enlarged'}
    
    vessel_size = i_sz{1};
    apex_dir = [nailfoldroot 'data/' rsa_dir '/apexes/' vessel_size '/'];
    contour_dir = [nailfoldroot 'data/' rsa_dir '/vessel_contours/' vessel_size '/'];

    %Get list of vessel names
    v_list = dir([ori_dir, vessel_size '*.mat']);
    %
    for i_ve = 1:length(v_list)%1:20%
    
        v_count = v_count + 1;
        
        if rand_prior_scores_all(v_count) < 120 || rand_prior_scores_all(v_count) > 150
            continue;
        end
        
        %load vessel
        apex_name = v_list(i_ve).name(length(vessel_size)+(1:8));
        contour_struc = load([contour_dir apex_name '_vessel_contour.mat']);
        apex_struc = load([apex_dir apex_name '.mat']);
        vessel_struc = u_load([apex_dir apex_name '_vessel.mat']);
        
        prob_patch = u_load([prob_dir v_list(i_ve).name]);
        prob_patch = conv2(g', g, prob_patch, 'same');
        ori_patch = u_load([ori_dir v_list(i_ve).name]);    

        %Resample the vessel pts to be equalled spaced, 1 pixel apart
        v_pts = spline_contour(contour_struc.vessel_centre, [], 1);


        %Correct apex coordinates frame
        apex_xy(:,1) = apex_struc.apex_xy(:,1) +...
            apex_struc.apex_properties.sc - vessel_struc.vessel_properties.sc;
        apex_xy(:,2) = apex_struc.apex_xy(:,2) +...
            apex_struc.apex_properties.sr - vessel_struc.vessel_properties.sr;

        %Find the vessel pt nearets the apex        
        dists = sum(bsxfun(@minus, v_pts, mean(apex_xy)).^2, 2);
        [~, apex_idx] = min(dists);

        rand_idx = rand_idxs(v_count);
        
        figure;

        subplot(1,2,1); imgray(complex2rgb(ori_patch));
        plot(v_pts(:,1), v_pts(:,2));
        plot(v_pts(apex_idx,1), v_pts(apex_idx,2), 'rx');
        plot(v_pts(rand_idx,1), v_pts(rand_idx,2), 'go');

        title(num2str(apex_prior_scores_all(v_count)));

        subplot(1,2,2); imgray(prob_patch);
        plot(v_pts(:,1), v_pts(:,2));
        plot(v_pts(apex_idx,1), v_pts(apex_idx,2), 'rx');
        plot(v_pts(rand_idx,1), v_pts(rand_idx,2), 'go');

        title(num2str(rand_prior_scores_all(v_count)));

    end
    %
end
%%
v_count = 0;
for i_sz = {'normal', 'enlarged'}
    
    vessel_size = i_sz{1};
    apex_dir = [nailfoldroot 'data/' rsa_dir '/apexes/' vessel_size '/'];
    contour_dir = [nailfoldroot 'data/' rsa_dir '/vessel_contours/' vessel_size '/'];

    %Get list of vessel names
    v_list = dir([ori_dir, vessel_size '*.mat']);
    %
    for i_ve = 1:length(v_list)%1:20%
    
        v_count = v_count + 1;
        
        if apex_corr_scores_all(v_count) > 0.4
            continue;
        end
        
        %load vessel
        apex_name = v_list(i_ve).name(length(vessel_size)+(1:8));
        contour_struc = load([contour_dir apex_name '_vessel_contour.mat']);
        apex_struc = load([apex_dir apex_name '.mat']);
        vessel_struc = u_load([apex_dir apex_name '_vessel.mat']);
        
        prob_patch = u_load([prob_dir v_list(i_ve).name]);
        prob_patch = conv2(g', g, prob_patch, 'same');
        ori_patch = u_load([ori_dir v_list(i_ve).name]);    

        %Resample the vessel pts to be equalled spaced, 1 pixel apart
        v_pts = spline_contour(contour_struc.vessel_centre, [], 1);


        %Correct apex coordinates frame
        apex_xy(:,1) = apex_struc.apex_xy(:,1) +...
            apex_struc.apex_properties.sc - vessel_struc.vessel_properties.sc;
        apex_xy(:,2) = apex_struc.apex_xy(:,2) +...
            apex_struc.apex_properties.sr - vessel_struc.vessel_properties.sr;

        %Find the vessel pt nearets the apex        
        dists = sum(bsxfun(@minus, v_pts, mean(apex_xy)).^2, 2);
        [~, apex_idx] = min(dists);

        rand_idx = rand_idxs(v_count);
        
        figure;

        subplot(1,2,1); imgray(complex2rgb(ori_patch));
        plot(v_pts(:,1), v_pts(:,2));
        plot(v_pts(apex_idx,1), v_pts(apex_idx,2), 'rx');
        plot(v_pts(rand_idx,1), v_pts(rand_idx,2), 'go');

        title(num2str(apex_prior_scores_all(v_count)));

        subplot(1,2,2); imgray(prob_patch);
        plot(v_pts(:,1), v_pts(:,2));
        plot(v_pts(apex_idx,1), v_pts(apex_idx,2), 'rx');
        plot(v_pts(rand_idx,1), v_pts(rand_idx,2), 'go');

        title(num2str(rand_prior_scores_all(v_count)));

    end
    %
end
%%
figure; 
imgray(curv_stats.l.map_curv_by_step.hist);
title('Curvature computed from orientation map, left limb');
xlabel('Curvature');
ylabel('Step length from apex');

figure; 
imgray(curv_stats.r.map_curv_by_step.hist);
title('Curvature computed from orientation map, right limb');
xlabel('Curvature');
ylabel('Step length from apex');

figure; 
imgray(curv_stats.l.marked_curv_by_step.hist);
title('Curvature computed from marked vessel contour, left limb');
xlabel('Curvature');
ylabel('Step length from apex');

figure; 
imgray(curv_stats.r.marked_curv_by_step.hist);
title('Curvature computed from marked vessel contour, right limb');
xlabel('Curvature');
ylabel('Step length from apex');

figure; 
imgray(curv_stats.l.pred_curv_by_step.hist);
title('Curvature computed from predicted orientations along marked vessel contour, left limb');
xlabel('Curvature');
ylabel('Step length from apex');

figure; 
imgray(curv_stats.r.pred_curv_by_step.hist);
title('Curvature computed from predicted orientations along marked vessel contour, right limb');
xlabel('Curvature');
ylabel('Step length from apex');

figure; 
imgray(curv_stats.l.marked_ori_by_step.hist);
title('Orientations along marked vessel contour, right limb, left limb');
xlabel('Orientation');
ylabel('Step length from apex');

figure; 
imgray(curv_stats.r.marked_ori_by_step.hist);
title('Orientations along marked vessel contour, right limb');
xlabel('Orientation');
ylabel('Step length from apex');

figure; 
imgray(curv_stats.l.pred_ori_by_step.hist);
title('Predicted orientations along marked vessel contour, right limb, left limb');
xlabel('Orientation');
ylabel('Step length from apex');

figure; 
imgray(curv_stats.r.pred_ori_by_step.hist);
title('Predicted orientations along marked vessel contour, right limb');
xlabel('Orientation');
ylabel('Step length from apex');

figure; imgray(curv_stats.l.pred_error_by_step.hist);
plot([90 90], [0 100], 'r');
title('Error in predicted orientations along marked vessel contour, left limb');
xlabel('Orientation');
ylabel('Step length from apex');

figure; imgray(curv_stats.r.pred_error_by_step.hist);
plot([90 90], [0 100], 'r');
title('Error in predicted orientations along marked vessel contour, right limb');
xlabel('Orientation');
ylabel('Step length from apex');

figure; hold on; 
plot(angle(curv_stats.l.pred_error_by_step.dist)); plot(angle(curv_stats.r.pred_error_by_step.dist), 'r');
title('Mean angle of predicted orientation errors as a function of step length from axis');
ylabel('Orientation error');
xlabel('Step length from apex');
legend({'Left limb', 'Right limb'});
%%

DECOMP_TYPE="{'g2da','h2da'}" WIN_SIZE=3 PROBABILITY_DIR="88126/" MAKE_RESAMPLING_MAPS=0.75 MAX_N_IMAGES=10 SHIFT_IMAGES=10 OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" NUM_TREES=10 IMAGE_TYPE="real" IMAGE_ROOT="retinograms/STARE/training" NUM_SAMPLES=100000 MODEL_ROOT="models/vessel" qsub -V -l highmem -t 1-20 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="{'g2da','h2da'}" PROBABILITY_DIR="222836/" MAKE_RESAMPLING_MAPS=0.75 NUM_ANGLES=6 WIN_SIZE=3 SIGMA_RANGE="[4 8 16]" OUTPUT_TYPE="detection" PREDICTION_TYPE="rf_classification" MASK_DIR="vessel_centre_masks" NUM_TREES=10 IMAGE_TYPE="real" DATA_ROOT="scratch/nailfold/" IMAGE_ROOT="data/set12" NUM_SAMPLES=50000 MODEL_ROOT="models/vessel" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DATA_ROOT="scratch/nailfold/models/" MODEL_ROOT="vessel/detection/rf_classification" MODEL_PATH="'246035'" MAKE_SAMPLED_MAPS=1 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
%%
rsa_dir = 'rsa_study/';

model_dir = [nailfoldroot 'data/' rsa_dir 'models/apex_templates/'];
prob_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/detection/rf_classification/222836/'];
ori_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/orientation/rf_regression/222835/'];
centre_dir = [nailfoldroot 'data/' rsa_dir 'test/vessel_centres/'];
fov_mask_dir = ['C:\isbe\nailfold\data\' rsa_dir 'test\fov_masks\'];

markup_names = dir([nailfoldroot 'data/' rsa_dir 'markup/']);
markup_dirs = cell(0,1);
for i_dir = 3:length(markup_names)
    if markup_names(i_dir).isdir
        markup_dirs{end+1,1} = [nailfoldroot 'data/' rsa_dir 'markup/' markup_names(i_dir).name '/'];
    end
end

%Create smoothing kernel for vessel probs
vessel_prob_smoothing_sigma = 2;
g = gaussian_filters_1d(vessel_prob_smoothing_sigma);
g = g / sum(g);

load([model_dir 'apex_curv_stats.mat']);
ori_hist_r_smoothed = imfilter(curv_stats.r.pred_ori_by_step.hist, conv(g,g), 'circular');
ori_hist_r_smoothed = imfilter(ori_hist_r_smoothed, conv(g,g)', 'replicate');

ori_hist_l_smoothed = imfilter(curv_stats.l.pred_ori_by_step.hist, conv(g,g), 'circular');
ori_hist_l_smoothed = imfilter(ori_hist_l_smoothed, conv(g,g)', 'replicate');

apex_prior = cat(3, ori_hist_l_smoothed, ori_hist_r_smoothed);
%
test_list = dir([centre_dir '*_vc.mat']);

apex_corr_sections = zeros(length(test_list),1);
distal_count = zeros(length(test_list),1);
remaining_pts = zeros(length(test_list),1);

num_y_bins = 100;
y_bins = linspace(-1000, 1000, num_y_bins);
distal_y_hist = zeros(length(test_list),num_y_bins);

num_streams = 100;
%
for i_im = 22:30%6:10%length(test_list)
    
    im_num = test_list(i_im).name(1:6);
    
    load([centre_dir im_num '_vc.mat']);
    
    %First up, discard points from the side of the image or any long
    %horizontal lines
    %vessel_prob = u_load([prob_dir im_num '_pred.mat']);
    
    %[nrows ncols] = size(vessel_prob);
    vessel_centre_mask = false(nrows, ncols);
    vessel_centre_idx = sub2ind([nrows ncols], vessel_centre_y, vessel_centre_x);
    vessel_centre_mask(vessel_centre_idx) = 1;
    %
    %figure; imgray(vessel_prob);
    clear vessel_prob;
    
    centre_labels = bwlabel(vessel_centre_mask, 8);
    centre_labels = centre_labels(vessel_centre_mask);
    
    discard_points = false(length(vessel_centre_x),1);
    apex_corr_max_points = false(length(vessel_centre_x),1);
    
    
    for i_lab = 1:max(centre_labels(:))
        label_mask = centre_labels == i_lab;     
        label_corr = vessel_centre_corr(label_mask);
            
        %plot(label_x, label_y, 'b.');

        if any(label_corr > 0.4)
            discard_points(label_mask) = 1;
            %plot(label_x, label_y, 'm.');
            apex_corr_sections(i_im) = apex_corr_sections(i_im) + 1;
            label_apex_corr_max_points = false(size(label_corr));
            [~,max_idx] = max(label_corr);
            label_apex_corr_max_points(max_idx) = 1;
            apex_corr_max_points(label_mask) = label_apex_corr_max_points;         
        end
        
    end
    max_corr_apex_x = vessel_centre_x(apex_corr_max_points);
    max_corr_apex_y = vessel_centre_y(apex_corr_max_points);
    
    %[counts bins] = hist(max_apex_y, 1:50:nrows);
    %[~,mi] = max(counts);
    %mid_y = bins(mi);
    mid_y = median(max_corr_apex_y);
    
    %Now select only probable, highly curved points    
    discard_points = discard_points | (vessel_centre_y < mid_y-250) | ...
        (vessel_centre_y > mid_y+200) | (vessel_centre_prob < 0.4) | ...
        (abs(vessel_centre_curv) < 0.01);

    remaining_pts(i_im) = sum(~discard_points);
    display(['Original num pts = ' num2str(length(vessel_centre_x)) ' remaining pts = ' num2str(remaining_pts(i_im))]);
    
    %Discard these points
    apex_mask = vessel_centre_mask;   
    apex_mask(vessel_centre_mask) = ~discard_points;
    apex_x = vessel_centre_x(~discard_points);
    apex_y = vessel_centre_y(~discard_points);
    clear vessel_centre_*;
         
    vessel_prob = u_load([prob_dir im_num '_pred.mat']);
    vessel_ori = u_load([ori_dir im_num '_pred.mat']);
    vessel_prob = conv2(g', g, vessel_prob, 'same');
    
    tic;
    [apex_prior_scores] = apex_prior_prob_track_mult(...
                    vessel_prob, vessel_ori,...
                    apex_y, apex_x,...
                    'num_streams', num_streams, 'stopping_prob', 0.1,...
                    'apex_prior', apex_prior);
    toc;
    
    %
    apex_prior_scores_sum = squeeze(sum(sum(apex_prior_scores(:,:,:),2)));
    potential_apices = apex_prior_scores_sum > 120;
    %
    %Now take only the highest curvature from each connected section
    centre_labels = bwlabel(apex_mask, 8);
    centre_labels = centre_labels(apex_mask);
    apex_points = false(length(apex_x),1);
    
    for i_lab = 1:max(centre_labels(:))
        label_mask = centre_labels == i_lab;
        label_prior_sum = apex_prior_scores_sum(label_mask);
        label_bw = false(sum(label_mask),1);
        [~, max_idx] = max(label_prior_sum);
        label_bw(max_idx) = 1;
        apex_points(label_mask) = label_bw;
    end
    
    max_prior_apex_x = apex_x(apex_points & potential_apices);
    max_prior_apex_y = apex_y(apex_points & potential_apices);
%  
    figure; imgray(vessel_prob);
    plot(apex_x, apex_y, 'b.')
    plot(apex_x(potential_apices), apex_y(potential_apices), 'y.')
    plot(max_prior_apex_x, max_prior_apex_y, 'ms')
    plot(max_corr_apex_x, max_corr_apex_y, 'mo')
    
% 
    %
    for i_dir = 1:length(markup_dirs)
        vessel_markup_list_i = dir([markup_dirs{i_dir} '*' im_num '*.txt']);
        if ~isempty(vessel_markup_list_i)
            markup_dir = markup_dirs{i_dir};
        	vessel_markup_list = vessel_markup_list_i;
        end
    end
    if isempty(vessel_markup_list)
        distal_count(i_im) = -1;
    else
        
        %Read in markup
        vessel_markup = read_markup_from([markup_dir vessel_markup_list(end).name]);
        
        %Loop through each marked vessel apex
        num_vessels = length(vessel_markup.vessels);
        for i_v = 1:num_vessels

            %Check this is a valid vessel
            anchor_xy = vessel_markup.vessels(i_v).anchor;

            if isempty(anchor_xy); continue; end

            %Check if distal
            is_distal = vessel_markup.vessels(i_v).ncm_vessel_properties.is_distal;

            if is_distal
                
                num_apices = length(vessel_markup.vessels(i_v).apices);
                for i_a = 1:num_apices
                    if isempty(vessel_markup.vessels(i_v).apices(i_a).inner_point)                    
                        %Use the anchor
                        plot(anchor_xy(:,1), anchor_xy(:,2), 'ro');

                    else
                        %Compute the centre of the apex
                        apex_xy =  ...
                            [ vessel_markup.vessels(i_v).apices(i_a).outer_point;...
                              vessel_markup.vessels(i_v).apices(i_a).inner_point];
                         
                        apex_width = sqrt(sum(diff(apex_xy).^2));
                        
                        distal_count(i_im) = distal_count(i_im) + 1;
                        distal_y_hist(i_im,:) = distal_y_hist(i_im,:) +...
                            hist(mean(apex_xy(:,2)) - mid_y, y_bins);
                          
                        if apex_width < 50
                            plot(apex_xy(:,1), apex_xy(:,2), 'g-*')
                        else
                            plot(apex_xy(:,1), apex_xy(:,2), 'r-*')
                        end
                        

                    end
                end
            else
                %Mask out the region around the anchor
                plot(anchor_xy(:,1), anchor_xy(:,2), 'r^');
                
            end
        end
    end
    %
    display(['Num marked distal apexes = ' num2str(distal_count(i_im)) ', num apex xcorr sections = ' num2str(apex_corr_sections(i_im))]);
    
end
