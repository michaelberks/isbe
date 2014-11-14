load('C:\isbe\nailfold\data\rsa_study\apexes\giant\apex0314_vessel.mat');
%%
image_dir = 'C:\isbe\nailfold\data\rsa_study\training\images\';
pred_dir = 'C:\isbe\nailfold\data\rsa_study\training\predictions\detection\rf_classification\182262\';

i_list = dir([image_dir '*.mat']);
p_list = dir([pred_dir '*.mat']);
%%
for i_ve = 1:20
    vessel_im = u_load([image_dir i_list(i_ve).name]);
    vessel_pred = u_load([pred_dir p_list(i_ve).name]);
    figure; 
    subplot(1,2,1); imgray(vessel_im);
    subplot(1,2,2); imgray(vessel_pred);
end
%%
image_dir = 'C:\isbe\nailfold\data\rsa_study\training\images\';
pred_dir = 'C:\isbe\nailfold\data\rsa_study\training\predictions\detection\rf_classification\182321\';
ori_dir = 'C:\isbe\nailfold\data\rsa_study\training\predictions\orientation\rf_regression\182263\';
width_dir = 'C:\isbe\nailfold\data\rsa_study\training\predictions\width\rf_regression\182367\';


i_list = dir([image_dir '*.mat']);
p_list = dir([pred_dir '*.mat']);
o_list = dir([ori_dir '*.mat']);
w_list = dir([width_dir '*.mat']);
g = gaussian_filters_1d(2);
g = g / sum(g);
%%
%clims = zeros(20,2);
%clim = mean(clims);
for i_ve = 1:20
    vessel_im = u_load([image_dir i_list(i_ve).name]);
    vessel_pred = u_load([pred_dir p_list(i_ve).name]);
    vessel_ori = u_load([ori_dir o_list(i_ve).name]);
    vessel_width = u_load([width_dir w_list(i_ve).name]);
    

    vessel_ori = vessel_ori ./ (abs(vessel_ori) + 1e-6);
    [potential_map scale_map, rotation_map, mean_width_map, mean_ori_map, mean_curv_map] = ...
        estimate_vessel_pose3(vessel_pred, vessel_width, vessel_ori, 'patch_size', 33); 
    ori_diff_map = mean(diff(abs(mean_ori_map),1,3),3);
    
    [vy vx] = find(potential_map);
    %clims(i_ve,:) = [min(mean_curv_map(potential_map)) max(mean_curv_map(potential_map))];
    figure; 
    subplot(1,2,1); imgray(two_map_colour(potential_map, mean_curv_map(:,:,1)));
    subplot(1,2,2); imgray(two_map_colour(potential_map, ori_diff_map, [], [-0.05 0.27]));

    
end
%%
for i_ve = 81:102
    vessel_im = u_load([image_dir i_list(i_ve).name]);
    vessel_pred = u_load([pred_dir p_list(i_ve).name]);
    vessel_ori = u_load([ori_dir o_list(i_ve).name]);
    
    valid_mask = ~isnan(vessel_pred);
    vessel_pred(~valid_mask) = 0;
    valid_weights = conv2(g', g, double(valid_mask), 'same');
    vessel_pred_smoothed = conv2(g', g, vessel_pred, 'same') ./ valid_weights;

    vessel_nms = mb_non_maximal_supp(vessel_pred_smoothed, angle(vessel_ori)/2);
    
    strong_vessels = vessel_nms > 0.4;
    if any(strong_vessels(:))
        [rstrong cstrong] = find(strong_vessels);
        combined_vessels = bwselect(vessel_nms > 0, cstrong, rstrong, 8);
    else
        combined_vessels = strong_vessels;
    end
    
    figure; 
    subplot(1,2,1); imgray(complex2rgb(vessel_pred_smoothed.*vessel_ori));
    subplot(1,2,2); imgray(combined_vessels);
end
%%
for i_ve = 1:20
    vessel_im = u_load([image_dir i_list(i_ve).name]);
    vessel_pred = u_load([pred_dir p_list(i_ve).name]);
    vessel_width = u_load([width_dir w_list(i_ve).name]);
    figure; 
    subplot(1,2,1); imgray(vessel_im);
    subplot(1,2,2); imgray(two_map_colour(vessel_pred, vessel_width, [], [10 40]));
end
%%
for i_ve = 21:40
    vessel_im = u_load([image_dir i_list(i_ve).name]);
    vessel_pred = u_load([pred_dir p_list(i_ve).name]);
    vessel_ori = u_load([ori_dir o_list(i_ve).name]);
    figure; 
    subplot(1,2,1); imgray(vessel_im);
    subplot(1,2,2); imgray(complex2rgb(vessel_pred.*vessel_ori));
end
%%
%Check how many vessel and centre points we have in the training data
vessel_mask_dir = 'C:\isbe\nailfold\data\rsa_study\training\vessel_masks\';
vessel_cmask_dir = 'C:\isbe\nailfold\data\rsa_study\training\vessel_centre_masks\';

cmask_list = dir([vessel_cmask_dir '*.mat']);
vmask_list = dir([vessel_mask_dir '*.mat']);

for i_ve = 1:20
    c_mask = u_load([vessel_cmask_dir cmask_list(i_ve).name]);
    v_mask = u_load([vessel_mask_dir vmask_list(i_ve).name]);
    c_mask2 = bwmorph(v_mask, 'thin', 'inf');
    figure;
    subplot(1,3,1); imgray(v_mask);
    subplot(1,3,2); imgray(c_mask);
    subplot(1,3,3); imgray(c_mask2);
end
%%
%--------------------------------------------------------------------------
%Make lists of file locations used in the training data
data_lists_dir = 'C:\isbe\nailfold\data\rsa_study\data_lists\';
prob_dir = 'C:\isbe\nailfold\data\rsa_study\training\predictions\detection\rf_classification\182321\';
ori_dir = 'C:\isbe\nailfold\data\rsa_study\training\predictions\orientation\rf_regression\182263\';
width_dir = 'C:\isbe\nailfold\data\rsa_study\training\predictions\width\rf_regression\182367\';
create_folder(data_lists_dir);

v_list = dir('C:\isbe\nailfold\data\rsa_study\training\images\*.mat');
num_images = length(v_list);

vessel_files = cell(num_images,1);
contour_files = cell(num_images,1);
apex_files = cell(num_images,1);
vessel_prob_files = cell(num_images,1);
vessel_ori_files = cell(num_images,1);
vessel_width_files = cell(num_images,1);

counter = 1;
vessel_sizes = {'normal', 'enlarged'};
for i_sz = 1:2
    vessel_dir = ['C:\isbe\nailfold\data\rsa_study\apexes\' vessel_sizes{i_sz} '\'];
    contour_dir = ['C:\isbe\nailfold\data\rsa_study\vessel_contours\' vessel_sizes{i_sz} '\'];
    v_files = dir([contour_dir '*contour.mat']);
    
    num_vessels = length(v_files);
       
    for i_v = 1:length(v_files)
        %load data
        vessel_files{counter} = [vessel_dir v_files(i_v).name(1:8) '_vessel.mat'];
        contour_files{counter} = [contour_dir v_files(i_v).name];
        apex_files{counter} = [vessel_dir v_files(i_v).name(1:8) '.mat'];    
        vessel_prob_files{counter} = [prob_dir vessel_sizes{i_sz} v_files(i_v).name(1:8) '_vessel_pred.mat'];
        vessel_ori_files{counter} = [ori_dir vessel_sizes{i_sz} v_files(i_v).name(1:8) '_vessel_pred.mat'];
        vessel_width_files{counter} = [width_dir vessel_sizes{i_sz} v_files(i_v).name(1:8) '_vessel_pred.mat'];
        counter = counter + 1;
    end
end
save([data_lists_dir 'training_data_lists.mat'], '*_files');
%%
clear;
data_lists_dir = 'C:\isbe\nailfold\data\rsa_study\data_lists\';
load([data_lists_dir 'training_data_lists.mat'], '*_files');
[unaligned_shape_data aligned_shape_data] = vessel_detection_shape_alignment(... 
    'contour_files', contour_files,...
    'vessel_files', vessel_files,...
    'model_dir', 'C:\isbe\nailfold\data\rsa_study\models\apex_templates\',...
    'num_shape_pts', 31,...
    'plot', 1);
%%
apex_template = vessel_detection_template_building( ...
    'apex_shapes', unaligned_shape_data.apex_shapes, ...
    'apex_widths', unaligned_shape_data.apex_widths, ...
    'a_scales', aligned_shape_data.a_scales,...
    'a_rots', aligned_shape_data.a_rots,...
    'model_dir', 'C:\isbe\nailfold\data\rsa_study\models\apex_templates\',...
    'vessel_files', vessel_files, ...
    'vessel_prob_files', vessel_prob_files,...
    'template_size', 49,...
    'do_intensity', 1,...
    'do_g1', 1,...
    'do_g2', 1,...
    'do_vessel_prob', 1,...
    'plot', 0);
figure;
subplot(2,2,1); imgray(mean_apex_aligned_i);
subplot(2,2,2); imgray(mean_apex_aligned_g1);
subplot(2,2,3); imgray(mean_apex_aligned_g2);
subplot(2,2,4); imgray(mean_apex_aligned_p);

%%

model_dir = 'C:\isbe\nailfold\data\rsa_study\models\apex_templates\';
image_dir = 'C:\isbe\nailfold\data\rsa_study\images\';
mask_dir = 'C:\isbe\nailfold\data\rsa_study\test\fov_masks\';
pred_dir = 'C:\isbe\nailfold\data\rsa_study\test\predictions\detection\rf_classification\182321_corrected\';
ori_dir = 'C:\isbe\nailfold\data\rsa_study\test\predictions\orientation\rf_regression\182263\';
width_dir = 'C:\isbe\nailfold\data\rsa_study\test\predictions\width\rf_regression\182367\';
markup_dir = 'C:\isbe\nailfold\data\rsa_study\markup\tmoore\';

load([model_dir 'estimated_width_transforms.mat'], 'width_m', 'width_c');
load('C:\isbe\nailfold\data\rsa_study\models\apex_templates\apex_templates.mat');

i_list = dir([image_dir '*.mat']);
p_list = dir([pred_dir '*.mat']);
o_list = dir([ori_dir '*.mat']);
w_list = dir([width_dir '*.mat']);

nearest_nms = [];
nearest_nms_vals = [];
nearest_c_vals = [];

max_r = 800;
max_c = 800;
%
for i_im = 1:10%length(o_list)
    
    im_name = o_list(i_im).name(1:end-9);
    
    vessel_markup_list = dir([markup_dir '*' im_name '*.txt']);
    
    if isempty(vessel_markup_list) ||...
        ~exist([pred_dir o_list(i_im).name], 'file') ||...
        ~exist([width_dir o_list(i_im).name], 'file');
        continue;
    end 
    
    vessel_markup = read_markup_from([markup_dir vessel_markup_list(end).name]);
    vessel_pred = u_load([pred_dir o_list(i_im).name]);
    
    [rows cols] = size(vessel_pred);
    if rows <= max_r
        keep_rows = 1:rows;
        off_y = 0;
    else
        off_y = floor((rows-max_r)/2);
        keep_rows = off_y+(1:max_r);
    end
    if cols <= max_c
        keep_cols = 1:cols;
        off_x = 0;
    else
        off_x = floor((cols-max_c)/2);
        keep_cols = off_x+(1:max_c);
    end
    vessel_pred = vessel_pred(keep_rows, keep_cols);
    
    vessel_ori = u_load([ori_dir o_list(i_im).name]);
    vessel_ori = vessel_ori(keep_rows, keep_cols);
    vessel_width = u_load([width_dir o_list(i_im).name]);
    vessel_width = vessel_width(keep_rows, keep_cols);
    
    vessel_im = imread([image_dir im_name '.png']);
    vessel_im = vessel_im(keep_rows, keep_cols);
    vessel_im_eq = equalise_nailfold_intensity(vessel_im);
    
%     figure; imgray(vessel_pred);
    
    [potential_map scale_map, rotation_map] = ...
        estimate_vessel_pose(vessel_pred, vessel_width, vessel_ori,...
        'width_m', width_m, 'width_c', width_c, 'patch_size', 101);
    
    potential_map = imdilate(potential_map, strel('disk', 1));
    
    vessel_pred_nms = vessel_pred;
    vessel_pred_nms(~potential_map) = 0;
    
    gp = prctile(double(vessel_im(:)), 0:100);
    [~, ip] = max(diff(gp(51:end)));
    g_min = gp(2);
    g_max = gp(ip+50);
    
    figure; 
    a1 = subplot(1,2,1); imgray(vessel_im); caxis([g_min g_max]);
    a2 = subplot(1,2,2); imgray(vessel_pred_nms);
    linkaxes([a1 a2]);
    
%     C2 = targetted_template_matching(...
%         vessel_im_eq, potential_map, scale_map, rotation_map, apex_template.intensity);      
%     
%     if i_im <= 10
%         %g_lims = prctile(double(vessel_im(:)), [2 98]);
%         figure; 
%         a1 = subplot(2,1,1); imgray(ind2rgb(vessel_im, gray(256)));
%         a2 = subplot(2,1,2); imgray(C2); colormap(jet(256)); caxis([0 1]);
%         linkaxes([a1 a2]);
%     end
    
    [poss_vessel_pts_y poss_vessel_pts_x] = find(potential_map);
    num_vessels = length(vessel_markup.vessels);
        
    for i_v = 1:num_vessels

        %Check this is a valid vessel
        anchor_xy = vessel_markup.vessels(i_v).anchor;

        if isempty(anchor_xy); continue; end          
        
        if mean(anchor_xy(:,1))>(max_c+off_x) || mean(anchor_xy(:,2))>(max_r+off_y)
            continue;
        end

        %Check if distal
        is_distal = vessel_markup.vessels(i_v).ncm_vessel_properties.is_distal;

        if is_distal
            num_apices = length(vessel_markup.vessels(i_v).apices);
            for i_a = 1:num_apices
                if isempty(vessel_markup.vessels(i_v).apices.inner_point)
                    %plot the anchor
                    plot(a1, -off_x+anchor_xy(1), -off_y+anchor_xy(2), 'go');
                    plot(a2, -off_x+anchor_xy(1), -off_y+anchor_xy(2), 'go');
                else
                    %plot the apex
                    apex_xy = ...
                        [ vessel_markup.vessels(i_v).apices(i_a).outer_point;...
                          vessel_markup.vessels(i_v).apices(i_a).inner_point];
                    plot(a1, -off_x+apex_xy(:,1), -off_y+apex_xy(:,2), 'r-x');
                    plot(a2, -off_x+apex_xy(:,1), -off_y+apex_xy(:,2), 'r-x');
                    
                    ax = mean(apex_xy(:,1));
                    ay = mean(apex_xy(:,2));
                    dists = (poss_vessel_pts_x-ax).^2 + (poss_vessel_pts_y-ay).^2;
                    [min_dist, min_idx] = min(dists);
                    nearest_nms(end+1,1) = sqrt(min_dist); %#ok
                    %nearest_nms_vals(end+1,1) =...
                    %    vessel_nms(poss_vessel_pts_y(min_idx), poss_vessel_pts_x(min_idx)); %#ok
                end
            end
        else
            %plot the anchor
            plot(a1, -off_x+anchor_xy(1), -off_y+anchor_xy(2), 'g+');
            plot(a2, -off_x+anchor_xy(1), -off_y+anchor_xy(2), 'g+');
        end
    end
end
%%
load([model_dir 'apex_shape_data_prealignment.mat'], 'apex_widths', 'apex_shapes', 'apex_names', 'apex_theta');
load([model_dir 'apex_shape_data.mat'], 'a_scales', 'a_rots', 'apex_theta', 'a_shapes', 'a_trans');
figure; plot(1./a_scales, apex_widths, 'rx');
axis([0 0.4 0 70]);
hold on;
plot([0 0.4], [0 70], 'k--');

image_dir = 'C:\isbe\nailfold\data\rsa_study\training\images\';
pred_dir = 'C:\isbe\nailfold\data\rsa_study\training\predictions\detection\rf_classification\182321\';
ori_dir = 'C:\isbe\nailfold\data\rsa_study\training\predictions\orientation\rf_regression\182263\';
width_dir = 'C:\isbe\nailfold\data\rsa_study\training\predictions\width\rf_regression\182367\';

g = gaussian_filters_1d(2);
g = g / sum(g);
%
measured_widths = zeros(102,1);
measured_oris = zeros(102,1);
for i_ve = 1:102

    vessel_pred = u_load([pred_dir apex_names{i_ve} '_vessel_pred.mat']);
    vessel_width = u_load([width_dir apex_names{i_ve} '_vessel_pred.mat']);
    vessel_ori = u_load([ori_dir apex_names{i_ve} '_vessel_pred.mat']);
    
%     valid_mask = ~isnan(vessel_pred);
%     vessel_pred(~valid_mask) = 0;
%     valid_weights = conv2(g', g, double(valid_mask), 'same');
%     vessel_pred_smoothed = conv2(g', g, vessel_pred, 'same') ./ valid_weights;
    
    valid_mask = ~isnan(vessel_width) & ~isnan(vessel_pred) & ~isnan(vessel_ori);
    vessel_width(~valid_mask) = 0;
    vessel_pred(~valid_mask) = 0;
    vessel_ori(~valid_mask) = 0;
    
    pred_weights = conv2(g', g, double(valid_mask), 'same');
    vessel_pred_smoothed = conv2(g', g, vessel_pred, 'same') ./ pred_weights;
    
    vessel_weights = conv2(g', g, vessel_pred, 'same');
    
    vessel_width_smoothed = conv2(g', g, vessel_pred.*vessel_width, 'same') ./ vessel_weights;
    vessel_ori_smoothed = conv2(g', g, vessel_pred.*vessel_ori, 'same') ./ vessel_weights;
            
    apex_xy = reshape(apex_shapes(i_ve, :), 31, 2);
    apex_idx = sub2ind(size(vessel_pred), round(apex_xy(:,2)), round(apex_xy(:,1)));
    
    if i_ve <= 10
        figure; 
        subplot(1,2,1); imgray(two_map_colour(vessel_pred_smoothed, vessel_width_smoothed, [], [10 40]));
        plot(apex_xy(:,1), apex_xy(:,2), 'yx', 'markersize', 10);
        plot(apex_xy(16,1), apex_xy(16,2), 'bx', 'markersize', 10);
        subplot(1,2,2); imgray(complex2rgb(vessel_pred_smoothed .* vessel_ori_smoothed));
        plot(apex_xy(:,1), apex_xy(:,2), 'yx', 'markersize', 10);
        plot(apex_xy(16,1), apex_xy(16,2), 'bx', 'markersize', 10);
    end
    measured_widths(i_ve) = vessel_width_smoothed(round(apex_xy(16,2)), round(apex_xy(16,1)));
    measured_oris(i_ve) = mod(angle(mean(vessel_ori_smoothed(apex_idx)))/2,pi);
end

figure; plot(apex_widths, measured_widths, 'rx');
axis([0 70 0 45]);
hold on;
plot([0 70], [0 45], 'k--');

figure; plot(measured_widths, 1./a_scales, 'rx');
axis([0 max(measured_widths) 0 max(1./a_scales)]);
hold on;
plot([0 max(measured_widths)], [0 max(1./a_scales)], 'k--');
[~, width_m, width_c] = regression(measured_widths', 1./a_scales');
plot([0 max(measured_widths)], width_m*[0 max(measured_widths)]+width_c, 'g');

a_thetas = pi/2 + atan2(squeeze(a_rots(1,2,:)), squeeze(a_rots(1,1,:)));
figure; plot(measured_oris, a_thetas, 'rx');
axis equal; axis([0 pi 0 pi]); hold on;
plot([0 pi], [0 pi], 'k--');

%%
model_dir = 'C:\isbe\nailfold\data\rsa_study\models\apex_templates\';
image_dir = 'C:\isbe\nailfold\data\rsa_study\training\images\';
pred_dir = 'C:\isbe\nailfold\data\rsa_study\training\predictions\detection\rf_classification\182321\';
ori_dir = 'C:\isbe\nailfold\data\rsa_study\training\predictions\orientation\rf_regression\182263\';
width_dir = 'C:\isbe\nailfold\data\rsa_study\training\predictions\width\rf_regression\182367\';

load([model_dir 'apex_templates.mat']);

i_list = dir([image_dir '*.mat']);
p_list = dir([pred_dir '*.mat']);
o_list = dir([ori_dir '*.mat']);
w_list = dir([width_dir '*.mat']);
g = gaussian_filters_1d(2);
g = g / sum(g);

i_ve = 1;
%%
vessel_im = u_load([image_dir i_list(i_ve).name]);
vessel_pred = u_load([pred_dir p_list(i_ve).name]);
vessel_ori = u_load([ori_dir o_list(i_ve).name]);
vessel_width = u_load([width_dir w_list(i_ve).name]);

valid_mask = ~isnan(vessel_pred);
vessel_pred(~valid_mask) = 0;
valid_weights = conv2(g', g, double(valid_mask), 'same');
vessel_pred_smoothed = conv2(g', g, vessel_pred, 'same') ./ valid_weights;

vessel_nms = mb_non_maximal_supp(vessel_pred_smoothed, angle(vessel_ori)/2);
vessel_width(isnan(vessel_width)) = mean(vessel_width(~isnan(vessel_width) & (vessel_pred > 0.25)));

strong_vessels = vessel_nms > 0.4;
if any(strong_vessels(:))
    [rstrong cstrong] = find(strong_vessels);
    combined_vessels = bwselect(vessel_nms > 0, cstrong, rstrong, 8);
else
    combined_vessels = strong_vessels;
end
tic;
C1 = mb_normxcorr2(apex_template.intensity, vessel_im);
toc;
tic;
C2 = targetted_template_matching(...
    vessel_im, combined_vessels, ones(size(vessel_im)), zeros(size(vessel_im)), apex_template.intensity);
toc;
%%
clear;
data_lists_dir = 'C:\isbe\nailfold\data\rsa_study\data_lists\';
load('C:\isbe\nailfold\data\rsa_study\models\apex_templates\unaligned_apex_shape_data.mat');
load('C:\isbe\nailfold\data\rsa_study\models\apex_templates\aligned_apex_shape_data.mat');
load([data_lists_dir 'training_data_lists.mat'], '*_files');

measured_widths = zeros(102,1);
measured_oris = zeros(102,1);
measured_ori_diffs = zeros(102,1);

mean_scale = mean(aligned_shape_data.a_scales);
for i_ve = 1:102

    vessel_pred = u_load(vessel_prob_files{i_ve});
    vessel_ori = u_load(vessel_ori_files{i_ve});
    vessel_width = u_load(vessel_width_files{i_ve});
    
    vessel_pred(isnan(vessel_pred)) = 1;
    vessel_width(isnan(vessel_width)) = mean(vessel_width(~isnan(vessel_width) & (vessel_pred > 0.25)));
    vessel_ori(isnan(vessel_ori)) = mean(vessel_ori(~isnan(vessel_ori) & (vessel_pred > 0.25)));
    [potential_map, ~, rotation_map, mean_width_map, mean_ori_map] = ...
        estimate_vessel_pose3(vessel_pred, vessel_width, vessel_ori, 'patch_size', 33);
            
    apex_xy = unaligned_shape_data.apex_shapes(i_ve, [16 47]);
    [vy vx] = find(potential_map);
    dists = (vx-apex_xy(1)).^2 + (vy-apex_xy(2)).^2;
    [~, min_i] = min(dists);
    v_apex_x = vx(min_i);
    v_apex_y = vy(min_i);
    
    if ismember(i_ve, [34 68 69 70 71 87 96])
        figure; 
        subplot(1,2,1); imgray(complex2rgb(potential_map.*mean_ori_map(:,:,1)));
        plot(apex_xy(1), apex_xy(2), 'gx');
        plot(v_apex_x, v_apex_y, 'go');
        
        subplot(1,2,2); imgray(two_map_colour(potential_map, mean_width_map(:,:,1), []));
        plot(apex_xy(1), apex_xy(2), 'gx');
        plot(v_apex_x, v_apex_y, 'go');
    end
    
    ori_diff_map = mean(diff(abs(mean_ori_map),1,3),3);
    measured_ori_diffs(i_ve) = ori_diff_map(v_apex_y, v_apex_x);
    measured_widths(i_ve) = mean_width_map(v_apex_y, v_apex_x);
    measured_oris(i_ve) = rotation_map(v_apex_y, v_apex_x);
end

figure; plot(unaligned_shape_data.apex_widths, measured_widths, 'rx');
axis([0 70 0 45]);
hold on;
plot([0 70], [0 45], 'k--');

i_scales = mean_scale ./ aligned_shape_data.a_scales;
figure; plot(measured_widths, i_scales, 'rx');
axis([0 max(measured_widths) 0 max(i_scales)]);
hold on;
plot([0 max(measured_widths)], [0 max(i_scales)], 'k--');
[reg_score_widths, width_m, width_c] = regression(measured_widths', i_scales');
plot([0 max(measured_widths)], width_m*[0 max(measured_widths)]+width_c, 'g');

a_thetas = atan2(squeeze(aligned_shape_data.a_rots(1,2,:)), squeeze(aligned_shape_data.a_rots(1,1,:)));
figure; plot(measured_oris, a_thetas, 'rx');
axis equal; axis([-pi/2 pi/2 -pi/2 pi/2]); hold on;
plot([-pi/2 pi/2], [-pi/2 pi/2], 'k--');
[reg_score_oris, ori_m, ori_c] = regression(measured_oris', a_thetas');
%%
data_lists_dir = 'C:\isbe\nailfold\data\rsa_study\data_lists\';
load('C:\isbe\nailfold\data\rsa_study\models\apex_templates\unaligned_apex_shape_data.mat');
load('C:\isbe\nailfold\data\rsa_study\models\apex_templates\aligned_apex_shape_data.mat');
load([data_lists_dir 'training_data_lists.mat'], '*_files');

%Make initial sampling points
template_size = 49;
lim = (template_size-1)/2;
x = repmat(-lim:lim, template_size, 1);
y = x';
xy = [x(:) y(:)];
mean_scale = mean(aligned_shape_data.a_scales);
%%
for i_ve = [34 68 69 70 71 87 96]

    vessel_struc = u_load(vessel_files{i_ve});
    vessel_pred = u_load(vessel_prob_files{i_ve});
    vessel_ori = u_load(vessel_ori_files{i_ve});
    vessel_width = u_load(vessel_width_files{i_ve});
    
    vessel_pred(isnan(vessel_pred)) = 1;
    vessel_width(isnan(vessel_width)) = mean(vessel_width(~isnan(vessel_width) & (vessel_pred > 0.25)));
    vessel_ori(isnan(vessel_ori)) = mean(vessel_ori(~isnan(vessel_ori) & (vessel_pred > 0.25)));
    [potential_map scale_map, rotation_map] = ...
        estimate_vessel_pose(vessel_pred, vessel_width, vessel_ori,...
        'width_m', width_m, 'width_c', width_c, 'patch_size', 101);
            
    apex_xy = unaligned_shape_data.apex_shapes(i_ve, [16 47]);
    [vy vx] = find(potential_map);
    dists = (vx-apex_xy(1)).^2 + (vy-apex_xy(2)).^2;
    [~, min_i] = min(dists);
    v_apex_x = vx(min_i);
    v_apex_y = vy(min_i);
    
    estimated_scale = scale_map(v_apex_y, v_apex_x);
    estimated_theta = rotation_map(v_apex_y, v_apex_x);
    estimated_rot = rotation_matrix(estimated_theta);
    
    xya = xy * estimated_rot * estimated_scale;
    xa = reshape(xya(:,1) + v_apex_x, template_size, template_size);
    ya = reshape(xya(:,2) + v_apex_y, template_size, template_size);
    
    vessel_patch = equalise_nailfold_intensity(vessel_struc.vessel_patch);
    intensity_patch = interp2(vessel_patch, xa, ya, '*bilinear');
    
    if 1%i_ve <= 10
        figure; 
        subplot(1,2,1); imgray(vessel_patch);
        plot(xa, ya, 'b.', 'markersize', 4);
        subplot(1,2,2); imgray(intensity_patch);
    end
    
    
    %measured_dists(i_ve) = mean_dist_map(v_apex_y, v_apex_x);
end
%%
load('C:\isbe\nailfold\data\rsa_study\models\apex_templates\apex_templates.mat');
load('C:\isbe\nailfold\data\rsa_study\models\apex_templates\estimated_width_transforms.mat');
for i_ve = 1:20%81:102

    vessel_struc = u_load(vessel_files{i_ve});
    vessel_pred = u_load(vessel_prob_files{i_ve});
    vessel_ori = u_load(vessel_ori_files{i_ve});
    vessel_width = u_load(vessel_width_files{i_ve});
    
    vessel_pred(isnan(vessel_pred)) = 1;
    vessel_width(isnan(vessel_width)) = mean(vessel_width(~isnan(vessel_width) & (vessel_pred > 0.25)));
    vessel_ori(isnan(vessel_ori)) = mean(vessel_ori(~isnan(vessel_ori) & (vessel_pred > 0.25)));
    vessel_patch = equalise_nailfold_intensity(vessel_struc.vessel_patch);
    
    [potential_map scale_map, rotation_map, mean_width_map] = ...
        estimate_vessel_pose(vessel_pred, vessel_width, vessel_ori,...
        'width_m', width_m, 'width_c', width_c, 'patch_size', 101);
    
    mean_widths = mean(mean_width_map(potential_map));
    %sigmas = [0.5 0.75 1 1.25 1.5]*4;
    %mag_2d = zeros(size(vessel_pred,1), size(vessel_pred,2), 3);
    %for i_si = 1:5
    %    [mag_2d(:,:,i_si)] = gaussian_2nd_derivative_line(vessel_patch, i_si);
    %end
    %filter_scales_map = map_to_nearest_val(scale_map, [0.5 0.75 1 1.25 1.5]);
    %
    %C1 = mb_normxcorr2(apex_template.g2, mag_2d(:,:,3));
    %C2 = targetted_template_matching(...
    %    mag_2d, potential_map, scale_map, rotation_map, apex_template.g2,...
    %    'filter_scales_map', filter_scales_map);
    
    mag_2d = gaussian_2nd_derivative_line(vessel_patch, 4);
    C1 = mb_normxcorr2(apex_template.g2, mag_2d);
    C2 = targetted_template_matching(...
        mag_2d, potential_map, scale_map, rotation_map, apex_template.g2,...
        'filter_scales_map', []);
    C1_masked = C1;
    C1_masked(~potential_map) = 0;
    
    if 1%i_ve <= 10
        apex_xy = unaligned_shape_data.apex_shapes(i_ve, [16 47]);
    
        figure; 
        subplot(1,3,1); imgray(C1); colormap(jet(256)); caxis([0 1]);
        plot(apex_xy(1), apex_xy(2), 'kx');
        subplot(1,3,2); imgray(C1_masked); colormap(jet(256)); caxis([0 1]);
        plot(apex_xy(1), apex_xy(2), 'kx');
        subplot(1,3,3); imgray(C2); colormap(jet(256)); caxis([0 1]);
        plot(apex_xy(1), apex_xy(2), 'kx');
    end
    
    
    %measured_dists(i_ve) = mean_dist_map(v_apex_y, v_apex_x);
end
%%
load('C:\isbe\nailfold\data\rsa_study\models\apex_templates\apex_templates.mat');
load('C:\isbe\nailfold\data\rsa_study\models\apex_templates\estimated_width_transforms.mat');
for i_ve = 21:40%81:102

    vessel_struc = u_load(vessel_files{i_ve});
    vessel_pred = u_load(vessel_prob_files{i_ve});
    vessel_ori = u_load(vessel_ori_files{i_ve});
    vessel_width = u_load(vessel_width_files{i_ve});
    
    vessel_pred(isnan(vessel_pred)) = 1;
    vessel_width(isnan(vessel_width)) = mean(vessel_width(~isnan(vessel_width) & (vessel_pred > 0.25)));
    vessel_ori(isnan(vessel_ori)) = mean(vessel_ori(~isnan(vessel_ori) & (vessel_pred > 0.25)));
    vessel_patch = equalise_nailfold_intensity(vessel_struc.vessel_patch);
    
    [potential_map scale_map, rotation_map, mean_width_map] = ...
        estimate_vessel_pose(vessel_pred, vessel_width, vessel_ori,...
        'width_m', width_m, 'width_c', width_c, 'patch_size', 101);
    
    C1 = mb_normxcorr2(apex_template.vessel_prob, vessel_pred);
    C2 = targetted_template_matching(...
        vessel_pred, potential_map, scale_map, rotation_map, apex_template.vessel_prob);
    C1(~potential_map) = 0;
    
    if 1%i_ve <= 10
        apex_xy = unaligned_shape_data.apex_shapes(i_ve, [16 47]);
    
        figure; 
        subplot(1,2,1); imgray(C1); colormap(jet(256)); caxis([0 1]);
        subplot(1,2,2); imgray(C2); colormap(jet(256)); caxis([0 1]);
        plot(apex_xy(1), apex_xy(2), 'gx');
    end
    
    
    %measured_dists(i_ve) = mean_dist_map(v_apex_y, v_apex_x);
end
%%
orig_dir = 'C:\isbe\nailfold\data\rsa_study\test\predictions\detection\rf_classification\182321a\';
missing_dir = 'C:\isbe\nailfold\data\rsa_study\test\predictions\detection\rf_classification\182321\';
corrected_dir = 'C:\isbe\nailfold\data\rsa_study\test\predictions\detection\rf_classification\182321_corrected\';
fov_mask_dir = 'C:\isbe\nailfold\data\rsa_study\test\fov_masks\';
missing_mask_dir = 'C:\isbe\nailfold\data\rsa_study\test\missing_masks\';

image_list = dir([missing_dir '*.mat']);

for i_im = 1:length(image_list)
    if exist([orig_dir image_list(i_im).name], 'file')
        pred_image = u_load([orig_dir image_list(i_im).name]);
        pred_missing = u_load([missing_dir image_list(i_im).name]);

        fov_mask = u_load([fov_mask_dir image_list(i_im).name(1:6) '_f_mask.mat']);
        missing_mask = u_load([missing_mask_dir image_list(i_im).name(1:6) '_f_mask.mat']);

        pred_image(missing_mask) = pred_missing(missing_mask);
        pred_image(~fov_mask) = 0;
        save([corrected_dir  image_list(i_im).name], 'pred_image');
    else
        copyfile(...
            [missing_dir image_list(i_im).name],...
            [orig_dir image_list(i_im).name]);
    end
%     if i_im <= 10
%         figure; imgray(pred_image);
%     end
end
%%
for i_im = 641:length(image_list)
    copyfile(...
        [missing_dir image_list(i_im).name],...
        [orig_dir image_list(i_im).name]);
end
    

