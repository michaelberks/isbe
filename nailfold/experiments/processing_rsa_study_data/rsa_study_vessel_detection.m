%--------------------------------------------------------------------------
% ---------------- Script for detecting vessels in RSA data ---------------
%--------------------------------------------------------------------------
%%
% Steps in the process:
% 1)
% 2)
% 3)....
%%-------------------------------------------------------------------------

%% 1) Make lists of file locations used in the training data
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
%-------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% 2) Build shape model of shapes and apex templates 
clear;
data_lists_dir = 'C:\isbe\nailfold\data\rsa_study\data_lists\';
load([data_lists_dir 'training_data_lists.mat'], '*_files');
[unaligned_shape_data aligned_shape_data] = vessel_detection_shape_alignment(... 
    'contour_files', contour_files,...
    'vessel_files', vessel_files,...
    'model_dir', 'C:\isbe\nailfold\data\rsa_study\models\apex_templates\',...
    'num_shape_pts', 31,...
    'plot', 1);

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
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% 3) Model relationship between estimated and actual pose
clear;
data_lists_dir = 'C:\isbe\nailfold\data\rsa_study\data_lists\';
model_dir = 'C:\isbe\nailfold\data\rsa_study\models\apex_templates\';
load('C:\isbe\nailfold\data\rsa_study\models\apex_templates\unaligned_apex_shape_data.mat');
load('C:\isbe\nailfold\data\rsa_study\models\apex_templates\aligned_apex_shape_data.mat');
load([data_lists_dir 'training_data_lists.mat'], '*_files');

measured_widths = zeros(102,1);
measured_oris = zeros(102,1);

mean_scale = mean(aligned_shape_data.a_scales);
for i_ve = 1:102

    %Load data
    vessel_pred = u_load(vessel_prob_files{i_ve});
    vessel_ori = u_load(vessel_ori_files{i_ve});
    vessel_width = u_load(vessel_width_files{i_ve});
    
    %Check we havent got any NaNs in the data
    vessel_pred(isnan(vessel_pred)) = 1;
    vessel_width(isnan(vessel_width)) = mean(vessel_width(~isnan(vessel_width) & (vessel_pred > 0.25)));
    vessel_ori(isnan(vessel_ori)) = mean(vessel_ori(~isnan(vessel_ori) & (vessel_pred > 0.25)));
    
    %Estimate the vessel pose
    [potential_map, ~, rotation_map, mean_width_map, mean_ori_map] = ...
        estimate_vessel_pose(vessel_pred, vessel_width, vessel_ori, 'patch_size', 101);
         
    %Find the nearest vessel point to the apex
    apex_xy = unaligned_shape_data.apex_shapes(i_ve, [16 47]);
    [vy vx] = find(potential_map);
    dists = (vx-apex_xy(1)).^2 + (vy-apex_xy(2)).^2;
    [~, min_i] = min(dists);
    v_apex_x = vx(min_i);
    v_apex_y = vy(min_i);
    
    %Extract the estimated rotation and scale at this point
    measured_widths(i_ve) = mean_width_map(v_apex_y, v_apex_x);
    measured_oris(i_ve) = rotation_map(v_apex_y, v_apex_x);
    
    if ismember(i_ve, [34 68 69 70 71 87 96])
        figure; 
        subplot(1,2,1); imgray(complex2rgb(potential_map.*mean_ori_map));
        plot(apex_xy(1), apex_xy(2), 'gx');
        plot(v_apex_x, v_apex_y, 'go');
        
        subplot(1,2,2); imgray(two_map_colour(potential_map, mean_width_map, []));
        plot(apex_xy(1), apex_xy(2), 'gx');
        plot(v_apex_x, v_apex_y, 'go');
    end
    
    
end

%Display the estimated widths against the known widths
figure; plot(unaligned_shape_data.apex_widths, measured_widths, 'rx');
axis([0 70 0 45]);
hold on;
plot([0 70], [0 45], 'k--');
xlabel('Widths from MB annotation');
ylabel('Estimated widths');

%Display the estimated widths against the know tansform scales - the
%relationship we really want to model
i_scales = mean_scale ./ aligned_shape_data.a_scales;
figure; plot(measured_widths, i_scales, 'rx');
axis([0 max(measured_widths) 0 max(i_scales)]);
hold on;
plot([0 max(measured_widths)], [0 max(i_scales)], 'k--');
[reg_score_widths, width_m, width_c] = regression(measured_widths', i_scales');
plot([0 max(measured_widths)], width_m*[0 max(measured_widths)]+width_c, 'g');
xlabel('Estimated widths');
ylabel('Scales applied to apex sampling points');

%Save the relationship between estimated width and transform scale
save([model_dir 'estimated_width_transforms.mat'], 'width_m', 'width_c');

%Display the relationship between the known transform rotations and the
%estimated rotations
a_thetas = atan2(squeeze(aligned_shape_data.a_rots(1,2,:)), squeeze(aligned_shape_data.a_rots(1,1,:)));
figure; plot(measured_oris, a_thetas, 'rx');
axis equal; axis([-pi/2 pi/2 -pi/2 pi/2]); hold on;
plot([-pi/2 pi/2], [-pi/2 pi/2], 'k--');
[reg_score_oris, ori_m, ori_c] = regression(measured_oris', a_thetas');
xlabel('Estimated vessel orientation');
ylabel('Rotation applied to apex sampling points');