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

max_r = 900;
max_c = 900;
%
for i_im = 1:20%length(o_list)
    
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
    
    [potential_map scale_map, rotation_map] = ...
        estimate_vessel_pose(vessel_pred, vessel_width, vessel_ori,...
        'width_m', width_m, 'width_c', width_c, 'patch_size', 101);
    
    potential_map = imdilate(potential_map, strel('disk', 2));
    
    vessel_pred_nms = vessel_pred;
    vessel_pred_nms(~potential_map) = 0;
    
    gp = prctile(double(vessel_im(:)), 0:100);
    [~, ip] = max(diff(gp(51:end)));
    g_min = gp(2);
    g_max = gp(ip+50);
    
    figure; 
    %fpos = get(f1, 'position');
    a1 = subplot('Position',[0 0 .5 1], 'units', 'normalized');% subplot(1,2,1); 
    imgray(vessel_im); caxis([g_min g_max]);
    a2 = subplot('Position',[.5 0 .5 1], 'units', 'normalized');%subplot(1,2,2); 
    imgray(vessel_pred_nms);
    linkaxes([a1 a2]);
    
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
model_dir = 'C:\isbe\nailfold\data\rsa_study\models\apex_templates\';
image_dir = 'C:\isbe\nailfold\data\rsa_study\images\';
mask_dir = 'C:\isbe\nailfold\data\rsa_study\test\fov_masks\';
pred_dir = 'C:\isbe\nailfold\data\rsa_study\test\predictions\detection\rf_classification\182321\';
ori_dir = 'C:\isbe\nailfold\data\rsa_study\test\predictions\orientation\rf_regression\182263\';
width_dir = 'C:\isbe\nailfold\data\rsa_study\test\predictions\width\rf_regression\182367\';
markup_dir = 'C:\isbe\nailfold\data\rsa_study\markup\tmoore\';

load([model_dir 'estimated_width_transforms.mat'], 'width_m', 'width_c');
load('C:\isbe\nailfold\data\rsa_study\models\apex_templates\apex_templates.mat');

i_list = dir([image_dir '*.mat']);
p_list = dir([pred_dir '*.mat']);
o_list = dir([ori_dir '*.mat']);
w_list = dir([width_dir '*.mat']);
%%
nearest_nms = [];
nearest_nms_vals = [];
nearest_c_vals = [];

max_r = 900;
max_c = 900;
g = gaussian_filters_1d(2);
%
for i_im = 19%1:20%length(o_list)
    
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
    vessel_pred = conv2(g', g, vessel_pred, 'same');
    
    vessel_ori = u_load([ori_dir o_list(i_im).name]);
    vessel_ori = vessel_ori(keep_rows, keep_cols);
    vessel_width = u_load([width_dir o_list(i_im).name]);
    vessel_width = vessel_width(keep_rows, keep_cols);
    
    vessel_im = imread([image_dir im_name '.png']);
    vessel_im = vessel_im(keep_rows, keep_cols);
    
    gp = prctile(double(vessel_im(:)), 0:100);
    [~, ip] = max(diff(gp(51:end)));
    g_min = gp(2);
    g_max = gp(ip+50);
    
    figure; 
    %fpos = get(f1, 'position');
    a1 = subplot('Position',[0 0 .5 1], 'units', 'normalized');% subplot(1,2,1); 
    imgray(vessel_im); caxis([g_min g_max]);
    
    a2 = subplot('Position',[.5 0 .5 1], 'units', 'normalized');%subplot(1,2,2); 
    imgray(complex2rgb(vessel_pred.*exp(1i*angle(vessel_ori))));
    linkaxes([a1 a2]);

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
num_images = length(o_list);
vessel_curv_prctls = zeros(num_images, 101);
%%
vessel_curv_mm = zeros(num_images, 2);
for i_im = 1:num_images    
    display(['Processing image ' num2str(i_im)]);
    im_name = o_list(i_im).name(1:6);
    vessel_ori = u_load([ori_dir o_list(i_im).name]);
    vessel_mask = u_load([mask_dir im_name '_f_mask.mat']);
    %vessel_mask(1:2:end, 1:2:end) = 0;
    
    vessel_curv = abs(complex_gaussian_curvature(vessel_ori, 2));
    %vessel_curv_prctls(i_im,:) = prctile(vessel_curv(vessel_mask), 0:100);
    vessel_curv_mm(i_im,1) = min(vessel_curv(vessel_mask));
    vessel_curv_mm(i_im,2) = max(vessel_curv(vessel_mask));
end