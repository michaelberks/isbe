model_dir = 'C:\isbe\nailfold\data\rsa_study\models\apex_templates\';
image_dir = 'C:\isbe\nailfold\data\rsa_study\images\';
mask_dir = 'C:\isbe\nailfold\data\rsa_study\test\fov_masks\';
pred_dir = 'C:\isbe\nailfold\data\rsa_study\test\predictions\detection\rf_classification\182321\';
ori_dir = 'C:\isbe\nailfold\data\rsa_study\test\predictions\orientation\rf_regression\182263\';
width_dir = 'C:\isbe\nailfold\data\rsa_study\test\predictions\width\rf_regression\182367\';
markup_dir = 'C:\isbe\nailfold\data\rsa_study\markup\tmoore\';

load([model_dir 'estimated_width_transforms.mat'], 'width_m', 'width_c');
%%
i_list = dir([image_dir '*.mat']);
p_list = dir([pred_dir '*.mat']);
o_list = dir([ori_dir '*.mat']);
w_list = dir([width_dir '*.mat']);

max_r = 800;
max_c = 1200;
%%
i_im = 1;

im_name = o_list(i_im).name(1:end-9); 

vessel_markup_list = dir([markup_dir '*' im_name '*.txt']); 

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

g = gaussian_filters_1d(2);
g = g / sum(g);
vessel_pred_smooth = conv2(g', g, vessel_pred, 'same');

%Non-maximally supress vessel probs to locate potential centre lines
vessel_nms = mb_non_maximal_supp(vessel_pred_smooth, angle(vessel_ori)/2);

%Perform hysteris on the centrelines
strong_vessels = vessel_nms > 0.25;
if any(strong_vessels(:))
    [rstrong cstrong] = find(strong_vessels);
    potential_map = bwselect(vessel_nms > 0, cstrong, rstrong, 8);
else
    potential_map = strong_vessels;
end

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
            end
        end
    else
        %plot the anchor
        plot(a1, -off_x+anchor_xy(1), -off_y+anchor_xy(2), 'g+');
        plot(a2, -off_x+anchor_xy(1), -off_y+anchor_xy(2), 'g+');
    end
end
%%
orientation_map = angle(vessel_ori)/2;
feature_map = vessel_ori.*vessel_pred;
selection_map = potential_map;
weighting_map = vessel_pred;
%%
r_idx = 312;
c_idx = 448;

[py px] = find(selection_map);
%%
[model_sc_c] = compute_shape_contexts(feature_map, orientation_map, selection_map, r_idx, c_idx, 'weighting_map', []);
[model_sc_p model_pc_p] = shape_context_prob_track_mult(feature_map, vessel_pred_smooth, vessel_ori, r_idx, c_idx, 'weighting_map', []);

[sc_ci] = compute_shape_contexts(feature_map, orientation_map, selection_map, 253, 75);
[sc_pi pc_pi] = shape_context_prob_track_mult(feature_map, vessel_pred_smooth, vessel_ori, 253, 75, 'num_streams', 1e4);
%
figure; imgray(display_shape_context(model_sc_c, 10));
figure; imgray(display_shape_context(model_sc_p, 10));
figure; imgray(display_shape_context(sc_ci, 10));
figure; imgray(display_shape_context(sc_pi, 10));
%%
potential_skel = bwmorph(potential_map, 'skel', 'inf');
potential_ends = bwmorph(potential_skel, 'endpoints');
[ey ex] = find(potential_ends);
figure; imgray(potential_map);
plot(ex, ey, 'rx');
%%
potential_labels = bwlabel(potential_skel, 8);
potential_end_labels = potential_labels;
potential_end_labels(~potential_ends) = 0;
potential_end_labels = imdilate(potential_end_labels, strel('disk', 5));
label_max = max(potential_end_labels(:));

figure; imgray(potential_end_labels); colormap([1 1 1; lines(label_max)]);
%%
stream_connect_segments(vessel_pred_smooth, orientation_map, potential_skel, 1, 1e3, 20, 5);
%%
[sc_p] = shape_context_prob_track_mult(feature_map, vessel_pred_smooth, vessel_ori, py, px, 'num_streams', 1e3);
save C:\isbe\nailfold\shape_context_prob_tracks.mat sc_p px py ori_c
%%
pidx = sub2ind(size(vessel_ori),py,px);
ori_c = vessel_ori(pidx);
theta_i = angle(ori_c) / 2;
for i_p = 1:length(pidx)
    sc_p(:,:,i_p) = sc_p(:,:,i_p) ./ exp(-2*theta_i(i_p));
    sc_p(:,:,i_p) = sc_p(:,:,i_p) .* exp(-2i*theta_i(i_p));
end

%%
i_pt = find( (px==566)&(py==248) );
figure; imgray(display_shape_context(sc_p(:,:,i_pt), 10));
%%
for i_im = 1%:20

    im_name = o_list(i_im).name(1:end-9); 

    vessel_markup_list = dir([markup_dir '*' im_name '*.txt']); 

    vessel_markup = read_markup_from([markup_dir vessel_markup_list(end).name]);
    vessel_im = imread([image_dir im_name '.png']);

    [rows cols] = size(vessel_im);
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

    vessel_im = vessel_im(keep_rows, keep_cols);
    vessel_im_eq = equalise_nailfold_intensity(vessel_im);
    
    
    [maxima_pos, maxima_vals, corr_map] = template_match_apexes(vessel_im_eq, ...
        {'g2d', apex_template.g2}, 'threshold', 0.25);
    
    figure; imgray(vessel_im_eq);
    plot(maxima_pos(:,1), maxima_pos(:,2), 'cx');
    
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
                    plot(-off_x+anchor_xy(1), -off_y+anchor_xy(2), 'go');
                else
                    %plot the apex
                    apex_xy = ...
                        [ vessel_markup.vessels(i_v).apices(i_a).outer_point;...
                          vessel_markup.vessels(i_v).apices(i_a).inner_point];
                    plot(-off_x+apex_xy(:,1), -off_y+apex_xy(:,2), 'r-x');
                end
            end
        else
            %plot the anchor
            plot(-off_x+anchor_xy(1), -off_y+anchor_xy(2), 'g+');
        end
    end
end
