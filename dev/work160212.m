
vessel_pred = imread('N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\066wellcome\2015_04_15\L4_14_07_09\capillary_data\vessels_v_pred.png');
vessel_ori = imread('N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\066wellcome\2015_04_15\L4_14_07_09\capillary_data\vessels_o_pred.png');
vessel_width = imread('N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\066wellcome\2015_04_15\L4_14_07_09\capillary_data\vessels_w_pred.png');
load('C:\isbe\nailfold\data\wellcome_study\capillary_data\066wellcome_2015_04_15_L4_14_07_09_capillary_data.mat');

[sequence_data] = read_processed_sequence_from('C:\isbe\nailfold\data\wellcome_study\capillary_example\sequence_data.dat');
mosaic_w = sequence_data.final_mosaic_size(1);
mosaic_h = sequence_data.final_mosaic_size(2);

vessel_pred = imresize(vessel_pred, [mosaic_h mosaic_w]);
vessel_ori = imresize(vessel_ori, [mosaic_h mosaic_w]);
vessel_width = imresize(vessel_width, [mosaic_h mosaic_w]);
%%
load('N:\Nailfold Capillaroscopy\wellcome\flow_data\066wellcome_2015_04_15_L4_14_07_09_s02_v04.mat');
load('N:\Nailfold Capillaroscopy\wellcome\flow_results\066wellcome_2015_04_15_L4_14_07_09_s02_v04.mat');

[~, ~, dx, dy] = ...
    mosaic_limits(size(flow_results.flowPyramidEst{1}), vessel_transforms);
    
x_start = floor(x_min + find(any(~edge_mask,1),1) - (dx + 1));
y_start = floor(y_min + find(any(~edge_mask,2),1) - (dy + 1));
    
[flow_h flow_w] = size(flow_results.flowPyramidEst{1});
cropped_rows = y_start+(1:flow_h);
cropped_cols = x_start+(1:flow_w);
  
%Now compute information about these frames - first get vessel and
%orientation patches for the correctly translated bundign box
pred_rows = cropped_rows+floor(frames_offset(2));
pred_cols = cropped_cols+floor(frames_offset(1));
valid_pred_rows = pred_rows >= 1 & pred_rows <= size(vessel_ori,1);
valid_pred_cols = pred_cols >= 1 & pred_cols <= size(vessel_ori,2);
    
frame_ori = zeros(flow_h, flow_w, 3);
frame_ori(valid_pred_rows, valid_pred_cols, :) = vessel_ori(...
    pred_rows(valid_pred_rows),...
    pred_cols(valid_pred_cols),:);
frame_ori = rgb2complex(frame_ori, [], 1, [], 0);
frame_pred = zeros(flow_h, flow_w);
frame_pred(valid_pred_rows, valid_pred_cols) = double(vessel_pred(...
    pred_rows(valid_pred_rows),...
    pred_cols(valid_pred_cols))) / 100;

frame_width = zeros(flow_h, flow_w);
frame_width(valid_pred_rows, valid_pred_cols) = double(vessel_width(...
    pred_rows(valid_pred_rows),...
    pred_cols(valid_pred_cols)));

%Make a frame mask from the vessel predictions, selected pixels
%connected to the vessel apex
ax = apex_measures.distal.apex_xy(apex_idx,1)/resize_factor ...
    - (floor(frames_offset(1)) + x_start);
ay = apex_measures.distal.apex_xy(apex_idx,2)/resize_factor ...
    - (floor(frames_offset(2)) + y_start);

connect_thresh = 0.5;
frame_vessels_mask = frame_pred > connect_thresh;
frame_apex_mask = bwselect(frame_vessels_mask, ax, ay, 4);

g_prob = gaussian_filters_1d(2);
g_prob = g_prob / sum(g_prob);    
smooth_pred = conv2(g_prob', g_prob, frame_pred, 'same');

vessel_centre = mb_non_maximal_supp(smooth_pred, angle(frame_ori)/2);
vessel_centre(~frame_apex_mask) = 0;
vessel_centre_mask = vessel_centre > 0;
[vy vx] = find(vessel_centre);

smooth_flow = conv2(g_prob', g_prob, flow_results.flowPyramidEst{1}, 'same');
fxy = smooth_flow(vessel_centre_mask);

figure; show_flow_as('rgb', smooth_flow); hold all;
quiver(vx, vy, 10*real(fxy), 10*imag(fxy));
%%
figure;
subplot(2,3,1); imgray(frame_pred);
plot(ax, ay, 'rx');
subplot(2,3,2); imgray(frame_apex_mask);
plot(ax, ay, 'rx');
subplot(2,3,3); imgray(complex2rgb(frame_ori, [], [], [], 1));
plot(ax, ay, 'rx');
subplot(2,3,4); imgray(frame_width);
plot(ax, ay, 'rx');
subplot(2,3,5); imgray(vessel_centre);
plot(ax, ay, 'rx');
subplot(2,3,6); show_flow_as('rgb', flow_results.flowPyramidEst{1}); hold on;
plot(ax, ay, 'rx');
%%
frames = cropped_frames;
[fh fw num_frames] = size(frames);
g_lims = prctile(double(frames(:)), [1 99]);
g_range = g_lims(2) - g_lims(1);
frames = 255*(double(frames)-g_lims(1))/g_range;

frames = mb_pad(frames, [flow_h - fh, flow_w - fw, 0], 'replicate', 'post');

[~, aid] = min((vx-ax).^2 + (vy-ay).^2);
ix = vx(aid);
iy = vy(aid);
cell_radius2 = (frame_width(iy, ix) * 0.375).^2;

xx = repmat(1:flow_w, flow_h, 1);
yy = repmat((1:flow_h)', 1, flow_w);

temp_frames_dir = tempname;
mkdir(temp_frames_dir);
for i_fr = 1:num_frames;
    
    cell_mask = (xx-ix).^2 + (yy-iy).^2 < cell_radius2;
    frame_i = frames(:,:,i_fr);
    frame_i(cell_mask) = 255;    
    imwrite(uint8(frame_i),...
        [temp_frames_dir '\frame' zerostr(i_fr,4) '.bmp']);
    
    if ix > 1 && ix < flow_w && iy > 1 && iy < flow_h;
        fxyi = smooth_flow(round(iy), round(ix));
        ix = ix + real(fxyi);
        iy = iy + imag(fxyi);
    else
        ix = vx(aid);
        iy = vy(aid);
    end
end

video_path = 'C:\isbe\test_vid_fwd.mp4';
frame_rate = 30;
cmd = ['ffmpeg -y -r ' num2str(frame_rate) ' -i "' temp_frames_dir '\frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' video_path '"'];
system(cmd);

delete([temp_frames_dir '\*']);

ix = vx(aid);
iy = vy(aid);

temp_frames_dir = tempname;
mkdir(temp_frames_dir);
for i_fr = 1:num_frames;
    
    cell_mask = (xx-ix).^2 + (yy-iy).^2 < cell_radius2;
    frame_i = frames(:,:,num_frames - i_fr + 1);
    frame_i(cell_mask) = 255;    
    imwrite(uint8(frame_i),...
        [temp_frames_dir '\frame' zerostr(i_fr,4) '.bmp']);
    
    if ix > 1 && ix < flow_w && iy > 1 && iy < flow_h;
        fxyi = smooth_flow(round(iy), round(ix));
        ix = ix - real(fxyi);
        iy = iy - imag(fxyi);
    else
        ix = vx(aid);
        iy = vy(aid);
    end
end

video_path = 'C:\isbe\test_vid_rev.mp4';
frame_rate = 30;
cmd = ['ffmpeg -y -r ' num2str(frame_rate) ' -i "' temp_frames_dir '\frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' video_path '"'];
system(cmd);

delete([temp_frames_dir '\*']);

rmdir(temp_frames_dir);
%%
%Set max vessel path to 2 x the diagonal of the patch
max_dist = 2*sqrt(flow_h^2 + flow_w^2);
max_steps = ceil(max_dist);%flow_h*flow_w;
art_path = zeros(max_steps,3);
ven_path = zeros(max_steps,3);

%Go forwards from apex (i.e venous path)
ven_path(1,:) = [vx(aid) vy(aid) 0];
for i_step = 2:max_steps
    ix = ven_path(i_step-1,1);
    iy = ven_path(i_step-1,2);
    dist = ven_path(i_step-1,3);
    
    if ix > 1 && ix < flow_w && iy > 1 && iy < flow_h && dist < max_dist;
        fxyi = smooth_flow(round(iy), round(ix));
        fm = abs(fxyi) + 1e-6;
        ven_path(i_step,1) = ix + real(fxyi) / fm;
        ven_path(i_step,2) = iy + imag(fxyi) / fm;
        ven_path(i_step,3) = dist + 1;
    else
        ven_path(i_step:end,:) = [];
        break;
    end
end

%Go backwards from apex (i.e arterial path)
max_dist = -max_dist;
art_path(1,:) = [vx(aid) vy(aid) 0];
for i_step = 2:max_steps
    ix = art_path(i_step-1,1);
    iy = art_path(i_step-1,2);
    dist = art_path(i_step-1,3);
    
    if ix > 1 && ix < flow_w && iy > 1 && iy < flow_h && dist > max_dist;
        fxyi = smooth_flow(round(iy), round(ix));
        fm = abs(fxyi) + 1e-6;
        art_path(i_step,1) = ix - real(fxyi) / fm;
        art_path(i_step,2) = iy - imag(fxyi) / fm;
        art_path(i_step,3) = dist - 1;
    else
        art_path(i_step:end,:) = [];
        break;
    end
end

vessel_path = [flipud(art_path); ven_path(2:end,:)];
figure; show_flow_as('rgb', flow_results.flowPyramidEst{1}); hold on;
plot(vessel_path(:,1), vessel_path(:,2), 'k--');

figure; imgray(frame_apex_mask);
plot(vessel_path(:,1), vessel_path(:,2), 'r--');
%
[all_vy all_vx] = find(frame_apex_mask);
all_vd = griddata(vessel_path(:,1), vessel_path(:,2), vessel_path(:,3), ...
    all_vx, all_vy, 'nearest');
    
dist_map = zeros(flow_h, flow_w);
dist_map(frame_apex_mask) = all_vd;

figure; imgray(dist_map);
%%
frames = cropped_frames;
[fh fw num_frames] = size(frames);
g_lims = prctile(double(frames(:)), [1 99]);
g_range = g_lims(2) - g_lims(1);
frames = 255*(double(frames)-g_lims(1))/g_range;

frames = mb_pad(frames, [flow_h - fh, flow_w - fw, 0], 'replicate', 'post');

cell_radius2 = (frame_width(round(ay), round(ax)) * 0.2).^2;
cell_radius = sqrt(cell_radius2);
first_pt = ceil(cell_radius);

num_cells = 8;
cell_idx = round(linspace(first_pt, size(vessel_path,1), num_cells+1));
cell_idx(end) = [];
cell_dists = vessel_path(cell_idx,3);
% cell_pts = zeros(num_cells,2,2);
% for i_cell = 1:num_cells
%     candidate_cells = all_vd == cell_dists(i_cell);
%     
%     vxi = all_vx(candidate_cells);
%     vyi = all_vy(candidate_cells);
%     r_idx = ceil(length(vxi)*rand);
%     cell_pts(i_cell,1,:) = vxi(r_idx);
%     cell_pts(i_cell,2,:) = vyi(r_idx);
% end
cell_pts = repmat(vessel_path(cell_idx,1:2), 1, 1, 2);

xx = repmat(1:flow_w, flow_h, 1);
yy = repmat((1:flow_h)', 1, flow_w);

temp_frames_dir = tempname;
mkdir(temp_frames_dir);
scaling = [1.0 1.5];
for i_fr = 1:num_frames;
    
%     vid_frame = zeros(flow_h, 2*flow_w + 16, 3, 'uint8'); 
    vid_frame = zeros(flow_h, flow_w, 3, 'uint8'); 
    for i_dim = 1:2
        frame_i = frames(:,:,i_fr);
        for i_cell = 1:num_cells
            ix = cell_pts(i_cell,1,i_dim);
            iy = cell_pts(i_cell,2,i_dim);
            xy2 = (xx-ix).^2 + (yy-iy).^2;
            cell_mask = xy2 < cell_radius2;% & xy2 > (cell_radius2/4);    
            frame_i(cell_mask) = 255;    

            if ix > 1 && ix < flow_w && iy > 1 && iy < flow_h;
                fxyi = smooth_flow(round(iy), round(ix));
                cell_pts(i_cell,1,i_dim) = ix + scaling(i_dim)*real(fxyi);
                cell_pts(i_cell,2,i_dim) = iy + scaling(i_dim)*imag(fxyi);
            else
%                 candidate_cells = all_vd == vessel_path(first_pt, 3);
% 
%                 vxi = all_vx(candidate_cells);
%                 vyi = all_vy(candidate_cells);
%                 r_idx = ceil(length(vxi)*rand);
%                 cell_pts(i_cell,1,i_dim) = vxi(r_idx);
%                 cell_pts(i_cell,2,i_dim) = vyi(r_idx);
                cell_pts(i_cell,:,i_dim) = vessel_path(first_pt,1:2);
            end
        end
%         vid_frame(:, (1:flow_w) + (i_dim-1)*(flow_w + 16), 1) = uint8(frame_i);
%         vid_frame(:, (1:flow_w) + (i_dim-1)*(flow_w + 16), 2) = uint8(frames(:,:,i_fr));
%         vid_frame(:, (1:flow_w) + (i_dim-1)*(flow_w + 16), 3) = uint8(frames(:,:,i_fr));
        vid_frame(:, :, i_dim) = uint8(frame_i);        
    end
    vid_frame(:, :, 3) = uint8(frames(:,:,i_fr));
    imwrite(vid_frame,...
        [temp_frames_dir '\frame' zerostr(i_fr,4) '.bmp']);
end

video_path = 'C:\isbe\test_vid_pair10.mp4';
frame_rate = 15;
cmd = ['ffmpeg -y -r ' num2str(frame_rate) ' -i "' temp_frames_dir '\frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' video_path '"'];
system(cmd);

delete([temp_frames_dir '\*']);
%%
make_flow_comparison_video('066wellcome\2015_04_15\L4_14_07_09', 02, 04, 1.5)
%%
flow_data_dir = 'N:\Nailfold Capillaroscopy\wellcome\flow_data\';
flow_results_dir = 'N:\Nailfold Capillaroscopy\wellcome\flow_results\';
capillary_data_dir = 'C:\isbe\nailfold\data\wellcome_study\capillary_data\';
load('C:\isbe\nailfold\data\wellcome_study\sequence_names.mat');
sequence_names2 = sequence_names;
for i_seq = 1:numel(sequence_names)
    
    if isempty(sequence_names2{i_seq})
        continue;
    end
    
    seq_dir = [sequence_names2{i_seq} '\'];%sequence_dirs{i_seq};
    seq_name = seq_dir;
    dividers = seq_name == '\';
    pos = find(dividers, 4, 'last');
    seq_name(dividers) = '_';
    seq_name = seq_name(pos(1)+1:end-1);
    display(['processing sequence ' num2str(i_seq) ': ' seq_name]);
    
    flow_list = dir([flow_results_dir seq_name '*.mat']);
    
    if isempty(flow_list)
        continue;
    end
    
    load([capillary_data_dir seq_name '_capillary_data.mat'],...
        'apex_measures');
    
    apex_measures.distal.flow_names = cell(size(apex_measures.distal.width_at_apex));
    for i_f = 1:length(flow_list)
        load([flow_data_dir flow_list(i_f).name], 'apex_idx');
        if isempty(apex_measures.distal.flow_names{apex_idx})
           apex_measures.distal.flow_names{apex_idx} = {flow_list(i_f).name(end-10:end-4)};
        else
           apex_measures.distal.flow_names{apex_idx}(end+1) = {flow_list(i_f).name(end-10:end-4)};
        end
    end
    
    
    save([capillary_data_dir seq_name '_capillary_data.mat'],...
        'apex_measures', '-append');
end
%%
make_flow_comparison_video('066wellcome\2015_04_15\L4_14_07_09', 02, 04, 0.5,...
    'adjacent_videos', 0, ...
    'video_dir', 'C:\isbe\nailfold\data\wellcome_study\flow_videos\flow_comparisons\n8\',...
    'num_cells', 8, ...
    'frame_rate', 30,...
    'plot', 1)
%%
make_flow_comparison_video('066wellcome\2015_04_15\L4_14_07_09', 02, 04, 1.5,...
    'adjacent_videos', 0, ...
    'video_dir', 'C:\isbe\nailfold\data\wellcome_study\flow_videos\flow_comparisons\n8\',...
    'num_cells', 8)
make_flow_comparison_video('066wellcome\2015_04_15\L4_14_07_09', 02, 04, 0.5,...
    'adjacent_videos', 1, ...
    'video_dir', 'C:\isbe\nailfold\data\wellcome_study\flow_videos\flow_comparisons\n8\',...
    'num_cells', 8)
make_flow_comparison_video('066wellcome\2015_04_15\L4_14_07_09', 02, 04, 1.5,...
    'adjacent_videos', 1, ...
    'video_dir', 'C:\isbe\nailfold\data\wellcome_study\flow_videos\flow_comparisons\n8\',...
    'num_cells', 8);
%%
make_flow_comparison_video('066wellcome\2015_04_15\L4_14_07_09', 02, 04, 1.5,...
    'adjacent_videos', 1, ...
    'video_dir', 'C:\isbe\nailfold\data\wellcome_study\flow_videos\flow_comparisons\n4\',...
    'num_cells', 4);
make_flow_comparison_video('066wellcome\2015_04_15\L4_14_07_09', 02, 04, 1.5,...
    'adjacent_videos', 1, ...
    'video_dir', 'C:\isbe\nailfold\data\wellcome_study\flow_videos\flow_comparisons\n16\',...
    'num_cells', 16);
make_flow_comparison_video('066wellcome\2015_04_15\L4_14_07_09', 02, 04, 1.5,...
    'adjacent_videos', 0, ...
    'video_dir', 'C:\isbe\nailfold\data\wellcome_study\flow_videos\flow_comparisons\n4\',...
    'num_cells', 4);
make_flow_comparison_video('066wellcome\2015_04_15\L4_14_07_09', 02, 04, 1.5,...
    'adjacent_videos', 0, ...
    'video_dir', 'C:\isbe\nailfold\data\wellcome_study\flow_videos\flow_comparisons\n16\',...
    'num_cells', 16);
%%
make_flow_comparison_video('062wellcome\2015_04_13\L4_13_58_58', 13, 02, 0.5,...
    'adjacent_videos', 0, ...
    'video_dir', 'C:\isbe\nailfold\data\wellcome_study\flow_videos\flow_comparisons\n8\',...
    'num_cells', 8, ...
    'plot', 1)
%%
make_flow_comparison_video('066wellcome\2015_04_15\L4_14_07_09', 01, 01, 0.5,...
    'adjacent_videos', 0, ...
    'video_dir', 'C:\isbe\nailfold\data\wellcome_study\flow_videos\flow_comparisons\n8\',...
    'num_cells', 8, ...
    'plot', 1)
%%
make_flow_comparison_video('002wellcome\2015_02_27\L4_11_00_04', 'all', [], 1.5,...
    'adjacent_videos', 0, ...
    'video_dir', 'C:\isbe\nailfold\data\wellcome_study\flow_videos\flow_comparisons\n8\',...
    'num_cells', 8, ...
    'plot', 1)
%%
make_flow_comparison_video('002wellcome\2015_02_27\L4_11_00_04', 11, 03, 1.5,...
    'adjacent_videos', 0, ...
    'video_dir', 'C:\isbe\nailfold\data\wellcome_study\flow_videos\flow_comparisons\n8\',...
    'num_cells', 8, ...
    'plot', 1)
%%
%%
widths = interp2(frame_width, vessel_path(:,1), vessel_path(:,2))/resize_factor;
widths(1) = widths(2);
widths(end) = widths(end-1);
%%
normal_xy = compute_spline_normals(vessel_path);

inner_edge = vessel_path(:,1:2) + 0.5*bsxfun(@times, normal_xy, widths);
outer_edge = vessel_path(:,1:2) - 0.5*bsxfun(@times, normal_xy, widths);
figure; imgray(mean(cropped_frames,3));
plot(vessel_path(:,1), vessel_path(:,2), 'r');
plot(outer_edge(:,1), outer_edge(:,2), 'b');
plot(inner_edge(:,1), inner_edge(:,2), 'g');

directions = zeros(size(outer_edge));
directions(1,:) = vessel_path(2,1:2) - vessel_path(1,1:2);
directions(2:end-1,:) = vessel_path(3:end,1:2) - vessel_path(1:end-2,1:2);
directions(end,:) = vessel_path(end,1:2) - vessel_path(end-1,1:2);

[flowmap, mask] = ...
    create_flowmap(vessel_path(:,1:2), widths, directions, [flow_h flow_w]);

%%
cell_sz = min(widths(:));% * 0.8;
n_pixels_per_cell = cell_sz^2;
n_mask_pixels = sum(mask(:)>0);
n = 10;
n_cells = ceil(n * n_mask_pixels / n_pixels_per_cell);
samples_per_frame = 1;

cloud_add = noiseonf(max(flow_h, flow_w), 1.5);
cloud_add = 32 * normim(cloud_add(1:flow_h,1:flow_w), 'stretch_fixed');

cloud_mult = noiseonf(max(flow_h, flow_w), 1.5);
cloud_mult = 0.75 * normim(cloud_mult(1:flow_h,1:flow_w), 'stretch_fixed');
cloud_mult = 1 - cloud_mult;
%%
max_flow = max(abs(smooth_flow(:)));
cell_positions = generate_cell_positions(max_flow*flowmap, mask, ...
    vessel_path(:,1:2), widths, ...
    n_cells, num_frames+1, ...
    [], samples_per_frame);
%
temp_frames_dir = tempname;
mkdir(temp_frames_dir);

f = 1;
for i_fr = 1:num_frames
    vid_frame = make_frame(cell_positions(:,:,f:f+samples_per_frame), ...
                     [flow_h flow_w], cell_sz, ...
                     [], 128 + cloud_add, 128 * cloud_mult);
                 
    if ~rem(i_fr, 10)
        figure; imgray(vid_frame);
    end

	f = f+samples_per_frame;
    
    imwrite(uint8(vid_frame),...
            [temp_frames_dir '\frame' zerostr(i_fr,4) '.bmp']);
end

frame_rate = 30;
video_path = 'C:\isbe\test_synthetic.mp4';
cmd = ['ffmpeg -y -r ' num2str(frame_rate) ' -i "' temp_frames_dir '\frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' video_path '"'];
system(cmd);
delete([temp_frames_dir '\*']);
rmdir(temp_frames_dir);
%%
s = load('C:\isbe\nailfold\temp_cap_data.mat');
figure;
subplot(1,3,1); imgray(s.frame_pred);
plot(s.ax, s.ay, 'rx');
subplot(1,3,2); imgray(complex2rgb(s.frame_ori));
plot(s.ax, s.ay, 'kx');
subplot(1,3,3); imgray(complex2rgb(s.smooth_flow));
plot(s.ax, s.ay, 'kx');
%%
frame_apex_mask = bwselect(s.frame_pred > 0.5, s.ax, s.ay);
frame_theta = -angle(s.frame_ori) / 2;
frame_ori_half = complex(cos(frame_theta), sin(frame_theta));
frame_ori_half(~frame_apex_mask) = complex(0,0);
angle_diff = angle(frame_ori_half .* conj(s.smooth_flow));
angle_diff(~frame_apex_mask) = 0;

figure; 
subplot(1,3,1); imgray(complex2rgb(frame_ori_half)); colorbar;
subplot(1,3,2); imgray(angle_diff); colorbar;
subplot(1,3,3); imgray(abs(angle_diff) > pi/2); colorbar;
%%
figure; 
a1 = subplot(1,2,1); show_flow_as('quiver', s.frame_ori); axis ij;
a2 = subplot(1,2,2); show_flow_as('quiver', s.smooth_flow.^2); axis ij;
linkaxes([a1 a2]);
%%
figure; 
a1 = subplot(1,3,1); show_flow_as('quiver', frame_ori_half); axis ij;
a2 = subplot(1,3,2); show_flow_as('quiver', s.smooth_flow); axis ij;
a3 = subplot(1,3,3); show_flow_as('quiver', frame_ori_half .* conj(s.smooth_flow)); axis ij;
linkaxes([a1 a2 a3]);
%%
swap_mask = abs(angle_diff) > pi/2;
frame_theta2 = frame_theta;
frame_theta2(swap_mask) = frame_theta2(swap_mask) + pi;
frame_flow = complex(cos(frame_theta2), sin(frame_theta2));
frame_flow(~frame_apex_mask) = complex(0,0);
figure; 
subplot(1,2,1); imgray(complex2rgb(frame_flow));
subplot(1,2,2); imgray(complex2rgb(s.smooth_flow));
%%
[flow_h flow_w] = size(s.smooth_flow);
step_length = 2;
max_dist = 2*sqrt(flow_h^2 + flow_w^2);
max_steps = ceil(max_dist);%flow_h*flow_w;
art_path = zeros(max_steps,3);
ven_path = zeros(max_steps,3);

%Go forwards from apex (i.e venous path)
available_x = centre_x;
available_y = centre_y;
[~, aid] = min((available_x-s.ax).^2 + (available_y-s.ay).^2);
available_x(aid) = [];
available_y(aid) = [];
fxy_curr = complex(0,0);

ven_path(1,:) = [available_x(aid) available_y(aid) 0];
for i_step = 2:max_steps
    ix = ven_path(i_step-1,1);
    iy = ven_path(i_step-1,2);
    dist = ven_path(i_step-1,3);
    
    if ix > 1 && ix < flow_w && iy > 1 && iy < flow_h && dist < max_dist && ...
        frame_apex_mask(round(iy), round(ix))
        
        fxy_new = frame_flow1(round(iy), round(ix));
        if i_step == 2 || abs(angle(conj(fxy_new)*fxy_curr)) < pi/3;
            fxy_curr = fxy_new;
            available_xi = available_x;
            available_yi = available_y;
            snap = false;
            
        else
            fxy_new = frame_flow2(round(iy), round(ix));
            if i_step == 2 || abs(angle(conj(fxy_new)*fxy_curr)) < pi/3;
                fxy_curr = fxy_new;
                available_xi = available_x;
                available_yi = available_y;
                snap = false;
                display(['v ' num2str(i_step) ': using alternative direction!']); 
            else
                available_xi = available_x;
                available_yi = available_y;
                snap = false;
                display(['v ' num2str(i_step) ': neither direction works!']);
            end
        end
        jx = ix + step_length * real(fxy_curr) / fm;
        jy = iy + step_length * imag(fxy_curr) / fm;
        
        if snap            
            [~, aid] = min((available_xi-jx).^2 + (available_yi-jy).^2);

            ven_path(i_step,1) = available_xi(aid);
            ven_path(i_step,2) = available_yi(aid);
        else
            ven_path(i_step,1) = jx;
            ven_path(i_step,2) = jy;
        end
        ven_path(i_step,3) = dist + sqrt(sum(diff(ven_path(i_step+[-1 0],1:2)).^2,2));
        
        available_x(aid) = [];
        available_y(aid) = [];
    else
        ven_path(i_step:end,:) = [];
        break;
    end
end

%Go backwards from apex (i.e arterial path)
max_dist = -max_dist;
available_x = centre_x;
available_y = centre_y;
[~, aid] = min((available_x-s.ax).^2 + (available_y-s.ay).^2);
art_path(1,:) = [available_x(aid) available_y(aid) 0];
available_x(aid) = [];
available_y(aid) = [];
fxy_curr = complex(0,0);

for i_step = 2:max_steps
    ix = art_path(i_step-1,1);
    iy = art_path(i_step-1,2);
    dist = art_path(i_step-1,3);
    
    if ix > 1 && ix < flow_w && iy > 1 && iy < flow_h && dist > max_dist && ...
        frame_apex_mask(round(iy), round(ix))
        
        fxy_new = frame_flow1(round(iy), round(ix));
        if i_step == 2 || abs(angle(conj(fxy_new)*fxy_curr)) < pi/3;
            fxy_curr = fxy_new;
            available_xi = available_x;
            available_yi = available_y;
            snap = false;        
        else
            fxy_new = frame_flow2(round(iy), round(ix));
            if i_step == 2 || abs(angle(conj(fxy_new)*fxy_curr)) < pi/3;
                fxy_curr = fxy_new;
                available_xi = available_x;
                available_yi = available_y;
                snap = false;
                display(['a ' num2str(i_step) ': using alternative direction!']); 
            else
                available_xi = available_x;
                available_yi = available_y;
                snap = false;
                display(['a ' num2str(i_step) ': neither direction works!']);
            end
        end
        fm = abs(fxy_curr) + 1e-6;
        
        jx = ix - step_length * real(fxy_curr) / fm;
        jy = iy - step_length * imag(fxy_curr) / fm;
        
        if snap
            [~, aid] = min((available_xi-jx).^2 + (available_yi-jy).^2);

            art_path(i_step,1) = available_xi(aid);
            art_path(i_step,2) = available_yi(aid);
        else
            art_path(i_step,1) = jx;
            art_path(i_step,2) = jy;
        end
        art_path(i_step,3) = dist - sqrt(sum(diff(art_path(i_step+[-1 0],1:2)).^2,2));
        
        available_x(aid) = [];
        available_y(aid) = [];
    else
        art_path(i_step:end,:) = [];
        break;
    end
end

vessel_path = [flipud(art_path); ven_path(2:end,:)];
figure; show_flow_as('rgb', frame_flow); hold on;
plot(vessel_path(:,1), vessel_path(:,2), 'k--');
%%
g_prob = gaussian_filters_1d(2);
g_prob = g_prob / sum(g_prob);    
smooth_pred = conv2(g_prob', g_prob, s.frame_pred, 'same');

frame_nms = mb_non_maximal_supp(smooth_pred, -frame_theta);
[centre_y centre_x] = find(frame_nms > 0 & frame_apex_mask);
figure; show_flow_as('rgb', frame_flow); hold on;
plot(centre_x, centre_y, 'k.');
%%
flow_avg = imfilter(frame_flow, fspecial('disk', 8));
junction_mask = abs(flow_avg) < 0.5 & frame_apex_mask;
junction_mask2 = imopen(junction_mask, strel('disk', 4));

conn_comp = bwconncomp(junction_mask2);

figure;
subplot(1,2,1); hist(angle(frame_flow(conn_comp.PixelIdxList{1})), 36);
subplot(1,2,2); hist(angle(frame_flow(conn_comp.PixelIdxList{2})), 36);

xx = repmat(1:flow_w, flow_h, 1);
yy = repmat((1:flow_h)', 1, flow_w);

xm1 = mean(xx(conn_comp.PixelIdxList{1}));
xm2 = mean(xx(conn_comp.PixelIdxList{2}));
ym1 = mean(yy(conn_comp.PixelIdxList{1}));
ym2 = mean(yy(conn_comp.PixelIdxList{2}));
%%
vessel_dir = 'C:\isbe\nailfold\data\rsa_study\set12g_half\vessel_contours\';
mask_dir = 'C:\isbe\nailfold\data\rsa_study\set12g_half\vessel_masks\';
ori_dir = 'C:\isbe\nailfold\data\rsa_study\set12g_half\orientations3d\';
create_folder(ori_dir);

vessel_list = dir([vessel_dir '*.mat']);
for i_ve = 1:length(vessel_list)
    display(['Processing vessel ' num2str(i_ve)]);
    vessel_name = vessel_list(i_ve).name(1:end-7);
    mask = u_load([mask_dir vessel_name '_v_mask.mat']);
    mask_90 = u_load([mask_dir '90_' vessel_name '_v_mask.mat']);
    t = load([vessel_dir vessel_list(i_ve).name]);
    
    [~, ~, ~, ~, ori_map_3d] =...
        make_mask_from_contour(t.outer_edge, t.inner_edge, size(a.prediction_image,1), size(a.prediction_image,2));
    
    save([ori_dir vessel_name '_ori.mat'], 'ori_map_3d');
    
    mask_90a = rot90(mask, 1);
    mask_90b = rot90(mask, -1);
    if sum(mask_90(:) & mask_90a(:)) > sum(mask_90(:) & mask_90b(:))
        ori_map_3d = rot90(ori_map_3d, 1);
    else
        ori_map_3d = rot90(ori_map_3d, -1);
    end
    ori_map_3d = -ori_map_3d;
    save([ori_dir '90_' vessel_name '_ori.mat'], 'ori_map_3d');
end
%%
num_apexes = zeros(450,1);
for i_ve = 1:length(vessel_list)
    t = load([vessel_dir vessel_list(i_ve).name], 'apex_idx');
    num_apexes(i_ve) = length(t.apex_idx);

end
    
%%
t = load('C:\isbe\nailfold\data\rsa_study\set12g_half\vessel_contours\enlargedapex0015_vessel_vc.mat');
v_mask = u_load('C:\isbe\nailfold\data\rsa_study\set12g_half\vessel_masks\enlargedapex0015_vessel_v_mask.mat');
[flow_h flow_w] = size(v_mask);
flow_h = flow_h - rem(flow_h,4);
flow_w = flow_w - rem(flow_w,4);
widths = sqrt(sum((t.outer_edge - t.inner_edge).^2,2));
w1 = mean(widths(1:t.apex_idx));
w2 = mean(widths(t.apex_idx:end));

if w1 > w2
    widths = widths(end:-1:1,:);
    t.vessel_centre = t.vessel_centre(end:-1:1,:);
    t.outer_edge = t.outer_edge(end:-1:1,:);
    t.inner_edge = t.inner_edge(end:-1:1,:);
end
directions = zeros(size(t.vessel_centre));
directions(1,:) = t.vessel_centre(2,1:2) - t.vessel_centre(1,1:2);
directions(2:end-1,:) = t.vessel_centre(3:end,1:2) - t.vessel_centre(1:end-2,1:2);
directions(end,:) = t.vessel_centre(end,1:2) - t.vessel_centre(end-1,1:2);

[flowmap, mask] = ...
    create_flowmap(t.vessel_centre(:,1:2), widths, directions, [flow_h flow_w]);
%%
cell_sz = min(widths(:));% * 0.8;
n_pixels_per_cell = cell_sz^2;
n_mask_pixels = sum(mask(:)>0);
n = 10;
n_cells = ceil(n * n_mask_pixels / n_pixels_per_cell);
samples_per_frame = 1;

cloud_add = noiseonf(max(flow_h, flow_w), 1.5);
cloud_add = 32 * normim(cloud_add(1:flow_h,1:flow_w), 'stretch_fixed');

cloud_mult = noiseonf(max(flow_h, flow_w), 1.5);
cloud_mult = 0.75 * normim(cloud_mult(1:flow_h,1:flow_w), 'stretch_fixed');
cloud_mult = 1 - cloud_mult;
%%
max_flow = 3;
num_frames = 240;
cell_positions = generate_cell_positions(max_flow*flowmap, mask, ...
    t.vessel_centre, widths, ...
    n_cells, num_frames+1, ...
    [], samples_per_frame);
%
temp_frames_dir = tempname;
mkdir(temp_frames_dir);

f = 1;
for i_fr = 1:num_frames
    vid_frame = make_frame(cell_positions(:,:,f:f+samples_per_frame), ...
                     [flow_h flow_w], cell_sz, ...
                     [], 128 + cloud_add, 128 * cloud_mult);
                 
    if ~rem(i_fr, 10)
        figure; imgray(vid_frame);
    end

	f = f+samples_per_frame;
    
    imwrite(uint8(vid_frame),...
            [temp_frames_dir '\frame' zerostr(i_fr,4) '.bmp']);
end

frame_rate = 30;
video_path = 'C:\isbe\test_synthetic2.mp4';
cmd = ['ffmpeg -y -r ' num2str(frame_rate) ' -i "' temp_frames_dir '\frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' video_path '"'];
system(cmd);
delete([temp_frames_dir '\*']);
rmdir(temp_frames_dir);
%%
t = load('C:\isbe\nailfold\data\rsa_study\set12g_half\vessel_contours\enlargedapex0472_vessel_vc.mat');
[flowmap, mask, vessel_centre, widths] = ...
    create_flowmap_profile(t.vessel_centre, t.inner_edge, t.outer_edge, t.apex_idx, []);

[flow_h flow_w n_layers] = size(flowmap);
max_flow = 3;
num_frames = 240;
cell_sz = min(widths(:));% * 0.8;
n_pixels_per_cell = cell_sz^2;
n_mask_pixels = sum(mask(:)>0);
n = 10;
n_cells = ceil(n * n_mask_pixels / n_pixels_per_cell);
samples_per_frame = 1;

cell_positions = generate_cell_positions_layers(max_flow*flowmap, mask, ...
    vessel_centre, widths, ...
    n_cells, num_frames+1, ...
    [], samples_per_frame);
%%
%
temp_frames_dir = tempname;
mkdir(temp_frames_dir);
%
f = 1;
for i_fr = 1:num_frames
    vid_frame = make_frame(cell_positions(:,:,f:f+samples_per_frame), ...
                     [flow_h flow_w], cell_sz, ...
                     [], 128 + cloud_add, 128 * cloud_mult);
                 
    if ~rem(i_fr, 10)
        figure; imgray(vid_frame);
    end

	f = f+samples_per_frame;
    
    imwrite(uint8(vid_frame),...
            [temp_frames_dir '\frame' zerostr(i_fr,4) '.bmp']);
end

frame_rate = 30;
video_path = 'C:\isbe\test_synthetic3.mp4';
cmd = ['ffmpeg -y -r ' num2str(frame_rate) ' -i "' temp_frames_dir '\frame%04d.bmp" -c:v libx264 -preset slow -crf 18 -an "' video_path '"'];
system(cmd);
delete([temp_frames_dir '\*']);
rmdir(temp_frames_dir);
%%
make_synthetic_flow_comparison_video('enlargedapex0472_vessel', 3, 1.5);






