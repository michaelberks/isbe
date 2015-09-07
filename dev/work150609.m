

root_dir = nailfoldroot(1);

session_dir = ['N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\051wellcome\2015_03_30\'];
sequences_s = dir(session_dir);

num_seqs = length(sequences_s) - 2;
sequence_dir = cell(num_seqs,1);
for i_s = 1:num_seqs
    sequence_dir{i_s} = sequences_s(i_s+2).name;
end
%%
num_seqs = length(sequence_dir);
%sequences = cell(num_seqs,1);
for i_s = 1:num_seqs
    sequence = read_sequence_from([session_dir sequence_dir{i_s} '\sequence_properties.txt']);
    save([session_dir sequence_dir{i_s} '\sequence.mat'], 'sequence');
end
%%
root_dir = 'N:\nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\';
session_dir = '002wellcome\2015_02_27\';
sequence_dir = 'L4_11_00_04\';
create_folder([root_dir session_dir sequence_dir 'segments']);

sequence = u_load([root_dir session_dir sequence_dir 'sequences.mat']);

pixels_per_mm = 1000;

[segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness] = ...
    get_stationary_segments(sequence);

num_segs = length(segments_s);
discard = false(num_segs,1);
for i_seg = 1:num_segs    
    discard(i_seg) = length(segments_s{i_seg}) < 120;
end
segments_s(discard) = [];
num_segs = length(segments_s);
%    
hw = 320 / pixels_per_mm;
hh = 240 / pixels_per_mm;
frame_rect_x = [-hw hw hw -hw -hw];
frame_rect_y = [-hh -hh hh hh -hh];

figure; hold all; axis equal; a1 = gca;
plot(a1, motor_x, motor_y, 'b.', 'markersize', 4);
%
for i_seg = 1:num_segs

    frame_idx_i = segments_s{i_seg};
    num_frames = length(frame_idx_i);

    seg_x = motor_x(frame_idx_i);
    seg_y = motor_y(frame_idx_i);
    seg_z = motor_z(frame_idx_i);
    seg_s = sharpness(frame_idx_i);

    cx = median(seg_x);
    cy = median(seg_y);

    plot(a1, cx, cy, 'rx');
    plot(a1, cx+frame_rect_x, cy+frame_rect_y, 'r', 'linewidth', 2);  
end
%%
mosaic_lims_x = -310:309;
mosaic_lims_y = -230:229;

discard = false(num_segs,1);
for i_seg = 1:num_segs    
    discard(i_seg) = length(segments_s{i_seg}) < 120;
end
segments_s(discard) = [];
num_segs = length(segments_s);
   
im_folder = [root_dir session_dir sequence_dir];


frame_centres = zeros(num_segs, 2);
matched_counts = cell(num_segs,1);
segment_frames = zeros(length(mosaic_lims_y), length(mosaic_lims_x), num_segs);
segment_mosaics = cell(num_segs,1);
diff_imgs = zeros(480, 640, num_segs);
%
for i_seg = 1:num_segs
    
    segment_dir = [root_dir session_dir sequence_dir 'segments\segment' zerostr(i_seg,2)];
    create_folder(segment_dir);

    frame_idx = segments_s{i_seg};
    num_frames = length(frame_idx);

    seg_x = motor_x(frame_idx);
    seg_y = motor_y(frame_idx);
    seg_z = motor_z(frame_idx);
    seg_s = sharpness(frame_idx);

    frame_centres(i_seg,1) = seg_x(round(num_frames/2));
    frame_centres(i_seg,2) = seg_y(round(num_frames/2));

    [best_s, best_frames] = sort(seg_s, 'descend');

    num_frames = min(30,round(num_frames/4));

    frames = zeros(480,640,num_frames);
    for i_im = 1:num_frames
        frames(:,:,i_im) = rot90(imread([im_folder 'frame' zerostr(frame_idx(best_frames(i_im)), 5) '.bmp']),2);
    end

    [compound_transforms, matched_counts{i_seg}] = ...
        register_tiles_features(frames, ...
                            'theta_range', [0], ...
                            'offset_lim', 40, ...
                            'debug', false);

    [segment_mosaics{i_seg},~,compound_transforms] = ...
        create_mosaic(frames, compound_transforms);



    [r, c] = size(segment_mosaics{i_seg});
    cx = round(c/2);
    cy = round(r/2);
    segment_frames(:,:,i_seg) = segment_mosaics{i_seg}(cy+mosaic_lims_y, cx+mosaic_lims_x);

    [diff_imgs(:,:,i_seg)] = write_trans_tiles(frames, compound_transforms, ...
                             nan, [], [], segment_mosaic{i_seg});

    figure; 
    subplot(1,2,1); imgray(segment_frames(:,:,i_seg));
    subplot(1,2,2); imgray(diff_imgs(:,:,i_seg));
end

offset_centres = [diff(frame_centres*pixels_per_mm); 0 0];
offset_centres(:,2) = -offset_centres(:,2);
[segment_transforms] = ...
        register_tiles_features(segment_frames, ...
                            'theta_range', [0], ...
                            'offset_lim', 200, ...
                            'offset_centres', offset_centres,...
                            'debug', false);
[segment_mosaic, segment_mosaic_weights, segment_transforms] = ...
    create_mosaic(segment_frames, segment_transforms);

segment_mask = segment_mosaic_weights > 0.1;
gmin = min(segment_mosaic(segment_mask));
gmax = max(segment_mosaic(segment_mask));
segment_mosaic(~segment_mask) = gmax+1;
diff_img = mean(diff_imgs,3);

figure; imgray(segment_mosaic);

save([root_dir session_dir sequence_dir 'segments\segment_data.mat'],...
    'segment_frames', 'segment_mosaics', 'segment_mosaic', 'segment_mosaic_weights', 'segment_transforms', 'diff_img', 'segment_mask');

%%
for i_seg = 1:num_segs
    figure; imgray(segment_frames(:,:,i_seg));
end
%%
%focused_idx = cell(num_segs, 1);
for i_seg = 1:num_segs

    frame_idx = segments_s{i_seg};
    num_frames = length(frame_idx);

    seg_x = motor_x(frame_idx);
    seg_y = motor_y(frame_idx);
    seg_z = motor_z(frame_idx);
    seg_s = sharpness(frame_idx);

    seg_s30 = conv(seg_s, ones(1,30)/30, 'valid');
    
    limit_s = 0.95*max(seg_s30(seg_s30 < 20));
    valid_s = [seg_s30 > limit_s & seg_s30 < 20; false];
    pt0 = 0;
    max_length = 0;
    max_pt = 0;
    pt_length = 0;
    prev_pt = 0;
    for i_pt = 1:length(valid_s)
        if valid_s(i_pt)
            if prev_pt
                pt_length = pt_length+1;
            else
                pt_length = 1;
                pt0 = i_pt;
            end
        else
            if prev_pt && pt_length > max_length
                max_length = pt_length;
                max_pt = pt0;
            end
        end
        prev_pt = valid_s(i_pt);
    end
    valid_s(end) = [];
    
    focused_idx{i_seg} = max_pt:(max_pt+max_length+28);

    figure; 
    hold all;
    plot(1:num_frames, seg_s', 'r');
    plot(1:num_frames-29, seg_s30', 'b');
    plot([1 num_frames], [limit_s limit_s], 'g');
    plot(focused_idx{i_seg}, seg_s(focused_idx{i_seg}), 'gx');
end
%%
for i_seg = 1:num_segs

    frame_idx = segments_s{i_seg};
    num_frames = length(frame_idx);

    seg_x = motor_x(frame_idx);
    seg_y = motor_y(frame_idx);
    seg_z = motor_z(frame_idx);
    seg_s = sharpness(frame_idx);

    seg_s30 = conv(seg_s, ones(1,30)/30, 'valid');
    
    limit_s = 0.95*seg_s30(end);
    not_focused = [seg_s30 < limit_s & seg_s30 > 20; false];
    
    focus_start = find(not_focused(end:-1:1),1,'first');
    if isempty(focus_start)
        focus_start = 1;
    else
        focus_start = length(seg_s30) - focus_start;
    end
    
    focused_idx{i_seg} = focus_start:length(seg_s);

    figure; 
    hold all;
    plot(1:num_frames, seg_s', 'r');
    plot(1:num_frames-29, seg_s30', 'b');
    plot([1 num_frames], [limit_s limit_s], 'g');
    plot(focused_idx{i_seg}, seg_s(focused_idx{i_seg}), 'gx');
end
%%
for i_seg = 8:num_segs

    target_mosaic = segment_frames(:,:,i_seg);
    %
    frame_idx = segments_s{i_seg};
    focused_idx_i = frame_idx(focused_idx{i_seg});
    num_frames = length(focused_idx_i);
    frames = zeros(480,640,num_frames);
    for i_im = 1:num_frames
        frames(:,:,i_im) = ...
            double(rot90(imread([im_folder 'frame' zerostr(focused_idx_i(i_im), 5) '.bmp']),2)) - ...
            diff_image;
    end
    %
    [focused_transforms, focused_counts] = ...
        register_tiles_features(frames, ...
            'ref_type', 'mosaic',...
            'sigma', 6,...
            'mosaic', target_mosaic,...
            'compound_transforms', [],...
            'theta_range', [0], ...
            'offset_lim', 120, ...
            'debug', false);

    [focused_mosaic] = ...
        create_mosaic(frames, focused_transforms);
    figure; imgray(focused_mosaic);
    %
    reg_folder = [im_folder 'registered' zerostr(i_seg,2) '\'];
    g_lims = [min(focused_mosaic(:)) max(focused_mosaic(:))];
    [final_diff_img] = write_trans_tiles(frames, focused_transforms, ...
                             reg_folder, 'frame', g_lims, focused_mosaic);

    cmd = ['ffmpeg -y -r 120 -i "' reg_folder 'frame%04d.png" -c:v libx264 -preset slow -crf 18 -an "' reg_folder 'movie.mp4"'];
    system(cmd);
    %delete([reg_folder 'frame*.png']);
end
%%
rf_d = u_load('C:\isbe\nailfold\models\vessel\detection\rf_classification\296655\predictor.mat');
rf_d.tree_root = 'C:/isbe/nailfold/models/vessel/detection/rf_classification/';
args_d = u_load('C:\isbe\nailfold\models\vessel\detection\rf_classification\296655\job_args.mat');

pred_d = predict_image(... % non-strict mode
    'image_in', segment_mosaic,...
    'decomposition_args', args_d.decomposition_args,...
    'predictor', rf_d, ...
    'prediction_type', 'rf_classification',...
    'output_type', 'detection',...
    'use_probs', 0,...
    'mask', segment_mask,...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 128,...
    'incremental_results', 0);
%%
rfs = cell(3,1);
rf{1} = u_load('C:\isbe\nailfold\models\vessel\detection\rf_classification\296655\predictor.mat');
rf{2} = u_load('C:\isbe\nailfold\models\vessel\orientation\rf_regression\296621\predictor.mat');
rf{3} = u_load('C:\isbe\nailfold\models\vessel\width\rf_regression\297037\predictor.mat');
rf_args = u_load('C:\isbe\nailfold\models\vessel\detection\rf_classification\296655\job_args.mat');

pred_dow = predict_image(... % non-strict mode
    'image_in', segment_mosaic,...
    'decomposition_args', rf_args.decomposition_args,...
    'predictor', rfs, ...
    'prediction_type', {'rf_classification', 'rf_regression', 'rf_regression'},...
    'output_type', {'detection', 'orientation', 'width'},...
    'use_probs', 0,...
    'mask', segment_mask,...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 128,...
    'incremental_results', 0);
%%
for i_f = 301:310
    [line_strength, line_orientation] = gaussian_1st_derivative_gradient(frames(:,:,i_f), 6);
    line_nms = mb_non_maximal_supp(line_strength, line_orientation);
    [line_mask] = hysterisis(line_nms);
    [y x] = find(line_mask);
    
    figure; imgray(frames(:,:,i_f));
    plot(x, y, 'm.');
end
%%
[line_strength, line_orientation] = gaussian_1st_derivative_gradient(target_mosaic, 6);
    line_nms = mb_non_maximal_supp(line_strength, line_orientation);
    [line_mask] = hysterisis(line_nms);
    [y x] = find(line_mask);
    
    figure; imgray(target_mosaic);
    plot(x, y, 'm.');
%%
max_r = 480;
max_c = 640;
for i_seg = 1:num_segs
    [r c] = size(segment_mosaics{i_seg});
    if r > max_r
        max_r = r;
    end
    if r > max_c
        max_c = c;
    end
end
segment_frames_max = nan(max_r, max_c, num_segs);
for i_seg = 1:num_segs
    [r c] = size(segment_mosaics{i_seg});
    segment_frames_max(1:r, 1:c, i_seg) = segment_mosaics{i_seg};
end
%%    
[segment_transforms_max] = ...
    register_tiles_features(segment_frames_max, ...
        'theta_range', [0], ...
        'offset_lim', 200, ...
        'offset_centres', offset_centres,...
        'debug', false);
[segment_mosaic_max] = ...
    create_mosaic(segment_frames_max, segment_transforms_max, 'rect', ~isnan(segment_frames_max));
figure; imgray(segment_mosaic_max);
%%
[segment_mosaic_c] = ...
    create_mosaic(segment_frames, segment_transforms, 'circle');
figure; imgray(segment_mosaic_c);
%%
frame_means = zeros(size(segment_frames));
for i_seg = 1:14
    frame_dir = ['N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\002wellcome\2015_02_27\L4_11_00_04\registered' zerostr(i_seg,2) '\'];
    frame_list = dir([frame_dir 'frame*.png']);

    for i_f = 1:length(frame_list)
        frame = double(imread([frame_dir frame_list(i_f).name]));
        if i_f == 1
            frame_sum = frame;
        else
            frame_sum = frame_sum + frame;
        end
    end
    frame_sum = frame_sum / length(frame_list);
    [r, c] = size(frame_sum);
    cx = round(c/2);
    cy = round(r/2);
    frame_means(:,:,i_seg) = frame_sum(cy+mosaic_lims_y, cx+mosaic_lims_x);

    figure; 
    subplot(1,2,1); imgray(segment_frames(:,:,i_seg));
    subplot(1,2,2); imgray(frame_means(:,:,i_seg));
end
%%
[segment_transforms_max] = ...
    register_tiles_features(frame_means, ...
        'theta_range', [0], ...
        'offset_lim', 200, ...
        'offset_centres', offset_centres,...
        'debug', false);
[segment_mosaic_max] = ...
    create_mosaic(frame_means, segment_transforms_max, 'rect');
%%
figure; 
a1 = subplot(2,1,1); imgray(segment_mosaic);
a2 = subplot(2,1,2); imgray(segment_mosaic_max);
linkaxes([a1 a2]);
%%
figure; 
a1 = subplot(2,1,1); imgray(final_mosaic); caxis([100 150]);
a2 = subplot(2,1,2); imgray(segment_mosaic_max);
linkaxes([a1 a2]);
%%
[~, weights_rect] = ...
    create_mosaic(frame_means, segment_transforms_max, 'rect');
[~, weights_circ] = ...
    create_mosaic(frame_means, segment_transforms_max, 'circle');
[sm, weights_diam] = ...
    create_mosaic(segment_frames, segment_transforms, 'diamond');

figure; 
a1 = subplot(2,1,1); imgray(weights_rect);
a2 = subplot(2,1,2); imgray(weights_circ);
linkaxes([a1 a2]);
%%
figure; 
a1 = subplot(2,1,1); imgray(segment_mosaic);
a2 = subplot(2,1,2); imgray(sm);
linkaxes([a1 a2]);
%%
register_sequence(... % the user's input
    'sequence',                 [],...
    'sequence_data_path',       'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\002wellcome\2015_02_27\L4_11_00_04\sequence_frames_data.dat',...
    'sequence_dir',             'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\002wellcome\2015_02_27\L4_11_00_04\',...
    'pixels_per_mm',            1000,...
    'min_stationary_frames',    120,...
    'mosaic_lims_x',            -310:309,...
    'mosaic_lims_y',            -230:229,...
    'frame_h',                  480,...
    'frame_w',                  640,...
    'sigma',                    6,...
    'intra_segment_offset',     120,...
    'inter_segment_offset',     240,...
    'plot',                     1,...
    'debug',                     1);
%%
register_sequence(... % the user's input
    'sequence',                 [],...
    'sequence_data_path',       'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\002wellcome\2015_02_27\R4_10_32_15\sequences.mat',...
    'sequence_dir',             'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\002wellcome\2015_02_27\R4_10_32_15\',...
    'dirt_image',               s.dirt_image,...
    'pixels_per_mm',            1000,...
    'min_stationary_frames',    120,...
    'mosaic_lims_x',            -310:309,...
    'mosaic_lims_y',            -230:229,...
    'frame_h',                  480,...
    'frame_w',                  640,...
    'sigma',                    6,...
    'intra_segment_offset',     120,...
    'inter_segment_offset',     240,...
    'plot',                     1,...
    'debug',                    0);
%%
root_dir = 'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\';

session_dir = [root_dir '002wellcome\2015_02_27\'];
sequences_s = dir(session_dir);

num_seqs = length(sequences_s) - 2;
sequence_dir = cell(num_seqs,1);
for i_s = 1:num_seqs
    sequence_dir{i_s} = [session_dir sequences_s(i_s+2).name '\'];
end
for i_s = 1:10
    register_sequence(... % the user's input
        'sequence',                 [],...
        'sequence_data_path',       [sequence_dir{i_s} 'sequences.mat'],...
        'sequence_dir',             sequence_dir{i_s},...
        'dirt_image',               [],...
        'pixels_per_mm',            1000,...
        'min_stationary_frames',    120,...
        'mosaic_lims_x',            -310:309,...
        'mosaic_lims_y',            -230:229,...
        'frame_h',                  480,...
        'frame_w',                  640,...
        'sigma',                    6,...
        'intra_segment_offset',     120,...
        'inter_segment_offset',     240,...
        'plot',                     1,...
        'debug',                    0);
end
%%
%%
root_dir = 'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\';

session_dir = [root_dir '051wellcome\2015_03_30\'];
sequences_s = dir(session_dir);

num_seqs = length(sequences_s) - 2;
sequence_dir = cell(num_seqs,1);
for i_s = 1:num_seqs
    sequence_dir{i_s} = [session_dir sequences_s(i_s+2).name '\'];
end
for i_s = 1:num_seqs
    register_sequence(... % the user's input
        'sequence',                 [],...
        'sequence_data_path',       [sequence_dir{i_s} 'sequence.mat'],...
        'sequence_dir',             sequence_dir{i_s},...
        'dirt_image',               [],...
        'pixels_per_mm',            1000,...
        'min_stationary_frames',    120,...
        'mosaic_lims_x',            -310:309,...
        'mosaic_lims_y',            -230:229,...
        'frame_h',                  480,...
        'frame_w',                  640,...
        'sigma',                    6,...
        'intra_segment_offset',     120,...
        'inter_segment_offset',     240,...
        'plot',                     1,...
        'debug',                    0);
end
%%
i_s = 3;
s = load([sequence_dir{i_s} 'segments\segment_data.mat']);
register_sequence(... % the user's input
        'sequence',                 [],...
        'sequence_data_path',       [sequence_dir{i_s} 'sequences.mat'],...
        'sequence_dir',             sequence_dir{i_s},...
        'dirt_image',               s.dirt_image,...
        'pixels_per_mm',            1000,...
        'min_stationary_frames',    120,...
        'mosaic_lims_x',            -310:309,...
        'mosaic_lims_y',            -230:229,...
        'frame_h',                  480,...
        'frame_w',                  640,...
        'sigma',                    6,...
        'intra_segment_offset',     120,...
        'inter_segment_offset',     240,...
        'plot',                     1,...
        'debug',                    0);