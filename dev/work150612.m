%s = load('C:\isbe\nailfold\test_analysis.mat');
i_seg = 14;
segment_start_x = segment_transforms(1,3,i_seg);
segment_start_y = segment_transforms(2,3,i_seg);

segment_apexes_x = 2*s.candidate_xy(s.selected_distal,1) - segment_start_x;
segment_apexes_y = 2*s.candidate_xy(s.selected_distal,2) - segment_start_y;

figure; imgray(segment_frames(:,:,i_seg));
plot(segment_apexes_x, segment_apexes_y, 'rx');



in_segment = ...
    segment_apexes_x >= 1 & segment_apexes_x <= 640 &...
    segment_apexes_y >= 1 & segment_apexes_y <= 480;

bounding_boxes = 2*s.apex_measures.distal.bounding_box(:,:,in_segment);
bounding_boxes(:,1,:) = bounding_boxes(:,1,:) - segment_start_x;
bounding_boxes(:,2,:) = bounding_boxes(:,2,:) - segment_start_y;

figure; imgray(segment_mosaic);
plot(2*s.candidate_xy(s.selected_distal,1), 2*s.candidate_xy(s.selected_distal,2), 'rx');
plot(...
    [segment_start_x segment_start_x+640 segment_start_x+640 segment_start_x segment_start_x],...
    [segment_start_y segment_start_y segment_start_y+480 segment_start_y+480 segment_start_y]);


%plot(s.candidate_xy(in_segment,1), s.candidate_xy(in_segment,2), 'go');
%%
figure; imgray(segment_frames(:,:,i_seg));
for i_b = 1:sum(in_segment)
    tl_x = bounding_boxes(1,1,i_b);
    tl_y = bounding_boxes(1,2,i_b);
    
    bl_x = bounding_boxes(2,1,i_b);
    bl_y = bounding_boxes(2,2,i_b);
    
    tr_x = bounding_boxes(3,1,i_b);
    tr_y = bounding_boxes(3,2,i_b);  
    
    br_x = bounding_boxes(4,1,i_b);
    br_y = bounding_boxes(4,2,i_b);
    
    plot([tl_x tr_x br_x bl_x tl_x], [tl_y tr_y br_y bl_y tl_y]);
    plot(tl_x, tl_y, 'x');
    plot(tr_x, tr_y, '+');
    plot(bl_x, bl_y, 's');
    plot(br_x, br_y, 'o');
end
%%
figure; imgray(segment_frames(:,:,i_seg));
for i_b = 1:sum(in_segment)
    tl_x = bounding_boxes(1,1,i_b);
    tl_y = bounding_boxes(1,2,i_b);
    
    bl_x = bounding_boxes(2,1,i_b);
    bl_y = bounding_boxes(2,2,i_b);
    
    tr_x = bounding_boxes(3,1,i_b);
    tr_y = bounding_boxes(3,2,i_b);  
    
    br_x = bounding_boxes(4,1,i_b);
    br_y = bounding_boxes(4,2,i_b);
    
    x_lim = round(abs(tl_x - tr_x)) + 1;
    y_lim = round(abs(bl_y - tl_y)) + 1;

    ux = repmat(linspace(0, tr_x - tl_x, x_lim), y_lim, 1);
    uy = repmat(linspace(0, tr_y - tl_y, x_lim), y_lim, 1);
    vx = repmat(linspace(0, bl_x - tl_x, y_lim)', 1, x_lim);
    vy = repmat(linspace(0, bl_y - tl_y, y_lim)', 1, x_lim);

    x = tl_x + ux + vx;
    y = tl_y + uy + vy;
    
    plot(x, y, '.');
end
%%
[s_rows s_cols] = size(segment_frames(:,:,1));
sx = repmat((segment_start_x+(0:s_cols-1))/2, s_rows, 1);
sy = repmat((segment_start_y+(0:s_rows-1))'/2, 1, s_cols);
mosaic_vessels = interp2(s.vessel_predictions(:,:,1), sx, sy, '*linear', 0);

offset_yx = (size(segment_mosaics{i_seg}) - size(segment_frames(:,:,1)))/2;

frame_dir = 'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\002wellcome\2015_02_27\L4_11_00_04\registered07\';
frame_list = dir([frame_dir 'frame*.png']);
num_frames = length(frame_list);
%%
for i_b = 2:sum(in_segment)
    tl_x = bounding_boxes(1,1,i_b);
    tl_y = bounding_boxes(1,2,i_b);

    bl_x = bounding_boxes(2,1,i_b);
    bl_y = bounding_boxes(2,2,i_b);

    tr_x = bounding_boxes(3,1,i_b);
    tr_y = bounding_boxes(3,2,i_b);  

    br_x = bounding_boxes(4,1,i_b);
    br_y = bounding_boxes(4,2,i_b);

    x_lim = round(abs(tl_x - tr_x)) + 1;
    y_lim = round(abs(bl_y - tl_y)) + 1;

    if rem(x_lim,2)
        x_lim = x_lim + 1;
        tr_x = tr_x + 1;
        br_x = br_x + 1;
    end
    if rem(y_lim,2)
        y_lim = y_lim + 1;
        bl_y = bl_y + 1;
        br_y = br_y + 1;
    end

    ux = repmat(linspace(0, tr_x - tl_x, x_lim), y_lim, 1);
    uy = repmat(linspace(0, tr_y - tl_y, x_lim), y_lim, 1);
    vx = repmat(linspace(0, bl_x - tl_x, y_lim)', 1, x_lim);
    vy = repmat(linspace(0, bl_y - tl_y, y_lim)', 1, x_lim);

    x = tl_x + ux + vx;
    y = tl_y + uy + vy;    
    vessel_bg_mask = interp2(mosaic_vessels, x, y, '*linear', 0) < 0.8;

   
    
    x = x + offset_yx(2);
    y = y + offset_yx(1);
    
    vessel_dir = [frame_dir 'vessel' zerostr(i_b,2) '\'];
    create_folder(vessel_dir);
    vessel_patches = zeros(size(x,1), size(x,2), num_frames);
    
    for i_f = 1:num_frames
        frame = imread([frame_dir frame_list(i_f).name]);
        vessel_patch = interp2(double(frame), x, y, '*linear', 0);
        vessel_patches(:,:,i_f) = vessel_patch;
    end
    vessel_mosaic = mean(vessel_patches,3);
    vessel_bg = vessel_mosaic(vessel_bg_mask);
    
    figure; 
    subplot(1,2,1); imgray(vessel_bg_mask);
    subplot(1,2,2); imgray(vessel_mosaic);
    
    for i_f = 1:num_frames
        vessel_patch = vessel_patches(:,:,i_f);
        vessel_patch(vessel_bg_mask) = vessel_bg;
        imwrite(uint8(vessel_patch), [vessel_dir frame_list(i_f).name]);
    end   

    cmd = ['ffmpeg -y -r 120 -i "' vessel_dir 'frame%04d.png" -c:v libx264 -preset slow -crf 18 -an "' vessel_dir 'movie.mp4"'];
    system(cmd);

    delete([vessel_dir 'frame*.png']);
end
%%
root_dir = 'N:\nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\';
session_dir = '002wellcome\2015_02_27\';
sequence_dir = 'L4_11_00_04\';
im_folder = [root_dir session_dir sequence_dir];

sequence = u_load([im_folder 'sequences.mat']);

pixels_per_mm = 1000;

[segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness] = ...
    get_stationary_segments(sequence);

load([root_dir session_dir sequence_dir 'segments\segment_data.mat']);
%%
num_segs = length(segments_ns);
frame_sum = zeros(480,640);
total_frames = 0;
for i_seg = 1:num_segs

    frame_idx = segments_ns{i_seg};
    num_frames = length(frame_idx);
    total_frames = total_frames + num_frames;

    seg_x = motor_x(frame_idx)*pixels_per_mm;
    seg_y = motor_y(frame_idx)*pixels_per_mm; 

    x1 = segment_transforms(1,3,1);
    y1 = segment_transforms(2,3,1);

    offset_centres = [x1 y1; diff(seg_x) -diff(seg_y)];

    figure; imgray(segment_mosaic);
    plot(x1 + [0 640 640 0 0], y1+[0 0 480 480 0]);
    plot(cumsum(offset_centres(:,1)), cumsum(offset_centres(:,2)));

    for i_im = 1:num_frames
        frame_sum = frame_sum + double(rot90(imread([im_folder 'frame' zerostr(frame_idx(i_im), 5) '.bmp']),2));
    end
end
%%
figure; imgray(segment_mosaic);
initial_transforms = repmat(eye(3), 1, 1, num_frames);
initial_transforms(1,3,1) = x1;
initial_transforms(2,3,1) = y1;
    
for i_f = 2:num_frames
    initial_transforms(1,3,i_f) = initial_transforms(1,3,i_f-1) + seg_x(i_f) - seg_x(i_f-1);
    initial_transforms(2,3,i_f) = initial_transforms(2,3,i_f-1) - seg_y(i_f) + seg_y(i_f-1);
    
    plot(initial_transforms(1,3,i_f) + 320, initial_transforms(2,3,i_f) + 240, 'x');
end
%%
[moving_transforms, moving_counts] = ...
        register_tiles_features(frames, ...
            'ref_type', 'mosaic',...
            'mosaic', segment_mosaic,...
            'sigma', 6,...
            'compound_transforms', initial_transforms,...
            'offset_centres', offset_centres,...
            'theta_range', [0], ...
            'offset_lim', 120, ...
            'debug', false);
        
[moving_mosaic,~,moving_transforms] = ...
        create_mosaic(frames, moving_transforms);

[moving_diff_img] = write_trans_tiles(frames, moving_transforms, ...
                             nan, [], [], moving_mosaic);
                         
diff_image = frame_sum / total_frames;  
%%
s = load(segment_data_path);
v = load(vessel_data_path);
%%
for i_seg = 1:7
    frame_dir = ['N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\002wellcome\2015_02_27\L4_11_00_04\registered' zerostr(i_seg,2) '\'];
    segment_vessel_videos(i_seg, s, v, frame_dir);
end

    