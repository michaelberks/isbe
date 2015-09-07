function segment_vessel_videos(i_seg, s, v, frame_dir)

segment_start_x = s.segment_transforms(1,3,i_seg);
segment_start_y = s.segment_transforms(2,3,i_seg);

segment_apexes_x = 2*v.candidate_xy(v.selected_distal,1) - segment_start_x;
segment_apexes_y = 2*v.candidate_xy(v.selected_distal,2) - segment_start_y;

in_segment = ...
    segment_apexes_x >= 1 & segment_apexes_x <= 640 &...
    segment_apexes_y >= 1 & segment_apexes_y <= 480;

bounding_boxes = 2*v.apex_measures.distal.bounding_box(:,:,in_segment);
bounding_boxes(:,1,:) = bounding_boxes(:,1,:) - segment_start_x;
bounding_boxes(:,2,:) = bounding_boxes(:,2,:) - segment_start_y;

figure; imgray(s.segment_mosaic);
plot(2*v.candidate_xy(v.selected_distal,1), 2*v.candidate_xy(v.selected_distal,2), 'rx');
plot(...
    [segment_start_x segment_start_x+640 segment_start_x+640 segment_start_x segment_start_x],...
    [segment_start_y segment_start_y segment_start_y+480 segment_start_y+480 segment_start_y]);

%%
[s_rows s_cols] = size(s.segment_frames(:,:,1));
sx = repmat((segment_start_x+(0:s_cols-1))/2, s_rows, 1);
sy = repmat((segment_start_y+(0:s_rows-1))'/2, 1, s_cols);
mosaic_vessels = interp2(v.vessel_predictions(:,:,1), sx, sy, '*linear', 0);

offset_yx = (size(s.segment_mosaics{i_seg}) - size(s.segment_frames(:,:,1)))/2;
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
        imwrite(uint8(vessel_patch), [vessel_dir frame_list(i_f).name]);
    end
    cmd = ['ffmpeg -y -r 120 -i "' vessel_dir 'frame%04d.png" -c:v libx264 -preset slow -crf 18 -an "' vessel_dir 'movie.mp4"'];
    system(cmd);
    delete([vessel_dir 'frame*.png']);
    
    vessel_mosaic = mean(vessel_patches,3);
    vessel_bg = vessel_mosaic(vessel_bg_mask);
    
%     figure; 
%     subplot(1,2,1); imgray(vessel_bg_mask);
%     subplot(1,2,2); imgray(vessel_mosaic);
    
    for i_f = 1:num_frames
        vessel_patch = vessel_patches(:,:,i_f);
        vessel_patch(vessel_bg_mask) = vessel_bg;
        imwrite(uint8(vessel_patch), [vessel_dir frame_list(i_f).name]);
    end   

    cmd = ['ffmpeg -y -r 120 -i "' vessel_dir 'frame%04d.png" -c:v libx264 -preset slow -crf 18 -an "' vessel_dir 'movie_masked.mp4"'];
    system(cmd);

    delete([vessel_dir 'frame*.png']);
end

    