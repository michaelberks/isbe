root_dir = nailfoldroot(1);

session_dir = [root_dir 'camera_capture\testing_new_camera\bar\2015_01_16\'];
sequences_s = dir(session_dir);

num_seqs = length(sequences_s) - 2;
sequence_dir = cell(num_seqs,1);
for i_s = 1:num_seqs
    sequence_dir{i_s} = sequences_s(i_s+2).name;
end
%%
num_seqs = length(sequence_dir);
%sequences = cell(num_seqs,1);
for i_s = [15 17]%1:num_seqs
    [sequences{i_s}] = read_sequence_from([session_dir sequence_dir{i_s} '\sequence_properties.txt']);
end
%%

for i_s = [9 10 12 15 17]
    
    [segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness] = ...
        get_stationary_segments(sequences{i_s});
    
    xy_speed = [0; sqrt(diff(motor_x).^2 + diff(motor_y).^2)];
    z_speed = [0; diff(motor_z)];   
    frame_nums = 1:length(xy_speed);
    
    figure;
    a1 = subplot(3,1,1);hold on; 
    plot(xy_speed, 'b');
    plot(abs(z_speed), 'r');
    xlabel('Time');
    ylabel('Motor speed');
    legend({'|XY| speed', '|Z| speed'});
    set(a1, 'ylim', [ 0 0.04]);
    
    %plot(sharpness_change / (max(sharpness_change)-min(sharpness_change)), 'g--');
    a2 = subplot(3,1,2); hold on;
    xlabel('Time');
    ylabel('Z position');
    plot(0,0,'r-');
    plot(0,0,'b-');
    legend(a2, {'Stationary/moving in Z', '|XY| movement'});
    
    a3 = subplot(3,1,3); hold on;
    xlabel('Time');
    ylabel('Sharpness');
    plot(0,0,'r-');
    plot(0,0,'b-');
    legend(a3, {'Stationary/moving in Z', 'XY movement'});
    
    for i_sg = 1:length(segments_s)
        plot(a2, frame_nums(segments_s{i_sg}), motor_z(segments_s{i_sg}), 'r.');
        plot(a3, frame_nums(segments_s{i_sg}), sharpness(segments_s{i_sg}), 'r.');
    end
    
    for i_sg = 1:length(segments_ns)
        plot(a2, frame_nums(segments_ns{i_sg}), motor_z(segments_ns{i_sg}), 'b.');
        plot(a3, frame_nums(segments_ns{i_sg}), sharpness(segments_ns{i_sg}), 'b.');
    end
      
end
%%
for i_s = [9 10 12 15 17]
    sequence = sequences{i_s};
    [segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness] = ...
        get_stationary_segments(sequence);

    num_segs = length(segments_s);
    discard = false(num_segs,1);
    for i_seg = 1:num_segs    
        discard(i_seg) = length(segments_s{i_seg}) < 120;
    end
    segments_s(discard) = [];
    num_segs = length(segments_s);   
    im_folder = [root_dir 'camera_capture\testing_new_camera\bar\2015_01_16\' sequence.sequence_name '\'];

    %    
    hw = 320 / 1245;
    hh = 240 / 1245;
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

        if 0
            [best_s, best_frames] = sort(seg_s, 'descend');
            idx = best_frames(round(linspace(1, num_frames, 12)));

            figure;
            for i_f = 1:11

                frame = imread([im_folder 'frame' zerostr(frame_idx(idx(i_f)), 5) '.bmp']);
                subplot(3,4,i_f);
                imgray(frame);
                title(['F_i = ' num2str(round(num_frames *(i_f-1)/10)) ', Sharpness = ' num2str(seg_s(idx(i_f)))]);
            end

            subplot(3,4,12); hist(seg_s, 10);
        end
    end
end
%%
mosaic_lims_x = -310:309;
mosaic_lims_y = -230:229;

for i_s = [9 10 12 15 17]
    
    sequence = sequences{i_s};
    [segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness] = ...
        get_stationary_segments(sequence);

    num_segs = length(segments_s);
    discard = false(num_segs,1);
    for i_seg = 1:num_segs    
        discard(i_seg) = length(segments_s{i_seg}) < 120;
    end
    segments_s(discard) = [];
    num_segs = length(segments_s);   
    im_folder = [root_dir 'camera_capture\testing_new_camera\bar\2015_01_16\' sequence.sequence_name '\'];
    
    nailfold_mosaic = cell(num_segs,1);
    frame_centres = zeros(num_segs, 2);
    matched_counts = cell(num_segs,1);
    frames = zeros(length(mosaic_lims_y), length(mosaic_lims_x), num_segs);
    %
    for i_seg = 1:num_segs

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
            frames(:,:,i_im) = imread([im_folder 'frame' zerostr(frame_idx(best_frames(i_im)), 5) '.bmp']);
        end

        [compound_transforms, matched_counts{i_seg}] = ...
            register_tiles_features(frames, ...
                                'theta_range', [0], ...
                                'offset_lim', 40, ...
                                'debug', false);

        [nailfold_mosaic{i_seg}] = ...
            rot90(create_mosaic(frames, compound_transforms),2);
        
        

        [r, c] = size(nailfold_mosaic{i_seg});
        cx = round(c/2);
        cy = round(r/2);
        frames(:,:,i_seg) = nailfold_mosaic{i_seg}(cy+mosaic_lims_y, cx+mosaic_lims_x);
        
        %figure; imgray(nailfold_mosaic{i_seg});
        %figure; imgray(frames(:,:,i_seg));
        %figure; hist( matched_counts{i_seg}(2:end));
    end
    
    offset_centres = [diff(frame_centres*1245); 0 0];
    offset_centres(:,2) = -offset_centres(:,2);
    [compound_transforms] = ...
            register_tiles_features(frames, ...
                                'theta_range', [0], ...
                                'offset_lim', 20, ...
                                'offset_centres', offset_centres,...
                                'debug', false);
    [final_mosaic, mosaic_weights] = ...
        create_mosaic(frames, compound_transforms);
    
    mask = mosaic_weights > 0.1;
    gmin = min(final_mosaic(mask));
    gmax = max(final_mosaic(mask));
    final_mosaic(~mask) = gmax+1;
    
    figure; imgray(final_mosaic);
    write_im_from_colormap(final_mosaic, [im_folder 'static_mosaic.png'], gray(256), [gmin gmax]);
    
end
%%

figure; hold all; axis equal; a1 = gca;
for i_seg = 1:num_segs
  
    frame_idx = segments_s{i_seg};

    seg_x = motor_x(frame_idx);
    seg_y = motor_y(frame_idx);
    
    cx = median(seg_x);
    cy = median(seg_y);

    plot(a1, cx, cy, 'rx');
    plot(a1, cx+frame_rect_x, cy+frame_rect_y, 'b');  
end
for i_seg = 1
    frame_idx = segments_ns{i_seg};

    seg_x = motor_x(frame_idx);
    seg_y = motor_y(frame_idx);
    for i_f = 1:length(seg_x)
        
        cx = seg_x(i_f);
        cy = seg_y(i_f);
        plot(a1, cx+frame_rect_x, cy+frame_rect_y, 'c--');
    end
    
end

%%
num_nsegs = length(segments_ns);
frame_idx = [];
for i_seg = 1:num_nsegs 
    num_frames = length(segments_ns{i_seg});
    frame_idx_i = segments_ns{i_seg}(1:5:num_frames);
    frame_idx = [frame_idx; frame_idx_i]; %#ok
end
%%    
num_frames = length(frame_idx);
%
frames = zeros(480,640,num_frames);
motion_transforms = repmat(eye(3), 1, 1, num_frames);
for i_f = 1:num_frames
    frames(:,:,i_f) = rot90(imread([im_folder 'frame' zerostr(frame_idx(i_f), 5) '.bmp']),2);
    motion_transforms(1,3,i_f) = round((motor_x(frame_idx(i_f)) - motor_x(1))*1100);
    motion_transforms(2,3,i_f) = -round((motor_y(frame_idx(i_f)) - motor_y(1))*1100);
end

[motion_transforms_new, motion_counts] = ...
    register_tiles_features(frames, ...
        'ref_type', 'mosaic',...
        'mosaic', final_mosaic,...
        'compound_transforms', motion_transforms,...
        'theta_range', [0], ...
        'offset_lim', 120, ...
        'debug', false);

%%
[motion_mosaic, motion_weights] = ...
    create_mosaic(frames, motion_transforms_new, [], [], final_mosaic, mosaic_weights);

mask = motion_weights > 0.1;
gmin = min(motion_mosaic(mask));
gmax = max(motion_mosaic(mask));
motion_mosaic(~mask) = gmax+1;
    
figure; imgray(motion_mosaic);




