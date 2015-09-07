root_dir = nailfoldroot(1);   
[segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness] = ...
    get_stationary_segments(sequence);

xy_speed = [0; sqrt(diff(motor_x).^2 + diff(motor_y).^2)];
z_speed = [0; diff(motor_z)];   
frame_nums = 1:length(xy_speed);

pixels_per_mm = 1000;

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
      
%%
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
    
[segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness] = ...
    get_stationary_segments(sequence);

num_segs = length(segments_s);
discard = false(num_segs,1);
for i_seg = 1:num_segs    
    discard(i_seg) = length(segments_s{i_seg}) < 120;
end
segments_s(discard) = [];
num_segs = length(segments_s);   
im_folder = [root_dir sequence.image_dir(18:end) '/'];

nailfold_mosaic = cell(num_segs,1);
frame_centres = zeros(num_segs, 2);
matched_counts = cell(num_segs,1);
segment_frames = zeros(length(mosaic_lims_y), length(mosaic_lims_x), num_segs);
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
                            'offset_lim', 120, ...
                            'debug', false);

    [nailfold_mosaic{i_seg}] = ...
        rot90(create_mosaic(frames, compound_transforms),2);


    [r, c] = size(nailfold_mosaic{i_seg});
    cx = round(c/2);
    cy = round(r/2);
    segment_frames(:,:,i_seg) = nailfold_mosaic{i_seg}(cy+mosaic_lims_y, cx+mosaic_lims_x);

    %figure; imgray(nailfold_mosaic{i_seg});
    %figure; imgray(frames(:,:,i_seg));
    %figure; hist( matched_counts{i_seg}(2:end));
end
%%
offset_centres = [diff(frame_centres*pixels_per_mm); 0 0];
offset_centres(:,2) = -offset_centres(:,2);
[compound_transforms] = ...
        register_tiles_features(segment_frames, ...
                            'theta_range', [0], ...
                            'offset_lim', 200, ...
                            'offset_centres', offset_centres,...
                            'debug', false);
[final_mosaic, mosaic_weights] = ...
    create_mosaic(segment_frames, compound_transforms);

mask = mosaic_weights > 0.1;
gmin = min(final_mosaic(mask));
gmax = max(final_mosaic(mask));
final_mosaic(~mask) = gmax+1;

figure; imgray(final_mosaic);