function [] = register_sequence(varargin)
%REGISTER_SEQUENCE *Insert a one line summary here*
%   [] = register_sequence()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 15-Jun-2015
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
args = u_packargs(varargin, 0, ... % the user's input
    'sequence',                 [],...
    'sequence_data_path',       [],...
    'sequence_dir',             [],...
    'dirt_image',               [],...
    'pixels_per_mm',            1000,...
    'min_stationary_frames',    120,...
    'focus_smoothing_window',   30,... 
    'min_focus_pct',            0.95,...
    'max_focus_val',            50,...
    'bad_frame_thresh',         0.3,...
    'lower_g_lim',              0.25,...
    'upper_g_lim',              0.75,...
    'mosaic_lims_x',            -310:309,...
    'mosaic_lims_y',            -230:229,...
    'frame_h',                  480,...
    'frame_w',                  640,...
    'sigma',                    6,...
    'theta_range',              0,...
    'intra_segment_offset',     120,...
    'inter_segment_offset',     240,...
    'plot',                     0,...
    'debug',                    0);
clear varargin;

%Check if we've already got the sequence structure, if not load it in or
%compute from sequence properties (latter will be slooowwwww)
if isempty(args.sequence)
    if isempty(args.sequence_data_path)
         args.sequence = read_sequence_from([args.sequence_dir '\sequence_properties.txt']);
    else
        if strcmpi(args.sequence_data_path(end-2:end), 'mat')
            args.sequence = u_load(args.sequence_data_path);
        else
            args.sequence = load(args.sequence_data_path);
        end
        
    end
end

%--------------------------------------------------------------------------
%Compute stationary segments, discard any that were stationary for less
%than 120 frames (1 sec)
[segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness] = ...
    get_stationary_segments(args.sequence, args.min_stationary_frames);
num_segs = length(segments_s);


%--------------------------------------------------------------------------
%Compute backgorund dirt image from the average of the all the
%non-stationary frames
if isempty(args.dirt_image)
    frame_sum = zeros(args.frame_h, args.frame_w);
    total_frames = 0;
    for i_seg = 1:length(segments_ns)

        frame_idx = segments_ns{i_seg};
        num_frames = length(frame_idx);
        total_frames = total_frames + num_frames;

        for i_im = 1:num_frames
            frame = double(rot90(imread([args.sequence_dir 'frame' zerostr(frame_idx(i_im), 5) '.bmp']),2));
            if is_frame_ok(frame, args.bad_frame_thresh, args.lower_g_lim, args.upper_g_lim)
                frame_sum = frame_sum + frame;
            else
                total_frames = total_frames - 1;
            end
        end
    end
    args.dirt_image = frame_sum / total_frames;  
    
elseif ischar(args.dirt_image) && exist(args.dirt_image, 'file')
    s = load(args.dirt_image, 'dirt_image');
    if isfield(s, 'dirt_image')
        args.dirt_image = s.dirt_image;
    end
    clear s;
end
dirt_image = args.dirt_image; %#ok
create_folder([args.sequence_dir 'segments\']);

progress = 1; %#ok
save([args.sequence_dir 'segments\segment_data.mat'],...
    'dirt_image', 'progress');
clear dirt_image;
%--------------------------------------------------------------------------
%Make segment mosaics
frame_centres = zeros(num_segs, 2);
matched_counts = cell(num_segs,1);
segment_frames = zeros(length(args.mosaic_lims_y), length(args.mosaic_lims_x), num_segs);
segment_mosaics = cell(num_segs,1);
focused_idx = cell(num_segs, 1);
focused_transforms = cell(num_segs,1);
%
reject_segments = false(num_segs,1);
for i_seg = 1:num_segs
    
    frame_idx = segments_s{i_seg};
    num_frames = length(frame_idx);

    seg_x = motor_x(frame_idx);
    seg_y = motor_y(frame_idx);
    seg_z = motor_z(frame_idx);
    seg_s = sharpness(frame_idx);

    frame_centres(i_seg,1) = seg_x(round(num_frames/2));
    frame_centres(i_seg,2) = seg_y(round(num_frames/2));
    
    seg_s_cleaned = seg_s;
    seg_s_cleaned(seg_s > args.max_focus_val) = 0;
    
    [~, best_frames] = sort(seg_s_cleaned, 'descend'); %#ok

    target_frames = min(30,round(num_frames/4));

    frames = zeros(args.frame_h,args.frame_w,target_frames);
    frame_count = 0;
    i_im = 1;
    while frame_count < target_frames && i_im <= num_frames
        %frame = double(rot90(imread([args.sequence_dir 'frame' zerostr(frame_idx(best_frames(i_im)), 5) '.bmp']),2));       
        frame = double(rot90(imread([args.sequence_dir 'frame' zerostr(frame_idx(end+1-i_im), 5) '.bmp']),2));
        if is_frame_ok(frame,  args.bad_frame_thresh, args.lower_g_lim, args.upper_g_lim)
            frame_count = frame_count + 1;
            frames(:,:,frame_count) = frame  - args.dirt_image;       
        elseif args.debug
            %figure; imgray(frame);
        end
        %Increment frame index to load in next frame
        i_im = i_im + 1;
    end
    
    if frame_count == 0
        %No good frames - need to discard segment
        reject_segments(i_seg) = 1;
              
        if args.debug
            figure; imgray(frame);
            title(['Segment ' num2str(i_seg) ' failed.']);
            xlabel(args.sequence_dir);
        end
        continue;
        
    elseif frame_count < target_frames
        %Will happen if not enough ok frames to meet target
        %Discard unused zeros from frames container
        frames(:,:,frame_count+1:end) = [];
    end

    [compound_transforms, matched_counts{i_seg}] = ...
        register_tiles_features(frames, ...
                            'theta_range', args.theta_range, ...
                            'offset_lim', args.intra_segment_offset, ...
                            'sigma', args.sigma,...
                            'debug', false);

    [segment_mosaics{i_seg}] = ...
        create_mosaic(frames, compound_transforms);

    [r, c] = size(segment_mosaics{i_seg});
    cx = round(c/2);
    cy = round(r/2);
    segment_frames(:,:,i_seg) = segment_mosaics{i_seg}(cy+args.mosaic_lims_y, cx+args.mosaic_lims_x);
    
    %Also extract the largest consectutive sequence of focused frames
    %(measure relative to the highest 30 frame average sharpness score
    seg_s_smoothed = conv(seg_s, ones(1,args.focus_smoothing_window)/args.focus_smoothing_window, 'valid');   
    %focused_idx{i_seg} = get_max_run(seg_s, seg_s_smoothed, args.min_focus_pct, args.max_focus_val);
    
    dz = diff(seg_z);
    ni = find(dz(end:-1:1), 1);
    if isempty(ni)
        ni = num_frames;
    end
    focused_idx{i_seg} = num_frames - ni + (1:ni); 
    num_frames = length(focused_idx{i_seg});
    
    if args.debug
        limit_s = args.min_focus_pct*max(seg_s_smoothed);
        figure; 
        subplot(1,2,1); imgray(segment_frames(:,:,i_seg));
        subplot(1,2,2); hold all;
        plot(1:length(seg_s), seg_s', 'r');
        plot(1:length(seg_s_smoothed), seg_s_smoothed', 'b');
        plot([1 length(seg_s)], [limit_s limit_s], 'g');
        plot(focused_idx{i_seg}, seg_s(focused_idx{i_seg}), 'gx');
    end
    
    if num_frames < 120
        continue;
    end
    
    %Re-register all the focused frames
    focused_idx_i = frame_idx(focused_idx{i_seg});   
    frames = zeros(args.frame_h,args.frame_w,num_frames);
    reject_frames = false(num_frames,1);
    for i_im = 1:num_frames
        frame = ...
            double(rot90(imread([args.sequence_dir 'frame' zerostr(focused_idx_i(i_im), 5) '.bmp']),2));
     
        if 1 || is_frame_ok(frame, args.bad_frame_thresh, args.lower_g_lim, args.upper_g_lim)
            frames(:,:,i_im) = frame - args.dirt_image;
        else
            reject_frames(i_im) = 1;
        end
    end
    %Discard bad frames
    frames(:,:,reject_frames) = [];
    
    %What to do if no frames? Don't necessarily need to reject segment, but
    %it does mean we can't build a flow video
    if isempty(frames)
        continue;
    end
    
    %
    [focused_transforms{i_seg}] = ...
        register_tiles_features(frames, ...
            'ref_type', 'mosaic',...
            'sigma', args.sigma,...
            'mosaic', segment_frames(:,:,i_seg),...
            'compound_transforms', [],...
            'theta_range', args.theta_range, ...
            'offset_lim', args.intra_segment_offset, ...
            'debug', false);

    [segment_mosaics{i_seg}] = ...
        create_mosaic(frames, focused_transforms{i_seg});
    
    [r, c] = size(segment_mosaics{i_seg});
    cx = round(c/2);
    cy = round(r/2);
    segment_frames(:,:,i_seg) = segment_mosaics{i_seg}(cy+args.mosaic_lims_y, cx+args.mosaic_lims_x);

    %Finally, make a video of the segment
    reg_folder = [args.sequence_dir 'registered' zerostr(i_seg,2) '\'];
    g_lims = [min(segment_mosaics{i_seg}(:)) max(segment_mosaics{i_seg}(:))];
    delete([reg_folder 'frame*']);
    write_trans_tiles(frames, focused_transforms{i_seg}, ...
                             reg_folder, 'frame', g_lims, segment_mosaics{i_seg});

    cmd = ['ffmpeg -y -r 120 -i "' reg_folder 'frame%04d.png" -c:v libx264 -preset slow -crf 18 -an "' reg_folder 'movie.mp4"'];
    system(cmd);
    delete([reg_folder 'frame*.png']);
end

%--------------------------------------------------------------------------
%Make complete mosaic
progress = 2; %#ok
save([args.sequence_dir 'segments\segment_data.mat'],...
    'segment_frames', 'segment_mosaics', 'focused_idx', 'focused_transforms',...
    'reject_segments', 'progress',  '-append');

%First discard any rejected segments
if any(reject_segments)
    segment_frames(:,:,reject_segments) = [];
    segment_mosaics(reject_segments) = []; %#ok
    frame_centres(reject_segments,:) = [];
end

if size(frame_centres, 1) < 2
    return;
end

%Use motor positions to initialise alignment of frames
offset_centres = [diff(frame_centres*args.pixels_per_mm); 0 0];
offset_centres(:,2) = -offset_centres(:,2);

[segment_transforms] = ...
        register_tiles_features(segment_frames, ...
                            'theta_range', args.theta_range, ...
                            'offset_lim', args.inter_segment_offset, ...
                            'offset_centres', offset_centres,...
                            'debug', false);
[segment_mosaic, segment_mosaic_weights, segment_transforms] = ...
    create_mosaic(segment_frames, segment_transforms, 'diamond'); %#ok
segment_mask = segment_mosaic_weights > 0.1; %#ok

progress = 3; %#ok
save([args.sequence_dir 'segments\segment_data.mat'],...
    'segment_mosaic', 'segment_mosaic_weights',...
    'segment_transforms', 'segment_mask', 'progress', '-append');

%--------------------------------------------------------------------------
%End of main function
if args.plot
    figure; imgray(segment_mosaic);
    plot_segment_traces(args.pixels_per_mm, segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness);
end
%--------------------------------------------------------------------------

%%
%**************************************************************************
%**************************************************************************
function [run_idx] = get_max_run(seg_s, seg_s_smoothed, min_pct, max_val) %#ok

limit_s = min_pct*seg_s_smoothed(end);
not_focused = seg_s_smoothed < limit_s | seg_s_smoothed > max_val; 
focus_start = find(not_focused(end:-1:1),1,'first');
if isempty(focus_start)
    focus_start = 1;
else
    focus_start = length(seg_s) - focus_start;
end

run_idx = focus_start:length(seg_s);
    
% pt0 = 0;
% max_length = 0;
% max_pt = 0;
% pt_length = 0;
% prev_pt = 0;
% for i_pt = 1:length(input_v)
%     if input_v(i_pt)
%         if prev_pt
%             pt_length = pt_length+1;
%         else
%             pt_length = 1;
%             pt0 = i_pt;
%         end
%     else
%         if prev_pt && pt_length > max_length
%             max_length = pt_length;
%             max_pt = pt0;
%         end
%     end
%     prev_pt = input_v(i_pt);
% end
% run_idx = max_pt:(max_pt+max_length+28);

%**************************************************************************
%**************************************************************************
function [is_valid] = is_frame_ok(frame, thresh_pct, lower_g_lim, upper_g_lim) %#ok
% min_f = min(frame(:));
% frame = (frame - min_f) / (max(frame(:))-min_f); 
% thresh = numel(frame)*thresh_pct;
% is_valid = sum(frame(:) < lower_g_lim | frame(:) > upper_g_lim) < thresh;

is_valid = sum(frame(:) > 200) < 1000;
if is_valid
    frame_r = imresize(frame, 1/16);
    [~, p] = HartigansDipSignifTest(frame_r(:), 500);
    is_valid = p > 0.0001;
end

%**************************************************************************
%**************************************************************************
function plot_segment_traces(pixels_per_mm, segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness)

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
    
hw = 320 / pixels_per_mm;
hh = 240 / pixels_per_mm;
frame_rect_x = [-hw hw hw -hw -hw];
frame_rect_y = [-hh -hh hh hh -hh];

figure; hold all; axis equal; a1 = gca;
plot(a1, motor_x, motor_y, 'b.', 'markersize', 4);
%
for i_seg = 1:length(segments_s)

    frame_idx_i = segments_s{i_seg};

    seg_x = motor_x(frame_idx_i);
    seg_y = motor_y(frame_idx_i);

    cx = median(seg_x);
    cy = median(seg_y);

    plot(a1, cx, cy, 'rx');
    plot(a1, cx+frame_rect_x, cy+frame_rect_y, 'r', 'linewidth', 2);  
end