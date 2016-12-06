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
    'skip_dirt_frames',         10,...
    'selected_segments',        [],...
    'pixels_per_mm',            1000,...
    'min_stationary_frames',    120,...
    'num_segment_tile_frames',  30,... 
    'min_focus_pct',            0.95,...
    'max_focus_val',            50,...
    'bad_frame_thresh',         0.1,...
    'lower_g_lim',              50,...
    'upper_g_lim',              200,...
    'mosaic_lims_x',            -310:309,...
    'mosaic_lims_y',            -230:229,...
    'frame_h',                  480,...
    'frame_w',                  640,...
    'sigma',                    6,...
    'theta_range',              0,...
    'intra_segment_offset',     120,...
    'inter_segment_offset',     240,...
    'overwrite',                0,...
    'make_videos',              0,...
    'delete_frames',            0,...
    'plot',                     0,...
    'debug',                    0);
clear varargin;

tic;

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
    tic;
    display('Constructing dirt image');
    frame_sum = zeros(args.frame_h, args.frame_w);
    frame_diff = zeros(args.frame_h, args.frame_w);
    frame_prev = [];
    total_frames = 0;
    for i_seg = 1:length(segments_ns)

        frame_idx = segments_ns{i_seg};
        num_frames = length(frame_idx);
        
        for i_im = 1:args.skip_dirt_frames:num_frames
            frame = double(rot90(imread([args.sequence_dir 'frame' zerostr(frame_idx(i_im), 5) '.bmp']),2));
            if is_frame_ok(frame, args.bad_frame_thresh, args.lower_g_lim, args.upper_g_lim)
                frame_sum = frame_sum + frame;
                
                if isempty(frame_prev)
                    frame_prev = frame;
                else
                    frame_diff = frame_diff + abs(frame-frame_prev);
                end
                total_frames = total_frames + 1;
            end
        end
    end
    args.dirt_image = frame_sum / total_frames;  
    toc;
    
elseif ischar(args.dirt_image) && exist(args.dirt_image, 'file')
    s = load(args.dirt_image, 'dirt_image');
    if isfield(s, 'dirt_image')
        args.dirt_image = s.dirt_image;
    end
    clear s;
end
dirt_image = args.dirt_image;
dirt_image = dirt_image - mean(dirt_image(:)); 

h = fspecial('sobel');
dirt_grad_x = imfilter(dirt_image, h, 'replicate');
dirt_grad_y = imfilter(dirt_image, h', 'replicate');
dirt_grad_mag = sqrt(dirt_grad_x.^2 + dirt_grad_y.^2);
grad_threshold = 2*mean(dirt_grad_mag(:));
dirt_mask = dirt_grad_mag < grad_threshold;

create_folder([args.sequence_dir 'segments\']);

progress = 1; %#ok
if ~args.overwrite && exist([args.sequence_dir 'segments\segment_data.mat'], 'file')
    save([args.sequence_dir 'segments\segment_data.mat'],...
        'dirt_image', 'progress', '-append');
else
    save([args.sequence_dir 'segments\segment_data.mat'],...
        'dirt_image', 'progress');
end
clear dirt_image;
%--------------------------------------------------------------------------
%Make segment mosaics
if isempty(args.selected_segments)
    args.selected_segments =  1:num_segs;
    
    segment_frames = zeros(length(args.mosaic_lims_y), length(args.mosaic_lims_x), num_segs);
    segment_masks = false(size(segment_frames));
    segment_mosaics = cell(num_segs,1);
    focused_idx = cell(num_segs, 1);
    focused_transforms = cell(num_segs,1);
    reject_segments = false(num_segs,1);
    frame_centres = zeros(num_segs, 2);
    matched_counts = cell(num_segs,1);
else
    load([args.sequence_dir 'segments\segment_data.mat'],...
        'segment_frames', 'segment_mosaics', 'focused_idx', 'focused_transforms',...
        'reject_segments', 'frame_centres', 'matched_counts');%
end

n_segs_str = num2str(length(args.selected_segments));
for i_seg = args.selected_segments
    
    display(['Processing segment ' num2str(i_seg) ' of ' n_segs_str]);
    
    frame_idx = segments_s{i_seg};
    num_frames = length(frame_idx);

    seg_x = motor_x(frame_idx);
    seg_y = motor_y(frame_idx);
    seg_z = motor_z(frame_idx);
    seg_s = sharpness(frame_idx);

    frame_centres(i_seg,1) = seg_x(round(num_frames/2));
    frame_centres(i_seg,2) = seg_y(round(num_frames/2));
    
   %Extract the largest consectutive sequence of focused frames   
    dz = diff(seg_z);
    ni = find(dz(end:-1:1), 1);
    if isempty(ni)
        ni = num_frames;
    end
    focused_idx{i_seg} = num_frames - ni + (1:ni); 
    num_frames = length(focused_idx{i_seg});
    
    if 1==1 || num_frames < args.min_stationary_frames
        continue;
    end
    
    %Re-register all the focused frames
    focused_idx_i = frame_idx(focused_idx{i_seg});   
    frames = zeros(args.frame_h,args.frame_w,num_frames);
    reject_frames = false(num_frames,1);
    for i_im = 1:num_frames
        frame = ...
            double(rot90(imread([args.sequence_dir 'frame' zerostr(focused_idx_i(i_im), 5) '.bmp']),2));
        frames(:,:,i_im) = frame - args.dirt_image;
        
        if ~is_frame_ok(frame, args.bad_frame_thresh, args.lower_g_lim, args.upper_g_lim)
            reject_frames(i_im) = 1;
        end
    end
    
    %We can't just discard frames or it will mess up the flow videos
    %However, if we don't have enough good frames to make a video with,
    %don't both
    good_frames = sum(~reject_frames);
    target_frames = min(args.num_segment_tile_frames, round(num_frames/4));
    
    if good_frames > args.min_stationary_frames
        build_video = true;
        
    elseif good_frames >= target_frames
        build_video = false;
        frames = frames(:,:,~reject_frames);
        
    else %not enough good frames to meet target - potentially discard segment
      
        reject_segments(i_seg) = 1;
        
        %May as well just load in the final 30 frames to make some form of
        %mosaic
        for i_im = 1:target_frames
            frame = double(rot90(imread([args.sequence_dir 'frame'...
                zerostr(frame_idx(end+1-i_im), 5) '.bmp']),2));
            frames(:,:,i_im) = frame  - args.dirt_image;
        end
        
        if args.debug
            figure; imgray(mean(frames,3));
            title(['Segment ' num2str(i_seg) ' failed.']);
            xlabel(args.sequence_dir);
        end
    end

    tile_masks = repmat(dirt_mask, 1, 1, size(frames,3));
    [compound_transforms, matched_counts{i_seg}] = ...
        register_tiles_features(frames, ...
                            'theta_range', args.theta_range, ...
                            'offset_lim', args.intra_segment_offset, ...
                            'sigma', args.sigma,...
                            'tile_masks', tile_masks,...
                            'debug', false);

    tic;
    [segment_mosaics{i_seg}] = ...
        create_mosaic(frames, compound_transforms);
    toc;
    tic;
    [segment_mosaics{i_seg}, ~, ~, cleaned_frames, edge_mask] = ...
        create_mosaic(frames, compound_transforms);
    toc;
    
    mask_mosaic = create_mosaic(double(tile_masks), compound_transforms);

    [r, c] = size(segment_mosaics{i_seg});
    cx = round(c/2);
    cy = round(r/2);
    segment_frames(:,:,i_seg) = segment_mosaics{i_seg}(cy+args.mosaic_lims_y, cx+args.mosaic_lims_x);
    segment_masks(:,:,i_seg) = mask_mosaic(cy+args.mosaic_lims_y, cx+args.mosaic_lims_x);
    
    
    
    %
    [focused_transforms{i_seg}] = ...
        register_tiles_features(frames, ...
            'ref_type', 'mosaic',...
            'sigma', args.sigma,...
            'mosaic', segment_frames(:,:,i_seg),...
            'compound_transforms', [],...
            'theta_range', args.theta_range, ...
            'offset_lim', args.intra_segment_offset, ...
            'tile_masks', repmat(dirt_mask, 1, 1, size(frames,3)),...
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
    if args.make_videos
        system(cmd);
    else
        fid1 = fopen([reg_folder 'make_movie.bat'], 'wt');
        fprintf(fid1,'%s', cmd); 
        fclose(fid1);
    end
    if args.delete_frames
        delete([reg_folder 'frame*.png']);
    end
end

%--------------------------------------------------------------------------
%Make complete mosaic
progress = 2; %#ok
save([args.sequence_dir 'segments\segment_data.mat'],...
    'segment_frames', 'segment_mosaics', 'focused_idx', 'focused_transforms',...
    'reject_segments', 'frame_centres', 'matched_counts', 'progress',  '-append');

%Check what to do with rejected segments
required_segments = true(num_segs,1);
if any(reject_segments)
    frame_w = args.frame_w / args.pixels_per_mm;
    frame_h = args.frame_h / args.pixels_per_mm;
    
    %Do we have a completely connected mosaic without this segment   
    for i_r = find(reject_segments)
        required_segments(i_r) = 0;
        is_connected = is_mosaic_connected(frame_centres(required_segments,:), frame_w, frame_h);
        
        %If mosaic is not connected without this segment, we don't discard
        %it
        if ~is_connected
            required_segments(i_r) = 1;
        end
    end
    %segment_frames(:,:,reject_segments) = [];
    %segment_mosaics(reject_segments) = []; %#ok
    %frame_centres(reject_segments,:) = [];
end

if size(frame_centres, 1) < 2
    return;
end

%Use motor positions to initialise alignment of frames
offset_centres = [diff(frame_centres(required_segments,:)*args.pixels_per_mm); 0 0];
offset_centres(:,2) = -offset_centres(:,2);

[segment_transforms] = ...
        register_tiles_features(segment_frames(:,:,required_segments), ...
                            'theta_range', args.theta_range, ...
                            'offset_lim', args.inter_segment_offset, ...
                            'offset_centres', offset_centres,...
                            'sigma', args.sigma,...
                            'tile_masks', segment_masks,...
                            'debug', false);
[segment_mosaic, segment_mosaic_weights, segment_transforms] = ...
    create_mosaic(segment_frames, segment_transforms, 'diamond'); %#ok
segment_mask = segment_mosaic_weights > 0.1; %#ok

progress = 3; %#ok
save([args.sequence_dir 'segments\segment_data.mat'],...
    'segment_mosaic', 'segment_mosaic_weights',...
    'segment_transforms', 'segment_mask', 'required_segments', 'progress', '-append');

display('Total time to register sequence');
toc;
%--------------------------------------------------------------------------
%End of main function
if args.plot
    figure; imgray(segment_mosaic);
    plot_segment_traces(segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness,...
        args.frame_w/args.pixels_per_mm, args.frame_h/args.pixels_per_mm);
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
[r c] = size(frame);
rows = round(r/4):round(3*r/4);
cols = round(c/4):round(3*c/4);
frame = frame(rows,cols);
thresh = numel(frame)*thresh_pct;

is_valid = sum(frame(:) > upper_g_lim) < thresh;
% if is_valid
%     frame_r = imresize(frame, 1/16);
%     [~, p] = HartigansDipSignifTest(frame_r(:), 500);
%     is_valid = p > 0.0001;
% end

%**************************************************************************
%**************************************************************************
function plot_segment_traces(segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness, frame_w, frame_h)

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
    
hw = 0.5*frame_w;
hh = 0.5*frame_h;
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

function is_connected = is_mosaic_connected(frame_centres, frame_w, frame_h)

num_segs = size(frame_centres,1);
x_diffs = zeros(num_segs);
y_diffs = zeros(num_segs);
for i_seg = 1:num_segs
    x_diffs(:,i_seg) = abs(frame_centres(:,1) - frame_centres(i_seg,1));
    y_diffs(:,i_seg) = abs(frame_centres(:,2) - frame_centres(i_seg,2));
end

num_components = 0;
adjacency_matrix = (x_diffs < frame_w) & (y_diffs < frame_h);

%Result array.
marked = zeros(num_segs, 1);

%Enumerate all vertices, 
for i_vi = 1:num_segs

    %if for vertex number i, marks[i] == 0 then
    if ~marked(i_vi)
        
        %Increment components
        num_components = num_components+1;

        %Put this vertex into queue, and 
        queue = i_vi;
        
        while (~isempty(queue)) 
            %Pop current vertex from queue
        	current = queue(end);
            queue(end) = [];

            %Add all adjacent and not currently marked vertices to the
            %queue
            for i_vj = 1:num_segs
                if (adjacency_matrix(current, i_vj) && ~marked(i_vj))
                    marked(i_vj) = num_components;
                    queue(end+1) = i_vj; %#ok
                end
            end
        end
    end
end

is_connected = num_components == 1;

