function [] = ...
    save_vessel_flow_patches(bounding_boxes, segment_apex_idx, segment_offset, frame_dims, frame_dir, save_path, do_plot)
%COMPUTE_VESSEL_FLOW *Insert a one line summary here*
%   [segment_flow_sum segment_flow_count segment_flow_max] =...
%           compute_vessel_flow(segment_data, vessel_data, frame_dir)
%
% Inputs:
%      segment_data - *Insert description of input variable here*
%
%      vessel_data - *Insert description of input variable here*
%
%
% Outputs:
%      mosaic_flow - *Insert description of input variable here*
%
%      vessel_data - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 18-Nov-2015
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

if ~exist('do_plot', 'var')
    do_plot = 0;
end

inner_frame_w = frame_dims(1);
inner_frame_h = frame_dims(2);

frame_list = dir([frame_dir 'frame*.png']);
if isempty(frame_list)
    return;
end

num_frames = length(frame_list);
discard = false(num_frames,1);
for i_f = 1:num_frames
    if length(frame_list(i_f).name) < 14
        discard(i_f) = 1;
    else
        frame = imread([frame_dir frame_list(i_f).name]);
        if i_f == 1
            [frame_h frame_w] = size(frame);
            frames = zeros(frame_h,frame_w,num_frames, 'uint8');
        end
        frames(:,:,i_f) = frame;
    end
end

if ~exist(frame_dir, 'dir') || ~exist([frame_dir 'segment_mask.png'], 'file')
    return;
end

segment_mask = imread([frame_dir 'segment_mask.png']);
[frame_h frame_w] = size(segment_mask);

%Adjust bounding boxes for difference between full frame size and cropped
%inner size
frames_offset(1) = segment_offset(1) - (frame_w - inner_frame_w)/2;
frames_offset(2) = segment_offset(2) - (frame_h - inner_frame_h)/2;

bounding_boxes(:,1,:) = bounding_boxes(:,1,:) - frames_offset(1);
bounding_boxes(:,2,:) = bounding_boxes(:,2,:) - frames_offset(2);

if any(discard)
    display(['Warning: old frame encountered in ' save_path]);
    frames(:,:,discard) = [];
    num_frames = size(frames,3); %#ok
end

if do_plot
    figure; imgray(mean(frames,3)); a1 = gca;
end

num_frames = 0;
for i_b = 1:size(bounding_boxes,3)
    
%     if ~exist([save_path '_v' zerostr(i_b,2) '.mat'], 'file')
%         continue;
%     end
%     f = load([save_path '_v' zerostr(i_b,2) '.mat'], 'x_min');
%     if isfield(f, 'x_min')
%         return;
%     end
%     if ~num_frames
%         load([save_path '_v' zerostr(i_b,2) '.mat'], 'frames_i')
%         num_frames = size(frames_i,3);
%     end
    
    x_min = round(max(bounding_boxes(1,1,i_b),1));
    y_min = round(max(bounding_boxes(1,2,i_b),1));

    x_max = round(min(bounding_boxes(3,1,i_b), frame_w));
    y_max = round(min(bounding_boxes(2,2,i_b), frame_h));

    r = rem(y_max - y_min + 1, 4);
    if r
        y_min = y_min + r;
    end
    r = rem(x_max - x_min + 1, 4);
    if r
        x_min = x_min + r;
    end

    if do_plot
        plot(a1, [x_min x_max x_max x_min x_min], [y_min y_min y_max y_max y_min], 'linewidth', 2);
    end
    
    apex_idx = segment_apex_idx(i_b); %#ok
    frames_i = frames(y_min:y_max,x_min:x_max,:);
    
    rows = any(frames_i(:,:,1),2);
    cols = any(frames_i(:,:,1),1);
    cropped_frames = frames_i(rows, cols, :);

    [vessel_transforms] = ...
        register_tiles_features(cropped_frames, ...
                            'ref_type', 'mosaic',...
                            'region', 'all',...
                            'theta_range', 0, ...
                            'offset_lim', 20, ...
                            'mosaic', mean(cropped_frames,3),...
                            'sigma', 6,...
                            'tile_masks', [],...
                            'debug', 0);
    [~, ~, ~, cleaned_frames, edge_mask] = ...
        create_mosaic(cropped_frames, vessel_transforms);
    
    rows = ~all(edge_mask,2);
    cols = ~all(edge_mask);
    cropped_frames = cleaned_frames(rows, cols, :); %#ok
    
    [~, ~, dx, dy] = ...
        mosaic_limits([y_max-y_min x_max-x_min]+1, vessel_transforms);   
    
    mask_row = any(~edge_mask,1);
    mask_col = any(~edge_mask,2);
    x_min_cropped = floor(x_min + find(mask_row,1) - dx);
    y_min_cropped = floor(y_min + find(mask_col,1) - dy);
    x_max_cropped = x_min_cropped + sum(mask_row) - 1; %#ok
    y_max_cropped = y_min_cropped + sum(mask_col) - 1; %#ok
    
    save([save_path '_v' zerostr(i_b,2) '.mat'], 'frames_i', 'apex_idx', ...
        'vessel_transforms', 'cropped_frames', 'edge_mask',...
        'x_min', 'x_max', 'y_min', 'y_max', 'frames_offset', 'num_frames', ...
        'x_min_cropped', 'x_max_cropped', 'y_min_cropped', 'y_max_cropped');
    
    
end

    