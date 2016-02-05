function [segment_flow_sum segment_flow_count segment_flow_max] = ...
    compute_vessel_flow(i_seg, segment_data, vessel_data, frame_dir)
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

frame_w = 640;
frame_h = 480;

segment_start_x = segment_data.segment_transforms(1,3,i_seg);
segment_start_y = segment_data.segment_transforms(2,3,i_seg);

distal_xy = vessel_data.apex_measures.distal.apex_xy / vessel_data.resize_factor;

segment_apexes_x = distal_xy(:,1) - segment_start_x;
segment_apexes_y = distal_xy(:,2) - segment_start_y;

in_segment = ...
    segment_apexes_x >= 1 & segment_apexes_x <= 640 &...
    segment_apexes_y >= 1 & segment_apexes_y <= 480;

bounding_boxes = vessel_data.apex_measures.distal.bounding_box(:,:,in_segment)/vessel_data.resize_factor;
bounding_boxes(:,1,:) = bounding_boxes(:,1,:) - segment_start_x;
bounding_boxes(:,2,:) = bounding_boxes(:,2,:) - segment_start_y;

figure; imgray(segment_data.segment_mosaic);
plot(distal_xy(:,1), distal_xy(:,2), 'rx');
plot(...
    [segment_start_x segment_start_x+frame_w segment_start_x+frame_w segment_start_x segment_start_x],...
    [segment_start_y segment_start_y segment_start_y+frame_h segment_start_y+frame_h segment_start_y]);

%%
%[s_rows s_cols] = size(segment_data.segment_frames(:,:,1));
%sx = repmat((segment_start_x+(0:s_cols-1)), s_rows, 1)*vessel_data.resize_factor;
%sy = repmat((segment_start_y+(0:s_rows-1))', 1, s_cols)*vessel_data.resize_factor;
%mosaic_vessels = interp2(vessel_data.vessel_predictions(:,:,1), sx, sy, '*linear', 0);
%mosaic_ori = interp2(vessel_data.vessel_predictions(:,:,2), sx, sy, '*linear', 0);

%offset_yx = (size(segment_data.segment_mosaics{i_seg}) - size(segment_data.segment_frames(:,:,1)))/2;
frame_list = dir([frame_dir 'frame*.png']);
num_frames = length(frame_list);


for i_f = 1:num_frames
    frame = imread([frame_dir frame_list(i_f).name]);
    if i_f == 1
        [frame_h frame_w] = size(frame);
        frames = zeros(frame_h,frame_w,num_frames);
    end
    frames(:,:,i_f) = frame;
end

figure; imgray(mean(frames,3)); a1 = gca;

segment_flow_sum = zeros(frame_h, frame_w);
segment_flow_max = zeros(frame_h, frame_w);
segment_flow_count = zeros(frame_h, frame_w);
%
for i_b = 1:sum(in_segment)
    
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

    
    plot(a1, [x_min x_max x_max x_min x_min], [y_min y_min y_max y_max y_min], 'linewidth', 2);
    
    vessel_dir = [frame_dir 'vessel' zerostr(i_b,2) '\'];
    [flowPyramidEst, flowConfidence] = ...
        estimate_flow_multilevel(frames(y_min:y_max,x_min:x_max,:), [], vessel_dir, 1:3);
    
    segment_flow_sum(y_min:y_max,x_min:x_max) = ...
        segment_flow_sum(y_min:y_max,x_min:x_max) + flowPyramidEst{1};
    
    segment_flow_max(y_min:y_max,x_min:x_max) = ...
        max(segment_flow_max(y_min:y_max,x_min:x_max), flowPyramidEst{1});
    
    segment_flow_count(y_min:y_max,x_min:x_max) = ...
        segment_flow_count(y_min:y_max,x_min:x_max) + 1;
    
    
    
end

    