function [mosaic_flow, vessel_data] = compute_mosaic_flow(segment_data, vessel_data, varargin)
%COMPUTE_MOSAIC_FLOW *Insert a one line summary here*
%   [mosaic_flow, vessel_data] = compute_mosaic_flow(segment_data, vessel_data)
%
% Inputs:
%      segment_data - data structure obtained by applying register_sequence
%
%      vessel_data - data structure obtained by applying detect_capillaries
%
%
% Outputs:
%      mosaic_flow - estimated flow across the whole mosaic
%
%      vessel_data - updated vessel data with new flow measures added
%
%
% Example:
%
% Notes:
%
% See also: REGISTER_SEQUENCE DETECT_CAPILLARIES COMPUTE_VESSEL_FLOW
%
% Created: 18-Nov-2015
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
args = u_packargs(varargin, '0',...
    'plot', 0);

[mosaic_h mosaic_w] = size(segment_data.segment_mosaic);
[seg_h seg_w num_segs] = size(segment_data.segment_frames);

if args.plot
    figure; imgray(segment_data.segment_mosaic); a1 = gca;
end

%Allocate output containers
mosaic_flow.sum = zeros(mosaic_h, mosaic_w);
mosaic_flow.count = zeros(mosaic_h, mosaic_w);
mosaic_flow.max = zeros(mosaic_h, mosaic_w);
mosaic_flow.padding = zeros(4,1);

%Loop through each segment
for i_seg = 1:num_segs
    
    %Check it is a segment with valid registered frames
    frame_dir = [segment_data.sequence_dir 'registered' zerostr(i_seg,2) '\'];
    if ~exist(frame_dir, 'dir')
        continue;
    end
        
    %Estimate flow for all capillaries in the segment
%     [segment_flow_sum segment_flow_count segment_flow_max] =...
%         compute_vessel_flow(i_seg, segment_data, vessel_data, frame_dir);
    
    segment_flow_sum = double(imread([frame_dir 'frame0001.png']));
    segment_flow_count = ones(size(segment_flow_sum));
    segment_flow_max = segment_flow_sum;
    [reg_h reg_w] = size(segment_flow_max);
    
    %Workout where these should be saved in the main mosaic
    segment_start_x = round(segment_data.segment_transforms(1,3,i_seg));
    segment_start_y = round(segment_data.segment_transforms(2,3,i_seg));
    
    offset_x = round((reg_w - seg_w) / 2);
    offset_y = round((reg_h - seg_h) / 2);
    
    reg_start_x = segment_start_x - offset_x;
    reg_start_y = segment_start_y - offset_y;
    
    reg_end_x = reg_start_x + reg_w - 1;
    reg_end_y = reg_start_y + reg_h - 1;
    
    %If necessary, we may need to extend the original mosaic - make sure we
    %keep a track of this
    if reg_start_x < 1
        pad_val = abs(reg_start_x)+1;
        mosaic_flow.sum = [zeros(mosaic_h,pad_val) mosaic_flow.sum];
        mosaic_flow.count = [zeros(mosaic_h,pad_val) mosaic_flow.count];
        mosaic_flow.max = [zeros(mosaic_h,pad_val) mosaic_flow.max];
        reg_start_x = reg_start_x + pad_val;
        reg_end_x = reg_end_x + pad_val;
        mosaic_w = mosaic_w + pad_val;
        mosaic_flow.padding(1) = mosaic_flow.padding(1) + pad_val;
    end
    if reg_start_y < 1
        pad_val = abs(reg_start_y)+1;
        mosaic_flow.sum = [zeros(pad_val, mosaic_w); mosaic_flow.sum];
        mosaic_flow.count = [zeros(pad_val,mosaic_w); mosaic_flow.count];
        mosaic_flow.max = [zeros(pad_val,mosaic_w); mosaic_flow.max];
        reg_start_y = reg_start_y + pad_val;
        reg_end_y = reg_end_y + pad_val;
        mosaic_h = mosaic_h + pad_val;
        mosaic_flow.padding(3) = mosaic_flow.padding(3) + pad_val;
    end
    if reg_end_x > mosaic_w
        mosaic_flow.padding(2) = mosaic_flow.padding(3) + reg_end_x - mosaic_w;
        mosaic_w = reg_end_x;
        mosaic_flow.sum(mosaic_h, mosaic_w) = 0;
        mosaic_flow.count(mosaic_h, mosaic_w) = 0;
        mosaic_flow.max(mosaic_h, mosaic_w) = 0;
    end
    if reg_end_y > mosaic_h
        mosaic_flow.padding(4) = mosaic_flow.padding(4) + reg_end_y - mosaic_h;
        mosaic_h = reg_end_y;
        mosaic_flow.sum(mosaic_h, mosaic_w) = 0;
        mosaic_flow.count(mosaic_h, mosaic_w) = 0;
        mosaic_flow.max(mosaic_h, mosaic_w) = 0;
    end
    
    %Update the full mosaic sum, count and max flow maps with the segment flow maps
    mosaic_flow.sum(reg_start_y:reg_end_y, reg_start_x:reg_end_x) = ...
        mosaic_flow.sum(reg_start_y:reg_end_y, reg_start_x:reg_end_x) + segment_flow_sum;
    mosaic_flow.count(reg_start_y:reg_end_y, reg_start_x:reg_end_x) = ...
        mosaic_flow.count(reg_start_y:reg_end_y, reg_start_x:reg_end_x) + segment_flow_count;
    mosaic_flow.max(reg_start_y:reg_end_y, reg_start_x:reg_end_x) = ...
        max(mosaic_flow.max(reg_start_y:reg_end_y, reg_start_x:reg_end_x), segment_flow_max);
    
    if args.plot
        plot(a1, ...
            [segment_start_x segment_start_x+seg_w segment_start_x+seg_w segment_start_x segment_start_x],...
            [segment_start_y segment_start_y segment_start_y+seg_h segment_start_y+seg_h segment_start_y]);
        plot(a1, ...
            [reg_start_x reg_start_x+reg_w reg_start_x+reg_w reg_start_x reg_start_x],...
            [reg_start_y reg_start_y reg_start_y+reg_h reg_start_y+reg_h reg_start_y], '--');
    end
end

    
    