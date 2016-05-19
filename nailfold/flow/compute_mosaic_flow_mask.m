%OBSOLETE
function [mosaic_flow_mask, vessel_data] = compute_mosaic_flow_mask(vessel_data, mosaic_flow, varargin)
%COMPUTE_MOSAIC_FLOW *Insert a one line summary here*
%   [mosaic_flow, vessel_data] = compute_mosaic_flow(sequence_dir, vessel_data)
%
% Inputs:
%      sequence_dir - directory containing the sequence to be processed
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
    'plot', 0,...
    'vessel_thresh', 0.5);



[sequence_data] = read_processed_sequence_from([sequence_dir 'sequence_data\sequence_data.dat']);
vessel_pred = imread([sequence_dir 'capillary_data\vessels_v_pred.png']);
vessel_ori = imread([seq_dir 'capillary_data\vessels_o_pred.png']);

mosaic_w = sequence_data.final_mosaic_size(1);
mosaic_h = sequence_data.final_mosaic_size(2);

vessel_pred = imresize(vessel_pred, [mosaic_h mosaic_w]);
vessel_ori = imresize(vessel_ori, [mosaic_h mosaic_w]);

vessel_pred = padarray(vessel_pred, mosaic_flow.padding([3 1]), 0, 'pre');
vessel_pred = padarray(vessel_pred, mosaic_flow.padding([4 2]), 0, 'post');

flow_w = mosaic_w + mosaic_flow.padding(1) + mosaic_flow.padding(2);
flow_h = mosaic_h + mosaic_flow.padding(3) + mosaic_flow.padding(4);
mosaic_flow_mask = false(flow_h, flow_w);

seg_w = sequence_data.inner_frame_dims(1);
seg_h = sequence_data.inner_frame_dims(2);
num_segs = sequence_data.num_segments;

%Allocate output containers
cxy = cumsum(sequence_data.segment_displacements);
oxy = min(cxy);
mosaic_displacements = bsxfun(@minus, cxy, oxy);

distal_xy = vessel_data.apex_measures.distal.apex_xy / ...
    vessel_data.resize_factor;
discard_apices = isnan(vessel_data.apex_measures.distal.apex_xy(:,1));
distal_xy(discard_apices,:) = [];

num_vessels = size(distal_xy,1);
flow_measurable = false(num_vessels,1);

bounding_boxes = vessel_data.apex_measures.distal.bounding_box / ...
    vessel_data.resize_factor;
bounding_boxes(:,:,discard_apices) = [];
    
%Loop through each segment
for i_seg = 1:num_segs
    
    if ~sequence_data.segment_required(i_seg)
        continue;
    end
    
    %Check it is a segment with valid registered frames
    frame_dir = [sequence_dir 'registered' zerostr(i_seg-1,2) '\'];
    if ~exist(frame_dir, 'dir') || ~exist([frame_dir 'segment_mask.png'], 'file')
        continue;
    end
    
    segment_mask = imread([frame_dir 'segment_mask.png']);
    [reg_h reg_w] = size(segment_mask);
    
    %Workout which apices lie in this segment, and get bound boxes
    %(adjusted to segment coordinate frame) for them 
    segment_start_x = mosaic_displacements(i_seg,1);
    segment_start_y = mosaic_displacements(i_seg,2);
    
    segment_apexes_x = distal_xy(:,1) - segment_start_x;
    segment_apexes_y = distal_xy(:,2) - segment_start_y;

    in_segment = ...
        segment_apexes_x >= 1 & segment_apexes_x <= seg_w &...
        segment_apexes_y >= 1 & segment_apexes_y <= seg_h;
   
    segment_boxes = ...
        bounding_boxes(:,:,in_segment);
    
    if isempty(segment_boxes)
        continue;
    end
    
    offset_x = round((reg_w - seg_w) / 2);
    offset_y = round((reg_h - seg_h) / 2);
    
    reg_start_x = segment_start_x - offset_x;
    reg_start_y = segment_start_y - offset_y;
    
    %Otherwise continue and load in the flow segments
    segment_boxes(:,1,:) = segment_boxes(:,1,:) - reg_start_x;
    segment_boxes(:,2,:) = segment_boxes(:,2,:) - reg_start_y;
    
    segment_xy = distal_xy(in_segment,:);
    segment_xy(:,1) = segment_xy(:,1) - reg_start_x;
    segment_xy(:,2) = segment_xy(:,2) - reg_start_y;
    
    num_segment_vessels = size(segment_boxes,3);
    segment_flow_vessels = false(num_segment_vessels,1);
    
    for i_b = 1:num_segment_vessels
        
        x_min = round(max(segment_boxes(1,1,i_b),1));
        y_min = round(max(segment_boxes(1,2,i_b),1));

        x_max = round(min(segment_boxes(3,1,i_b), reg_w));
        y_max = round(min(segment_boxes(2,2,i_b), reg_h));

        r = rem(y_max - y_min + 1, 4);
        if r
            y_min = y_min + r;
        end
        r = rem(x_max - x_min + 1, 4);
        if r
            x_min = x_min + r;
        end
        seg_rows = y_min:y_max;
        seg_cols = x_min:x_max;
        flow_rows = seg_rows + reg_start_y + mosaic_flow.padding(3) - 1;
        flow_cols = seg_cols + reg_start_x + mosaic_flow.padding(1) - 1;
        flow_mask_patch = segment_mask(seg_rows, seg_cols);
        vessel_mask_patch = vessel_pred(flow_rows, flow_cols) > args.vessel_thresh;
        
        cx = segment_xy(i_b,1) - x_min;
        cy = segment_xy(i_b,2) - y_min;
        
        mask_patch = flow_mask_patch & bwselect(vessel_mask_patch, cx, cy, 4);
        
        mosaic_flow_mask(flow_rows, flow_cols) = mosaic_flow_mask(flow_rows, flow_cols) | ...
            mask_patch;
        
        if any(mask_patch(:))
            segment_flow_vessels(i_b) = 1;
        end
    end
    flow_measurable(in_segment) = flow_measurable(in_segment) | segment_flow_vessels;
end
num_vessels = size(discard_apices,1);
vessel_data.apex_measures.distal.flow_measurable = false(num_vessels,1);
vessel_data.apex_measures.distal.flow_measurable(~discard_apices) = flow_measurable;
% 
if args.plot
    figure;
    a1 = subplot(2,1,1); imgray(mosaic_flow_mask); 
    a2 = subplot(2,1,2); show_flow_as('rgb', mosaic_flow.max); hold all;
    linkaxes([a1 a2]);
end

flow_measurable = vessel_data.apex_measures.distal.flow_measurable;
[flow_h flow_w] = size(mosaic_flow.max);

num_vessels = length(flow_measurable);
mean_flow = nan(num_vessels,1);
median_flow = nan(num_vessels,1);
for i_ve = 1:num_vessels
    if ~flow_measurable(i_ve)
        continue;
    end
    bounding_box = vessel_data.apex_measures.distal.bounding_box(:,:,i_ve) / ...
        vessel_data.resize_factor;
    
    x_min = round(max(bounding_box(1,1) + mosaic_flow.padding(1),1));
    y_min = round(max(bounding_box(1,2) + mosaic_flow.padding(2),1));

    x_max = round(min(bounding_box(3,1) + mosaic_flow.padding(1), flow_w));
    y_max = round(min(bounding_box(2,2) + mosaic_flow.padding(2), flow_h));
    
    vessel_rows = y_min:y_max;
    vessel_cols = x_min:x_max;
        
    flow_patch = mosaic_flow.max(vessel_rows, vessel_cols);
    mask_patch = mosaic_flow_mask(vessel_rows, vessel_cols);
    flow_velocity = abs(flow_patch(mask_patch));
    mean_flow(i_ve) = mean(flow_velocity);
    median_flow(i_ve) = median(flow_velocity);
    
    if args.plot
        plot(a1,...
            [x_min x_max x_max x_min x_min], [y_min y_min y_max y_max y_min]);
        plot(a2,...
            [x_min x_max x_max x_min x_min], [y_min y_min y_max y_max y_min]);
        text((x_min+x_max)/2, (y_min+y_max)/2, num2str(mean_flow(i_ve),3));
    end
end
    
vessel_data.apex_measures.distal.mean_flow = mean_flow;
vessel_data.apex_measures.distal.median_flow = median_flow;
    
    