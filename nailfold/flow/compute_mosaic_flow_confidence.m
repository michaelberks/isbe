function [mosaic_flow] = compute_mosaic_flow_confidence(sequence_dir, vessel_data, save_dir, varargin)
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
    'plot', 0);

[sequence_data] = read_processed_sequence_from([sequence_dir 'sequence_data\sequence_data.dat']);

mosaic_w = sequence_data.final_mosaic_size(1);
mosaic_h = sequence_data.final_mosaic_size(2);

seg_w = sequence_data.inner_frame_dims(1);
seg_h = sequence_data.inner_frame_dims(2);
num_segs = sequence_data.num_segments;

if args.plot
    mosaic = imread([sequence_dir 'sequence_data\full_mosaic.png']);
    count_map = imread([sequence_dir 'sequence_data\count_map.png']);
    mask = count_map > 0;
    mask([1 end],:) = 0;
    mask(:, [1 end]) = 0;
    mask = imerode(mask, strel('disk', 2));

    gmin = min(mosaic(mask));
    gmax = max(mosaic(mask));

    f1 = figure; 
    a1 = subplot(2,1,1); 
    imgray(mosaic); 
    caxis([gmin gmax]);
end

%Allocate output containers
mosaic_flow.confidence = zeros(mosaic_h, mosaic_w);
mosaic_flow.adjusted_confidence = zeros(mosaic_h, mosaic_w);
mosaic_flow.padding = zeros(4,1);

cxy = cumsum(sequence_data.segment_displacements);
oxy = min(cxy);
mosaic_displacements = bsxfun(@minus, cxy, oxy);

distal_xy = vessel_data.apex_measures.distal.apex_xy / vessel_data.resize_factor;
discard_apices = isnan(vessel_data.apex_measures.distal.apex_xy(:,1));
distal_xy(discard_apices,:) = [];

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
    
    if args.plot
        plot(a1, ...
            [segment_start_x segment_start_x+seg_w segment_start_x+seg_w segment_start_x segment_start_x],...
            [segment_start_y segment_start_y segment_start_y+seg_h segment_start_y+seg_h segment_start_y]);
    end
    
    segment_apexes_x = distal_xy(:,1) - segment_start_x;
    segment_apexes_y = distal_xy(:,2) - segment_start_y;

    in_segment = ...
        segment_apexes_x >= 1 & segment_apexes_x <= seg_w &...
        segment_apexes_y >= 1 & segment_apexes_y <= seg_h;

    bounding_boxes = vessel_data.apex_measures.distal.bounding_box;
    bounding_boxes(:,:,discard_apices) = [];
    bounding_boxes = ...
        bounding_boxes(:,:,in_segment)/vessel_data.resize_factor;
    
    if isempty(bounding_boxes)
        continue;
    end
    
    %Otherwise continue and load in the flow segments
    bounding_boxes(:,1,:) = bounding_boxes(:,1,:) - segment_start_x;
    bounding_boxes(:,2,:) = bounding_boxes(:,2,:) - segment_start_y;
        
    save_path = [save_dir '_s' zerostr(i_seg,2)];
    
%     Estimate flow for all capillaries in the segment
    [segment_flow_confidence segment_flow_adjusted_confidence] =...
        compute_vessel_flow_confidence(bounding_boxes, [reg_w reg_h], [seg_w seg_h], save_path, args.plot > 1);
    
    %Workout where these should be saved in the main mosaic
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
        mosaic_flow.adjusted_confidence = [zeros(mosaic_h,pad_val) mosaic_flow.adjusted_confidence];
        mosaic_flow.confidence = [zeros(mosaic_h,pad_val) mosaic_flow.confidence];
        reg_start_x = reg_start_x + pad_val;
        reg_end_x = reg_end_x + pad_val;
        mosaic_w = mosaic_w + pad_val;
        mosaic_flow.padding(1) = mosaic_flow.padding(1) + pad_val;
    end
    if reg_start_y < 1
        pad_val = abs(reg_start_y)+1;
        mosaic_flow.adjusted_confidence = [zeros(pad_val, mosaic_w); mosaic_flow.adjusted_confidence];
        mosaic_flow.confidence = [zeros(pad_val,mosaic_w); mosaic_flow.confidence];
        reg_start_y = reg_start_y + pad_val;
        reg_end_y = reg_end_y + pad_val;
        mosaic_h = mosaic_h + pad_val;
        mosaic_flow.padding(3) = mosaic_flow.padding(3) + pad_val;
    end
    if reg_end_x > mosaic_w
        mosaic_flow.padding(2) = mosaic_flow.padding(3) + reg_end_x - mosaic_w;
        mosaic_w = reg_end_x;
        mosaic_flow.adjusted_confidence(mosaic_h, mosaic_w) = 0;
        mosaic_flow.confidence(mosaic_h, mosaic_w) = 0;
    end
    if reg_end_y > mosaic_h
        mosaic_flow.padding(4) = mosaic_flow.padding(4) + reg_end_y - mosaic_h;
        mosaic_h = reg_end_y;
        mosaic_flow.adjusted_confidence(mosaic_h, mosaic_w) = 0;
        mosaic_flow.confidence(mosaic_h, mosaic_w) = 0;
    end
    
    %Update the full mosaic sum, count and max flow maps with the segment flow maps
    mosaic_flow.adjusted_confidence(reg_start_y:reg_end_y, reg_start_x:reg_end_x) = ...
        max(mosaic_flow.confidence(reg_start_y:reg_end_y, reg_start_x:reg_end_x), segment_flow_adjusted_confidence);
 
    mosaic_flow.confidence(reg_start_y:reg_end_y, reg_start_x:reg_end_x) = ...
        max(mosaic_flow.confidence(reg_start_y:reg_end_y, reg_start_x:reg_end_x), segment_flow_confidence);
    
    if args.plot
        plot(a1, ...
            [reg_start_x reg_start_x+reg_w reg_start_x+reg_w reg_start_x reg_start_x],...
            [reg_start_y reg_start_y reg_start_y+reg_h reg_start_y+reg_h reg_start_y], '--');
    end
end

if args.plot
    figure(f1);
    a2 = subplot(2,1,2);
    show_flow_as('rgb', mosaic_flow.confidence);
    linkaxes([a1 a2]);
end

    
    