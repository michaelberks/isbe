function [] = save_mosaic_flow_patches(sequence_dir, vessel_data, save_dir, varargin)
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

    figure; imgray(mosaic); 
    a1 = gca; caxis([gmin gmax]);
end

cxy = cumsum(sequence_data.segment_displacements);
oxy = min(cxy);
mosaic_displacements = bsxfun(@minus, cxy, oxy);

distal_xy = vessel_data.apex_measures.distal.apex_xy / vessel_data.resize_factor;
discard_apices = isnan(vessel_data.apex_measures.distal.apex_xy(:,1));

distal_xy(discard_apices,:) = [];
apex_idx = find(~discard_apices);



%Loop through each segment
for i_seg = 1:num_segs
    
    if ~sequence_data.segment_required(i_seg)
        continue;
    end
    
    %Check it is a segment with valid registered frames
    frame_dir = [sequence_dir 'registered' zerostr(i_seg-1,2) '\'];
    if ~exist(frame_dir, 'dir')
        continue;
    end
    
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
    
    segment_apex_idx = apex_idx(in_segment);

    bounding_boxes = vessel_data.apex_measures.distal.bounding_box;
    bounding_boxes(:,:,discard_apices) = [];
    bounding_boxes = ...
        bounding_boxes(:,:,in_segment)/vessel_data.resize_factor;
    
    if ~isempty(bounding_boxes)
        save_path = [save_dir '_s' zerostr(i_seg,2)];
        save_vessel_flow_patches(bounding_boxes, segment_apex_idx, [segment_start_x segment_start_y],...
            sequence_data.inner_frame_dims, frame_dir, save_path, args.plot);
    end
    
end

    
    