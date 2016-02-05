function [segment_flow_confidence segment_flow_adjusted_confidence] = ...
    compute_vessel_flow_confidence(bounding_boxes, frame_dims, inner_frame_dims, save_path, flow_poly, do_plot)
%COMPUTE_VESSEL_FLOW *Insert a one line summary here*
%   [segment_flow_confidence segment_flow_count segment_flow_adjusted_confidence] =...
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

frame_w = frame_dims(1);
frame_h = frame_dims(2);
inner_frame_w = inner_frame_dims(1);
inner_frame_h = inner_frame_dims(2);

segment_flow_confidence = zeros(frame_h, frame_w);
segment_flow_adjusted_confidence = zeros(frame_h, frame_w);

bounding_boxes(:,1,:) = bounding_boxes(:,1,:) + (frame_w - inner_frame_w)/2;
bounding_boxes(:,2,:) = bounding_boxes(:,2,:) + (frame_h - inner_frame_h)/2;
%
for i_b = 1:size(bounding_boxes,3)
    
    if ~exist([save_path '_v' zerostr(i_b,2) '.mat'], 'file')
        continue;
    end
    
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
    
    load([save_path '_v' zerostr(i_b,2) '.mat'], 'flow_results');
    
    flow_confidence = flow_results.flowPyramidEst{1};
    adjusted_confidence = polyval(flow_poly, num_frames);
    
    segment_flow_confidence(y_min:y_max,x_min:x_max) = ...
        max(segment_flow_adjusted_confidence(y_min:y_max,x_min:x_max), flow_confidence);
    
    segment_flow_adjusted_confidence(y_min:y_max,x_min:x_max) = ...
        max(segment_flow_adjusted_confidence(y_min:y_max,x_min:x_max), adjusted_confidence);
    
end

if do_plot
    figure;
    subplot(1,2,1); imgray(segment_flow_confidence);
    subplot(1,2,2); imgray(segment_flow_adjusted_confidence);
end

    