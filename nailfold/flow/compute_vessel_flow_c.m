function [segment_flow] = ...
    compute_vessel_flow_c(segment_mask, ...
    flow_data_dir, flow_results_dir, flow_metrics_dir, vessels_list, do_plot)
%COMPUTE_VESSEL_FLOW *Insert a one line summary here*
%   [segment_flow.sum segment_flow.count segment_flow.max] =...
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

[frame_h frame_w] = size(segment_mask);

segment_flow.sum = zeros(frame_h, frame_w);
segment_flow.max = zeros(frame_h, frame_w);
segment_flow.count = zeros(frame_h, frame_w);
segment_flow.mask = zeros(frame_h, frame_w);
%
for i_ve = 1:length(vessels_list)  
    
    %Load in data defining the location of the vessel frames and the flow
    %results
    load([flow_data_dir vessels_list(i_ve).name],...
        'x_min_cropped', 'y_min_cropped', 'x_max_cropped', 'y_max_cropped');
    load([flow_results_dir vessels_list(i_ve).name], 'flow_results');
    load([flow_metrics_dir vessels_list(i_ve).name], 'flow_metrics');
    
    %Workout where the re-registered frames are extracted from in the
    cropped_rows = y_min_cropped:y_max_cropped;
    cropped_cols = x_min_cropped:x_max_cropped;
    
    valid_rows = cropped_rows >= 1 & cropped_rows <= frame_h; 
    valid_cols = cropped_cols >= 1 & cropped_cols <= frame_w; 
    
    %Project the flow results in to the segment
    segment_flow.sum(cropped_rows(valid_rows), cropped_cols(valid_cols)) = ...
        segment_flow.sum(cropped_rows(valid_rows), cropped_cols(valid_cols))...
        + flow_results.flowPyramidEst{1}(valid_rows,valid_cols);
    
    segment_flow.max(cropped_rows(valid_rows), cropped_cols(valid_cols)) = ...
        max(segment_flow.max(cropped_rows(valid_rows), cropped_cols(valid_cols)),...
        flow_results.flowPyramidEst{1}(valid_rows,valid_cols));
    
    segment_flow.count(cropped_rows(valid_rows), cropped_cols(valid_cols)) = ...
        segment_flow.count(cropped_rows(valid_rows), cropped_cols(valid_cols)) + 1;
    
    segment_flow.mask(cropped_rows(valid_rows), cropped_cols(valid_cols)) = ...
        segment_flow.mask(cropped_rows(valid_rows), cropped_cols(valid_cols)) |...
        flow_metrics.frame_apex_mask(valid_rows,valid_cols);
    
end

if do_plot > 1
    figure;
    subplot(1,2,1); imgray(segment_flow.mask);
    subplot(1,2,2); show_flow_as('rgb', segment_flow.max);
end

    