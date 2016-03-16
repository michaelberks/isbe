function [segment_flow] = ...
    compute_vessel_flow_c(segment_mask, vessel_pred, vessel_ori, vessel_data,...
    data_dir, results_dir, flow_metrics_dir, vessels_list, do_plot)
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

connect_thresh = 0.2;
num_ori_bins = 36;

[frame_h frame_w] = size(segment_mask);

segment_flow.sum = zeros(frame_h, frame_w);
segment_flow.max = zeros(frame_h, frame_w);
segment_flow.count = zeros(frame_h, frame_w);
segment_flow.mask = zeros(frame_h, frame_w);
%
for i_ve = 1:length(vessels_list)  
    
    %Load in data defining the location of the vessel frames and the flow
    %results
    load([data_dir vessels_list(i_ve).name],...
        'x_min', 'x_max', 'y_min', 'y_max', 'edge_mask',...
        'vessel_transforms', 'apex_idx', 'frames_offset', 'cropped_frames');
    load([results_dir vessels_list(i_ve).name], 'flow_results');
    
    %Workout where the re-registered frames are extracted from in the
    %segment
    [~, ~, dx, dy] = ...
        mosaic_limits(size(flow_results.flowPyramidEst{1}), vessel_transforms);
    
    x_start = floor(x_min + find(any(~edge_mask,1),1) - (dx + 1));
    y_start = floor(y_min + find(any(~edge_mask,2),1) - (dy + 1));
    
    [flow_h flow_w] = size(flow_results.flowPyramidEst{1});
    cropped_rows = y_start+(1:flow_h);
    cropped_cols = x_start+(1:flow_w);
    
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
    
    %Now compute information about these frames - first get vessel and
    %orientation patches for the correctly translated bundign box
    pred_rows = cropped_rows+floor(frames_offset(2));
    pred_cols = cropped_cols+floor(frames_offset(1));
    valid_pred_rows = pred_rows >= 1 & pred_rows <= size(vessel_ori,1);
    valid_pred_cols = pred_cols >= 1 & pred_cols <= size(vessel_ori,2);
    
    frame_ori = zeros(flow_h, flow_w, 3);
    frame_ori(valid_pred_rows, valid_pred_cols, :) = vessel_ori(...
        pred_rows(valid_pred_rows),...
        pred_cols(valid_pred_cols),:);
    frame_ori = rgb2complex(frame_ori, [], 1, [], 0);
    frame_pred = zeros(flow_h, flow_w);
    frame_pred(valid_pred_rows, valid_pred_cols) = double(vessel_pred(...
        pred_rows(valid_pred_rows),...
        pred_cols(valid_pred_cols))) / 100;
    
    %Make a frame mask from the vessel predictions, selected pixels
    %connected to the vessel apex
    ax = vessel_data.apex_measures.distal.apex_xy(apex_idx,1)/vessel_data.resize_factor ...
        - (floor(frames_offset(1)) + x_start);
    ay = vessel_data.apex_measures.distal.apex_xy(apex_idx,2)/vessel_data.resize_factor ...
        - (floor(frames_offset(2)) + y_start);
    
    apex_xy = [ax ay]; %#ok
    save([results_dir vessels_list(i_ve).name], 'apex_xy', '-append');
    continue;
    
    frame_vessels_mask = frame_pred > connect_thresh;
    frame_apex_mask = bwselect(frame_vessels_mask, ax, ay, 4);
    
    %Project the frame mask into the segment mask
    segment_flow.mask(cropped_rows(valid_rows), cropped_cols(valid_cols)) = ...
        segment_flow.mask(cropped_rows(valid_rows), cropped_cols(valid_cols)) |...
        frame_vessels_mask(valid_rows,valid_cols); 
    
    %Compute average  flow rates inside and outside the predicted vessel.
    %If these are the same flow probably hasn't worked
    flow_velocity = abs(flow_results.flowPyramidEst{1});
    vessel_flow = mean(flow_velocity(frame_vessels_mask));
    background_flow = mean(flow_velocity(~frame_vessels_mask));
    
    %Now compute weighted sums of flow rates and orientation errors for the 
    %flow predictions    
    frame_pred_masked = frame_pred;
    frame_pred_masked(~frame_apex_mask) = 0;
    total_vessel_prob = sum(frame_pred_masked(:));

    flow_angle = exp(2i*angle(flow_results.flowPyramidEst{1}));
    flow_angle_diffs = frame_pred_masked .* frame_ori .* flow_angle;
    flow_angle_weights = abs(flow_angle_diffs);
    flow_angle_dirs = angle(flow_angle_diffs);

    weighted_flow_rate = sum(frame_pred_masked(:).*flow_velocity(:)) /...
        total_vessel_prob;

    mean_error = mean(abs(flow_angle_dirs(:)));
    mean_weighted_error = sum(flow_angle_weights(:) .* abs(flow_angle_dirs(:))) /...
        sum(flow_angle_weights(:));

    %Compute a histogram of flow direction errors
    g_theta_idx = flow_angle_dirs/(2*pi) + 0.5;
    g_theta_idx = g_theta_idx * num_ori_bins;

    %Convert theta vals in into integer indices - there may still be some 0s,
    %which should be moved to the final bin
    g_theta_idx_c = ceil(g_theta_idx);
    g_theta_w = g_theta_idx_c - g_theta_idx;

    g_theta_idx_c(~g_theta_idx_c) = num_ori_bins;   
    g_theta_idx_f = g_theta_idx_c - 1;
    g_theta_idx_f(~g_theta_idx_f) = num_ori_bins;

    flow_errors_hist = full(...
        sparse(1, g_theta_idx_c, (1-g_theta_w).*flow_angle_weights, 1, num_ori_bins) + ...
        sparse(1, g_theta_idx_f, g_theta_w.*flow_angle_weights, 1, num_ori_bins));

    % what was the width of the associated vessel
    weighted_width = vessel_data.apex_measures.distal.mean_weighted_width(apex_idx); %#ok
    
    if i_ve < do_plot
        figure;
        subplot(2,3,1); imgray(complex2rgb(frame_ori, [], [], [], 1));
        subplot(2,3,2); imgray(complex2rgb(conj(flow_results.flowPyramidEst{1}.^2), [], [], [], 1));
        title(['MFR = ' num2str(weighted_flow_rate, 3)]);
        subplot(2,3,3); imgray(complex2rgb(flow_angle_diffs, [], [], [], 1));
        title(['Flow ratio = ' num2str(vessel_flow / background_flow, 3)]);
        subplot(2,3,4); bar(flow_errors_hist);
        title(['ME = ' num2str(mean_error,3) ' MWE = ' num2str(mean_weighted_error, 3)]);
        subplot(2,3,5); imgray(mean(cropped_frames,3)); plot(ax, ay, 'rx');
        subplot(2,3,6); imgray(frame_pred_masked); plot(ax, ay, 'rx');
    end

    save([flow_metrics_dir vessels_list(i_ve).name],...
        'mean_error', 'mean_weighted_error', 'weighted_flow_rate',...
        'flow_errors_hist', 'total_vessel_prob', 'weighted_width',...
        'vessel_flow', 'background_flow');
            
end

if do_plot > 1
    figure;
    subplot(1,2,1); show_flow_as('rgb', segment_flow.sum);
    subplot(1,2,2); show_flow_as('rgb', segment_flow.max);
end

    