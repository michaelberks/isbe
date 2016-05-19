function [flow_metrics] = ...
    compute_vessel_flow_rf(varargin)
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
args = u_packargs(varargin, '0', ...
    'vessel_name',              [],...
    'frames',                   [],...  
    'flow_results',             [],...
    'rescale_factor',           0,...
    'flow_data_dir',            'N:\Nailfold Capillaroscopy\wellcome\flow_data\',...
    'flow_results_dir',         'N:\Nailfold Capillaroscopy\wellcome\flow_results\',...
    'model_root',               'C:\isbe\nailfold\models\',...
    'rf_ori',                   'vessel/orientation/rf_regression/296621/predictor.mat',...
    'rf_vessels',               'vessel/detection/rf_classification/296655/predictor.mat',...
    'rf_width',                 'vessel/width/rf_regression/297037/predictor.mat',...
    'rf_vessels_args',          'vessel/detection/rf_classification/296655/job_args.mat',...
    'vessel_mask_thresh',       0.25,...
    'max_dist',                 100,...
    'vessel_flow',              [],...
    'patch_contrast_scale',     70,...
    'vessels_sigma',            2,...
    'hog_class_forest_path',    'apex/classification/set12g_half_296655/rf.mat',...
    'hog_off_x_forest_path',    'apex/offset_x/set12g_half_296655/rf.mat',...
    'hog_off_y_forest_path',    'apex/offset_y/set12g_half_296655/rf.mat',...
    'separate_trees',           0,...
    'apex_class_thresh',        0.1,...
    'base_width',               20,...
    'num_cells',                8,...
    'cell_sz',                  8,... %Size of HoG cells in blocks
    'block_sz',                 [2 2],...%Size of blocks in cells
    'hog_ori_bins',             9,... %Number of bins in orientation histograms
    'norm_method',              'none',... %Method for local normalisation 'none'
    'block_spacing',            8,... %Separation between blocks in pixels - controls level of overlap
    'gradient_operator',        [-1 0 1],...
    'spatial_sigma',            0, ...
    'angle_wrap',               1,...
    'hist_ori_bins',            36,...
    'hist_flow_bins',           0:0.2:10,...
    'plot', 0);

%Load frames and flow predictions
if isempty(args.flow_results)
    load([args.flow_results_dir args.vessel_name], 'flow_results');
    patch_flow = flow_results.flowPyramidEst{1}; 
else
    patch_flow = args.flow_results;
end

[rows cols] = size(patch_flow);
    
%Re-predict vessel probs
if isempty(args.frames)
    load([args.flow_data_dir args.vessel_name],...
        'cropped_frames');
    patch = args.patch_contrast_scale*mean(cropped_frames,3);
else
    patch = args.patch_contrast_scale*mean(args.frames,3);
end
[p_rows p_cols] = size(patch);

if args.rescale_factor    
    patch = imresize(patch, 1 / args.rescale_factor, 'lanczos2');    
end

%Load forests
rfs = cell(3,1);
rfs{1} = u_load([args.model_root args.rf_vessels]);
rfs{2} = u_load([args.model_root args.rf_ori]);
rfs{3} = u_load([args.model_root args.rf_width]);
rf_args = u_load([args.model_root args.rf_vessels_args]);
    
vessel_predictions = predict_image(... % non-strict mode
    'image_in', patch,...
    'decomposition_args', rf_args.decomposition_args,...
    'predictor', rfs, ...
    'prediction_type', {'rf_classification', 'rf_regression', 'rf_regression'},...
    'output_type', {'detection', 'orientation', 'width'},...
    'use_probs', 0,...
    'mask', [],...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 1024,...
    'incremental_results', 0);   

patch_ves = vessel_predictions(:,:,1);
patch_ori = vessel_predictions(:,:,2);
patch_wid = vessel_predictions(:,:,3);

if args.rescale_factor
    patch_ves = imresize(patch_ves, [p_rows p_cols], 'lanczos2');
    patch_ori = imresize(patch_ori, [p_rows p_cols], 'lanczos2');
    patch_wid = args.rescale_factor*imresize(patch_wid, [p_rows p_cols], 'lanczos2');
end

patch_ves = mb_pad(patch_ves, [rows cols] - [p_rows p_cols], 'replicate', 'post');
patch_ori = mb_pad(patch_ori, [rows cols] - [p_rows p_cols], 'replicate', 'post');
patch_wid = mb_pad(patch_wid, [rows cols] - [p_rows p_cols], 'replicate', 'post');

%Smooth to get a continous ridge at vessel centre
g_prob = gaussian_filters_1d(args.vessels_sigma);
g_prob = g_prob / sum(g_prob);
smooth_ves = conv2(g_prob', g_prob, patch_ves, 'same');
smooth_ori = conv2(g_prob', g_prob, patch_ori, 'same');
smooth_wid = conv2(g_prob', g_prob, patch_wid, 'same');

%make a vessel mask connected to the apex
frame_vessels_mask = smooth_ves > args.vessel_mask_thresh;

%Get the connected centre mask by thinning to infinity, then thicken by 1
%dilation
thin_mask = bwmorph(frame_vessels_mask, 'thin', 'inf');
background_mask = ~frame_vessels_mask;

%Get coordinates for thin and medium mask
[cy cx] = find(thin_mask);

%Re-predict apex location
apex_class_rf = u_load([args.model_root args.hog_class_forest_path]);
apex_offset_x_rf = u_load([args.model_root args.hog_off_x_forest_path]);
apex_offset_y_rf = u_load([args.model_root args.hog_off_y_forest_path]);

% Get HoG args from main args
hog_args.cell_sz = [args.cell_sz args.cell_sz];
hog_args.block_sz = args.block_sz;
hog_args.num_ori_bins = args.hog_ori_bins;
hog_args.norm_method = args.norm_method;
hog_args.block_spacing = args.block_spacing;
hog_args.gradient_operator = args.gradient_operator;
hog_args.spatial_sigma = args.spatial_sigma;
hog_args.angle_wrap = args.angle_wrap;

% Get patch size and form template x,y coordinates for the patch
patch_sz = args.num_cells*args.cell_sz;
patch_sz = patch_sz + 2; %Account for padding
patch_sz2 = (patch_sz - 1)/2;

% Get hog size from the hog_args struct
% Set up x,y coordinates for template patch
hx = repmat(-patch_sz2:patch_sz2, patch_sz, 1);
hy = hx';
hxy = [hx(:) hy(:)];

vessel_centre.x = cx;
vessel_centre.y = cy;
vessel_centre.ori = smooth_ori(thin_mask);
vessel_centre.width = smooth_wid(thin_mask);

[apex_offset_map] = ...
    predict_apex_offsets(...
        'apex_class_rf', apex_class_rf,...
        'apex_offset_x_rf', apex_offset_x_rf,...
        'apex_offset_y_rf', apex_offset_y_rf,...
        'vessel_feature_im', smooth_ves, ...
        'vessel_centre', vessel_centre, ...
        'separate_trees', args.separate_trees,...
        'smoothing_sigma', 0,...
        'num_cells', args.num_cells,...
        'hog_args', hog_args,...
        'xy', hxy,...
        'apex_class_thresh', args.apex_class_thresh,...
        'base_width', args.base_width,...
        'include_pts', []);

[~, max_idx] = max(apex_offset_map(:));
[maxr maxc] = ind2sub([rows cols], max_idx);
apex_xy = [maxc maxr];
    
%Update the vessel mask given the new apex
vessel_mask_new = bwselect(frame_vessels_mask, apex_xy(1), apex_xy(2));
if any(vessel_mask_new(:) ~= frame_vessels_mask(:))
    frame_vessels_mask = vessel_mask_new;
    thin_mask = thin_mask & frame_vessels_mask;
    [cy cx] = find(thin_mask);
end    

%Find point on thin maks closest to apex and set this region by 0 distance
dists = (cx - apex_xy(1)).^2 + (cy - apex_xy(:,2)).^2;
[~, min_i] = min(dists);
x1 = cx(min_i);
y1 = cy(min_i);
apex_mask = false(rows, cols);
apex_mask(y1, x1) = 1;
apex_mask = imdilate(apex_mask, strel('disk', 2));

%Make a distance map for each point on the medium mask to the vessel apex
dist_mask_0 = false(rows, cols);
dist_mask_0(apex_mask) = 1;
dist_map = inf(rows, cols);
dist_map(apex_mask) = 0;
for i_d = 1:args.max_dist
    dist_mask_1 = imdilate(dist_mask_0, strel('disk', 1)) & frame_vessels_mask;
    dist_map(dist_mask_1 & ~dist_mask_0) = i_d;
    if all(dist_mask_0(:) == dist_mask_1(:))
        break;
    else
        dist_mask_0 = dist_mask_1;
    end
end
frame_apex_mask = dist_mask_1;

%Compute average  flow rates inside and outside the predicted vessel.
%If these are the same flow probably hasn't worked
flow_velocity = abs(patch_flow);
flow_metrics.vessel_flow = mean(flow_velocity(frame_apex_mask));
flow_metrics.background_flow = mean(flow_velocity(background_mask));

%Now compute weighted sums of flow rates and orientation errors for the 
%flow predictions    
frame_pred_masked = smooth_ves;
frame_pred_masked(~frame_apex_mask) = 0;
flow_metrics.total_vessel_prob = sum(frame_pred_masked(:));

flow_angle = exp(2i*angle(patch_flow));
flow_angle_diffs = frame_pred_masked .* smooth_ori .* flow_angle;
flow_angle_weights = abs(flow_angle_diffs);
flow_angle_dirs = angle(flow_angle_diffs);

flow_metrics.weighted_flow_rate = sum(frame_pred_masked(:).*flow_velocity(:)) /...
    flow_metrics.total_vessel_prob;

flow_metrics.mean_error = mean(abs(flow_angle_dirs(frame_apex_mask)));
flow_metrics.mean_weighted_error = sum(flow_angle_weights(:) .* abs(flow_angle_dirs(:))) /...
    sum(flow_angle_weights(:));

%Compute a histogram of flow direction errors
g_theta_idx = flow_angle_dirs/(2*pi) + 0.5;
g_theta_idx = g_theta_idx * args.hist_ori_bins;

%Convert theta vals in into integer indices - there may still be some 0s,
%which should be moved to the final bin
g_theta_idx_c = ceil(g_theta_idx);
g_theta_w = g_theta_idx_c - g_theta_idx;

g_theta_idx_c(~g_theta_idx_c) = args.hist_ori_bins;   
g_theta_idx_f = g_theta_idx_c - 1;
g_theta_idx_f(~g_theta_idx_f) = args.hist_ori_bins;

flow_metrics.flow_errors_hist = full(...
    sparse(1, g_theta_idx_c, (1-g_theta_w).*flow_angle_weights, 1, args.hist_ori_bins) + ...
    sparse(1, g_theta_idx_f, g_theta_w.*flow_angle_weights, 1, args.hist_ori_bins));

flow_metrics.flow_hist = hist(abs(patch_flow(frame_apex_mask)), args.hist_flow_bins);

% what was the width of the associated vessel
flow_metrics.mean_width = mean(smooth_wid(frame_apex_mask));

% shape score
flow_metrics.shape_score = abs(mean(smooth_ori(frame_vessels_mask)));

if args.plot
    figure;
    subplot(2,3,1); imgray(patch); plot(apex_xy(1), apex_xy(2), 'rx');
    title(['ME = ' num2str(flow_metrics.mean_error,3) ' MWE = ' num2str(flow_metrics.mean_weighted_error, 3)]);
    subplot(2,3,2); imgray(patch_ves); plot(apex_xy(1), apex_xy(2), 'rx');
    title(['MFR = ' num2str(flow_metrics.weighted_flow_rate, 3)]);        
    subplot(2,3,3); bar(flow_metrics.flow_hist);
    title(['Flow ratio = ' num2str(flow_metrics.vessel_flow / flow_metrics.background_flow, 3)]);
    
    subplot(2,3,4); imgray(complex2rgb(smooth_ori, [], [], [], 1));
    subplot(2,3,5); imgray(complex2rgb(conj(patch_flow.^2), [], [], [], 1));
    subplot(2,3,6); imgray(apex_offset_map);    
    
end

flow_metrics.vessel_pred = uint8(100*patch_ves);
flow_metrics.vessel_ori = uint8(255*complex2rgb(patch_ori,[],1));
flow_metrics.vessel_width = uint8(patch_wid);
flow_metrics.apex_xy = apex_xy;
flow_metrics.frame_apex_mask = frame_apex_mask;

    