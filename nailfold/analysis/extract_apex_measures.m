function [apex_measures] = ...
    extract_apex_measures(vessel_prob, vessel_ori, vessel_width, vessel_im, apex_xy, varargin)
%EXTRACT_VESSEL_CENTRES *Insert a one line summary here*
%   [] = extract_vessel_centres()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 18-Jun-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'prob_sigma',           2,...
    'ori_sigma',            2,...
    'width_sigma',            2,...
    'num_c_pts', 20,...
    'xy', [],...
    'patch_sz', [], ...
    'num_ori_bins', 36,...
    'connect_thresh', 0.5,...
    'base_width', 20, ...
    'width_predictor', [],...
    'plot', 0);
    
%Smooth the vessel probs
if args.prob_sigma
    g_prob = gaussian_filters_1d(args.prob_sigma);
    g_prob = g_prob / sum(g_prob);    
    vessel_prob = conv2(g_prob', g_prob, vessel_prob, 'same');
end
if args.ori_sigma
    g_ori = gaussian_filters_1d(args.ori_sigma);
    g_ori = g_ori / sum(g_ori);
    vessel_ori = conv2(g_ori', g_ori, vessel_ori, 'same');
end 
if args.width_sigma
    g_width = gaussian_filters_1d(args.width_sigma);
    g_width = g_width / sum(g_width);
    vessel_width = conv2(g_width', g_width, vessel_width, 'same');
end 

if isempty(args.xy)    
    %Set up x,y coordinates for template patch
    x_lim = (args.patch_sz(2)-1)/2;
    y_lim = (args.patch_sz(1)-1)/2;
    x = repmat(linspace(-x_lim, x_lim, args.patch_sz(2)), args.patch_sz(1), 1);
    y = repmat(linspace(-y_lim, y_lim, args.patch_sz(1))', 1, args.patch_sz(2));
    args.xy = [x(:) y(:)];
else
    x_lim = 1-min(args.xy(:,1));
    y_lim = 1-min(args.xy(:,2));
end

num_apexes = size(apex_xy, 1);
apex_measures.width_at_apex = zeros(num_apexes,1);
apex_measures.mean_width = zeros(num_apexes,1);
apex_measures.median_width = zeros(num_apexes,1);
apex_measures.max_width = zeros(num_apexes,1);
apex_measures.min_width = zeros(num_apexes,1);
apex_measures.std_width = zeros(num_apexes,1);
apex_measures.mean_weighted_width = zeros(num_apexes,1);
apex_measures.mean_weighted_prob = zeros(num_apexes,1);
apex_measures.total_prob = zeros(num_apexes,1);
apex_measures.orientation_hist = zeros(num_apexes,args.num_ori_bins);
apex_measures.base_orientation = zeros(num_apexes,1);
apex_measures.connected_orientation = zeros(num_apexes,1);
apex_measures.weighted_orientation = zeros(num_apexes,1);
apex_measures.apex_xy = apex_xy;

for i_pt = 1:num_apexes

    %Get predicted scale and orientation at this point
    ax = apex_xy(i_pt,1);
    ay = apex_xy(i_pt,2);
    a_ori = vessel_ori(round(ay), round(ax));
    a_theta = angle(a_ori) / 2;
    a_width = vessel_width(round(ay), round(ax));
    
    apex_measures.base_orientation(i_pt) = a_ori;    
    
    if isstruct(args.width_predictor)
        [predicted_width_ratio] = predict_apex_width_image(...
           vessel_im, [ax, ay], a_width, args.width_predictor,...
           'num_profile_pts', 50, 'scale_offsets', [0.5 1 2], 'theta', a_theta+pi/2);
        apex_measures.width_at_apex(i_pt) = a_width * predicted_width_ratio;
    else
        apex_measures.width_at_apex(i_pt) = a_width;
    end
    
    %Get scale relative to base width a make rotation matrix
    rot = [cos(a_theta) -sin(a_theta); sin(a_theta) cos(a_theta)];
    scale = a_width / args.base_width;

    %Transform points given scale and angle and translate to
    %candidate position
    xya = args.xy * rot * scale;
    xa = reshape(xya(:,1) + ax, args.patch_sz);
    ya = reshape(xya(:,2) + ay, args.patch_sz);

    %Sample vessel prob patch
    vessel_prob_patch = interp2(vessel_prob, xa, ya, '*linear', 0);
    vessel_ori_patch = interp2(vessel_ori, xa, ya, '*linear', 0);
    vessel_width_patch = interp2(vessel_width, xa, ya, '*linear', 0);
    
    
    [connectivity_map] = make_connectivity_map(vessel_prob_patch, x_lim, y_lim, args.num_c_pts);
    connectivity_mask = connectivity_map >= args.connect_thresh;
    
    %Compute the weighted mean width
    apex_measures.mean_weighted_width(i_pt) = ...
        sum(connectivity_map(:) .* vessel_width_patch(:)) / sum(connectivity_map(:));
    
    apex_measures.mean_weighted_prob(i_pt) = ...
        sum(connectivity_map(:) .* vessel_prob_patch(:)) / sum(connectivity_map(:));
    
    apex_measures.weighted_orientation(i_pt) = mean(vessel_ori_patch(:).*connectivity_map(:));
    
    %Now compute measure based on a fixed connectivity threshold
    if any(connectivity_mask(:))
        apex_measures.mean_width(i_pt) = mean(vessel_width_patch(connectivity_mask));
        apex_measures.median_width(i_pt) = median(vessel_width_patch(connectivity_mask));
        apex_measures.max_width(i_pt) = max(vessel_width_patch(connectivity_mask));
        apex_measures.min_width(i_pt) = min(vessel_width_patch(connectivity_mask));
        apex_measures.std_width(i_pt) = std(vessel_width_patch(connectivity_mask));
        apex_measures.total_prob(i_pt) = sum(vessel_prob_patch(connectivity_mask));
        apex_measures.connected_orientation(i_pt) = mean(vessel_ori_patch(connectivity_mask));
    end 
    
    %Create the orientation histogram indices
    
    %Derotate the ori patch by the base orientation - make sure we do this
    %after we've saved the base, connected and weighted orientation
    %measures
    vessel_ori_patch = vessel_ori_patch * conj(a_ori);
    
    %Compute angle on half-circle (will be [-pi/2 pi/2]) then scale to [0 1]
    g_theta_idx = angle(vessel_ori_patch)/(2*pi) + 0.5;
    g_theta_idx = g_theta_idx * args.num_ori_bins;
    
    %Convert theta vals in into integer indices - there may still be some 0s,
    %which should be moved to the final bin
    g_theta_idx_c = ceil(g_theta_idx);
    g_theta_w = g_theta_idx_c - g_theta_idx;
    
    g_theta_idx_c(~g_theta_idx_c) = args.num_ori_bins;   
    g_theta_idx_f = g_theta_idx_c - 1;
    g_theta_idx_f(~g_theta_idx_f) = args.num_ori_bins;
    
    %Now compute the orientation histograms weighted by the connectivity
    %map
    apex_measures.orientation_hist(i_pt,:) = full(...
        sparse(1,g_theta_idx_c,(1-g_theta_w).*connectivity_map,1,args.num_ori_bins) + ...
        sparse(1,g_theta_idx_f,g_theta_w.*connectivity_map,1,args.num_ori_bins));   
    
    if args.plot
        figure;
        subplot(2,3,1); imgray(vessel_prob_patch);
        subplot(2,3,2); imgray(complex2rgb(vessel_ori_patch));
        subplot(2,3,3); imgray(vessel_width_patch);
        title(num2str([apex_measures.mean_width(i_pt) apex_measures.std_width(i_pt)]));
        
        subplot(2,3,4); imgray(connectivity_map);
        subplot(2,3,5); imgray(connectivity_mask);
        subplot(2,3,6); bar(apex_measures.orientation_hist(i_pt,:));
    end     
    
end

%Compute the distance between each

    
 
    
    
