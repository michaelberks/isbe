function [apex_measures] = ...
    extract_apex_measures_newest(vessel_prob, vessel_ori, vessel_width, vessel_im, apex_xy, varargin)
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
    'max_dist', 120, ...
    'border_sz', 16, ...
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

num_apexes = size(apex_xy, 1);
apex_measures.width_at_apex = nan(num_apexes,1);
apex_measures.mean_width = nan(num_apexes,1);
apex_measures.median_width = nan(num_apexes,1);
apex_measures.max_width = nan(num_apexes,1);
apex_measures.min_width = nan(num_apexes,1);
apex_measures.std_width = nan(num_apexes,1);
apex_measures.mean_weighted_width = nan(num_apexes,1);
apex_measures.mean_weighted_prob = nan(num_apexes,1);
apex_measures.total_prob = nan(num_apexes,1);
apex_measures.orientation_hist = nan(num_apexes,args.num_ori_bins);
apex_measures.base_orientation = nan(num_apexes,1);
apex_measures.initial_connected_orientation = nan(num_apexes,1);
apex_measures.connected_orientation = nan(num_apexes,1);
apex_measures.weighted_orientation = nan(num_apexes,1);
apex_measures.bounding_box = nan(4,2,num_apexes);
apex_measures.apex_xy = apex_xy;

if args.plot
    figure; imgray(vessel_im); a1 = gca;
end

vessel_mask = vessel_prob > args.connect_thresh;
max_dist = args.max_dist;
border_size = args.border_sz;
[rows cols] = size(vessel_mask);

for i_pt = 1:num_apexes

    %Get predicted scale and orientation at this point
    ax = apex_xy(i_pt,1);
    ay = apex_xy(i_pt,2);
    a_ori = vessel_ori(round(ay), round(ax));
    a_theta = angle(a_ori) / 2;
    a_width = vessel_width(round(ay), round(ax));
    
    dist_mask_0 = false(rows, cols);
    dist_mask_0(ay, ax) = 1;
    dist_map = inf(rows, cols);
    dist_map(ay, ax) = 0;
    for i_d = 1:max_dist
        dist_mask_1 = imdilate(dist_mask_0, strel('disk', 1)) & vessel_mask;
        dist_map(dist_mask_1 & ~dist_mask_0) = i_d;
        dist_mask_0 = dist_mask_1;
    end
    [y x] = find(dist_mask_0);
    
    if isempty(x)
        apex_measures.apex_xy(i_pt,:) = nan;
        continue;
    end
    
    apex_measures.base_orientation(i_pt) = a_ori;    
    if isstruct(args.width_predictor)
        [predicted_width_ratio] = predict_apex_width_image(...
           vessel_im, [ax, ay], a_width, args.width_predictor,...
           'num_profile_pts', 50, 'scale_offsets', [0.5 1 2], 'theta', a_theta+pi/2);
        apex_measures.width_at_apex(i_pt) = a_width * predicted_width_ratio;
    else
        apex_measures.width_at_apex(i_pt) = a_width;
    end
    
    x_min = max(1, min(x) - border_size);
    y_min = max(1, min(y) - border_size);
    x_max = min(max(x) + border_size, cols);
    y_max = min(max(y) + border_size, rows);
    
    if args.plot
        plot(a1, [x_min x_max x_max x_min x_min], [y_min y_min y_max y_max y_min]);
    end
    
%     sz_x = x_max - x_min + 1;
%     sz_y = y_max - y_min + 1;
    
    x_lim = ax - x_min + 1;
    y_lim = ay - y_min + 1;

    %Sample vessel prob patch
    vessel_prob_patch = vessel_prob(y_min:y_max, x_min:x_max);% interp2(vessel_prob, xa, ya, '*linear', 0);
    vessel_ori_patch = vessel_ori(y_min:y_max, x_min:x_max);%interp2(vessel_ori, xa, ya, '*linear', 0);  
    vessel_width_patch = vessel_width(y_min:y_max, x_min:x_max);%interp2(vessel_width, xa, ya, '*linear', 0);
    connectivity_mask = dist_mask_0(y_min:y_max, x_min:x_max);%;
    vessel_prob_patch(~connectivity_mask) = 0;
    
    %Compute connectivity map and estimate vessel orientation
    [connectivity_map] = make_connectivity_map(vessel_prob_patch, x_lim, y_lim, args.num_c_pts);
    
    apex_measures.initial_connected_orientation(i_pt) = mean(vessel_ori_patch(:).*connectivity_map(:));
   
    %Compute the weighted mean width
    apex_measures.mean_weighted_width(i_pt) = ...
        sum(connectivity_map(:) .* vessel_width_patch(:)) / sum(connectivity_map(:));
    
    apex_measures.mean_weighted_prob(i_pt) = ...
        sum(connectivity_map(:) .* vessel_prob_patch(:)) / sum(connectivity_map(:));
    
    apex_measures.weighted_orientation(i_pt) = mean(vessel_ori_patch(:).*connectivity_map(:));
    
    apex_measures.bounding_box(:,:,i_pt) = [...
        x_min y_min
        x_min y_max
        x_max y_min
        x_max y_max];
    
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
    vessel_ori_patch = vessel_ori_patch * conj(-apex_measures.connected_orientation(i_pt));
    
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
    
    if args.plot > 1
        figure;
        subplot(2,3,1); imgray(vessel_prob_patch);
        subplot(2,3,2); imgray(complex2rgb(vessel_ori_patch));
        subplot(2,3,3); imgray(vessel_width_patch);
        title(num2str([apex_measures.mean_width(i_pt) apex_measures.std_width(i_pt)]));
        
        subplot(2,3,4); imgray(connectivity_map);
        subplot(2,3,5); imgray(connectivity_mask);
        subplot(2,3,6); bar(circshift(apex_measures.orientation_hist(i_pt,:), [0 args.num_ori_bins/2]));
        
%         centre_mask = bwmorph(connectivity_mask, 'thin', inf);
%         [cy cx] = find(centre_mask);
%         [vy vx] = find(connectivity_mask);
%         dist_transform = bwdist(~connectivity_mask);
%         cw = 2*double(dist_transform(centre_mask))-1;
%         v_widths = griddata(cx, cy, cw, vx, vy, 'nearest');
%         v_width_map = zeros(size(connectivity_mask));
%         v_width_map(connectivity_mask) = v_widths;
%         w_diff = v_width_map - vessel_width_patch;
%         w_diff(~connectivity_mask) = 0;
%         
%         figure;
%         subplot(1,2,1); imgray(v_width_map); colorbar;
%         subplot(1,2,2); imgray(w_diff); colormap(jet(256)); colorbar;       
        
    end     
    
end

%Compute the distance between each

    
 
    
    
