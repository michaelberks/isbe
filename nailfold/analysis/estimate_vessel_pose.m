function [potential_map scale_map, rotation_map, mean_width_map, mean_ori_map, mean_curv_map] = ...
    estimate_vessel_pose(vessel_prob_map, width_map, orientation_map, varargin)
%ESTIMATE_VESSEL_POSE *Insert a one line summary here*
%   [scale_map, rotation_map] = estimate_vessel_pose(potential_map, width_map, orientation_map)
%
% Inputs:
%      vessel_prob_map - *Insert description of input variable here*
%
%      width_map - *Insert description of input variable here*
%
%      orientation_map - *Insert description of input variable here*
%
%
% Outputs:
%      scale_map - *Insert description of input variable here*
%
%      rotation_map - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 01-May-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
args = u_packargs(varargin, '0', ...
    'patch_size', 101,...
    'smooth_vessel_probs', 2,...
    'any_vessel_thresh', 0.25,...
    'strong_vessel_thresh', 0.25,...
    'width_m', 1,...
    'width_c', 0,...
    'vessel_ori_offset', pi/2,...
    'plot', 0);
clear varargin;

[rows cols] = size(vessel_prob_map);

%prior to smoothing, compute a BW mask of all potential vessel pts
all_vessels_map = vessel_prob_map > args.any_vessel_thresh;

%Smooth the vessel probs
if args.smooth_vessel_probs
    g = gaussian_filters_1d(args.smooth_vessel_probs);
    g = g / sum(g);
    vessel_prob_map = conv2(g', g, vessel_prob_map, 'same');
end

%Compute curvature map from the orientation map
curvature_map = abs(complex_gaussian_curvature(orientation_map, 2));

%Non-maximally supress vessel probs to locate potential centre lines
vessel_nms = mb_non_maximal_supp(vessel_prob_map, angle(orientation_map)/2);

%Perform hysteris on the centrelines
strong_vessels = vessel_nms > args.strong_vessel_thresh;
if any(strong_vessels(:))
    [rstrong cstrong] = find(strong_vessels);
    potential_map = bwselect(vessel_nms > 0, cstrong, rstrong, 8);
else
    potential_map = strong_vessels;
end

%Make sure the all vessels map is 1 where the potential points are (this
%might not happen because we took the threshold before smoothing the vessel
%probs)
all_vessels_map(potential_map) = 1;

%Get template size and set up sampling pts matrix
pad_sz = (args.patch_size-1)/2;
patch_idx = -pad_sz:pad_sz;

%Pad the maps to avoid boundary condition
potential_map = padarray(potential_map, [pad_sz pad_sz], 0);
width_map = padarray(width_map, [pad_sz pad_sz], 0);
orientation_map = padarray(orientation_map, [pad_sz pad_sz], 0);
curvature_map = padarray(curvature_map, [pad_sz pad_sz], 0);
all_vessels_map = padarray(all_vessels_map, [pad_sz pad_sz], 0);

%Get indices and subscripts of potential points
potential_idx = find(potential_map);
[potential_r potential_c] = ind2sub(size(potential_map), potential_idx);

%Pre-allocate output
mean_width_map = zeros(rows, cols);
mean_ori_map = zeros(rows, cols);
mean_curv_map = zeros(rows, cols);

%Get map of connect components
connected_label_map = bwlabel(all_vessels_map);

for i_pt = 1:length(potential_idx)
    r_i = potential_r(i_pt);
    c_i = potential_c(i_pt);
    
    %Get patches - guaranteed to be in bounds because of padding
    %all_vessels_patch = all_vessels_map(r_i + patch_idx, c_i + patch_idx);
    connected_label_patch = connected_label_map(r_i + patch_idx, c_i + patch_idx);
    width_patch = width_map(r_i + patch_idx, c_i + patch_idx);
    orientation_patch = orientation_map(r_i + patch_idx, c_i + patch_idx);
    curvature_patch = curvature_map(r_i + patch_idx, c_i + patch_idx);
    
    %Workout which pixels are connected to the centre
    %all_vessel_patch = bwselect(all_vessels_patch, pad_sz+1, pad_sz+1);
    all_vessels_patch = connected_label_patch == connected_label_map(r_i,c_i);
    
    mean_width_map(r_i-pad_sz, c_i-pad_sz) = mean(width_patch(all_vessels_patch));
    mean_ori_map(r_i-pad_sz, c_i-pad_sz) = mean(orientation_patch(all_vessels_patch));
    mean_curv_map(r_i-pad_sz, c_i-pad_sz) = mean(curvature_patch(all_vessels_patch));
    
end

%Remove the padding from the potential map
potential_map = potential_map(pad_sz+(1:rows), pad_sz+(1:cols));
scale_map = args.width_m * mean_width_map + args.width_c;
rotation_map = mod(angle(mean_ori_map)/2,pi) - args.vessel_ori_offset;
    
    
    
    
    
