function [potential_map scale_map, rotation_map, mean_width_map, mean_ori_map] = ...
    estimate_vessel_pose2(vessel_prob_map, width_map, orientation_map, varargin)
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
    'patch_size', 65,...
    'smooth_vessel_probs', 2,...
    'any_vessel_thresh', 0.25,...
    'strong_vessel_thresh', 0.25,...
    'width_m', 1,...
    'width_c', 0,...
    'vessel_ori_offset', pi/2,...
    'plot', 0);
clear varargin;

[rows cols] = size(vessel_prob_map);

%Smooth the vessel probs
if args.smooth_vessel_probs
    g = gaussian_filters_1d(args.smooth_vessel_probs);
    g = g / sum(g);
    vessel_prob_map = conv2(g', g, vessel_prob_map, 'same');
end

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

%Pre-allocate output
mean_width_map = zeros(rows, cols);
mean_ori_map = zeros(rows, cols);

%Get map of connect components
connected_label_map = bwlabel(potential_map);
for i_la = 1:max(connected_label_map(:))
    label_mask = connected_label_map == i_la;
    mean_width_map(label_mask) = naNmean(width_map(label_mask));
    mean_ori_map(label_mask) = naNmean(orientation_map(label_mask));
end

%Remove the padding from the potential map
scale_map = args.width_m * mean_width_map + args.width_c;
rotation_map = mod(angle(mean_ori_map)/2,pi) - args.vessel_ori_offset;
    
    
    
    
    
