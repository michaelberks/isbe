function [shape_contexts] = compute_shape_contexts(feature_map, orientation_map, selection_map, r_idx, c_idx, varargin)
%SHAPE_CONTEXT *Insert a one line summary here*
%   [] = shape_context()
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
% Created: 07-Jun-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

args = u_packargs(varargin, 0,...
    'weighting_map', [],...
    'num_rdist', 5,...
    'num_theta', 12, ...
    'log_r', 1.5,...
    'min_rdist', 10);
clear varargin;

%Apply weighting to features if set
if ~isempty(args.weighting_map)    
    feature_map = feature_map.*args.weighting_map;
end

num_pts = length(r_idx);
shape_contexts = zeros(args.num_rdist, args.num_theta, num_pts);

%Get labelling of connections
connected_label_map = bwlabel(selection_map);

%Make radial masks
[~, radial_idx] = make_log_radial_masks(args.num_rdist, args.min_rdist, args.log_r, 0);
max_radial_mask = radial_idx > 0;
%Make map of angles in mask
mask_sz = size(radial_idx,1);
max_rdist = floor(mask_sz/2);
xy = repmat(-max_rdist:max_rdist, mask_sz, 1);
theta_map = atan2(-xy', xy);

for i_pt = 1:num_pts
    r_i = r_idx(i_pt);
    c_i = c_idx(i_pt);
    
    %Make shape context masks oriented to point's direction
    theta0 = orientation_map(r_i, c_i);
    theta_idx = ceil( .5*args.num_theta*mod(theta_map-theta0, 2*pi)/pi );
    
    %Get feature and selection patch about local point
    feature_patch = sample_window(feature_map, mask_sz, r_i, c_i, 0);
    connected_label_patch = sample_window(connected_label_map, mask_sz, r_i, c_i, false);
    
    selected_patch = (connected_label_patch == connected_label_map(r_i,c_i)) ...
        & max_radial_mask;
    
    %Get indices to radial and theta bins for each selected pixel
    radial_idx_i = radial_idx(selected_patch);
    theta_idx_i = theta_idx(selected_patch);
    
    %Use sparse indexing trick to compute histogram sum in each bin
    shape_contexts(:,:,i_pt) = ...
            full( sparse(radial_idx_i, theta_idx_i, feature_patch(selected_patch), args.num_rdist, args.num_theta) ) ...
            .* exp(-2*theta0);
        
end
            
            
    
    





