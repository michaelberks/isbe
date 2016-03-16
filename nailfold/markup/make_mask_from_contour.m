function [vessel_mask vessel_centre_mask width_map ori_map ori_map_3d] =...
    make_mask_from_contour(outer_edge, inner_edge, rows, cols, apex_idx)
%MAKE_MASK_FROM_CONTOUR *Insert a one line summary here*
%   [vessel_mask] = make_mask_from_contour(outer_edge, inner_edge, rows, cols)
%
% Inputs:
%      outer_edge, inner_edge - *Insert description of input variable here*
%
%      rows - *Insert description of input variable here*
%
%      cols - *Insert description of input variable here*
%
%
% Outputs:
%      vessel_mask - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 16-Apr-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
      
num_pts = size(outer_edge,1);
apex_idx = sort(apex_idx);

num_layers = length(apex_idx) + 1;
curr_layer = 1;

vessel_mask_layers = false(rows, cols, num_layers);
vessel_mask_i = false(rows, cols);
for i_pt = 1:num_pts-1
    segment_xy = [...
        outer_edge(i_pt,:); inner_edge(i_pt,:);...
        inner_edge(i_pt+1,:); outer_edge(i_pt+1,:)];
    segment_mask = poly2mask(segment_xy(:,1), segment_xy(:,2),...
        rows, cols);
    vessel_mask_i = vessel_mask_i | segment_mask;
    
    if curr_layer < num_layers && i_pt == apex_idx(curr_layer)
        vessel_mask_layers(:,:,curr_layer) = vessel_mask_i;
        vessel_mask_i = false(rows, cols);
        curr_layer = curr_layer + 1;
    end       
end
vessel_mask_layers(:,:,num_layers) = vessel_mask_i;
vessel_mask = any(vessel_mask_layers,3);

if nargout > 1
    vessel_centre_mask = false(rows, cols);
    vessel_centre = (inner_edge+outer_edge)/2;
    [vessel_centre_hi, dists, dists_i] = spline_contour(vessel_centre, [], 1);
    
    x = repmat(1:cols, rows, 1);
    y = repmat((1:rows)', 1, cols);
    num_pts_hi = size(vessel_centre_hi,1);
    for i_pt = 1:num_pts_hi
        pt_mask = ...
            ((x - vessel_centre_hi(i_pt,1)).^2 + (y - vessel_centre_hi(i_pt,2)).^2) < 1;
        vessel_centre_mask = vessel_centre_mask | pt_mask;
        
    end
end

if nargout > 2
   
    vessel_widths = sqrt(sum((inner_edge - outer_edge).^2,2));
    vessel_widths = interp1(dists, vessel_widths, dists_i);
    vessel_norms = compute_spline_normals(vessel_centre_hi);
    vessel_oris = complex(vessel_norms(:,2), vessel_norms(:,1)).^2;  
    
    apex_idx_hi = [1; round(num_pts_hi*apex_idx(:)/num_pts); num_pts_hi];

    ori_map_3d = zeros(rows, cols, num_layers);
    ori_map = exp(2i*rand(rows, cols)*pi);
    
    width_map = zeros(rows, cols);
    for layer = 1:num_layers
        layer_pts = apex_idx_hi(layer):apex_idx_hi(layer+1);
        
        [yi xi] = find(vessel_mask_layers(:,:,layer));
        width_map(vessel_mask_layers(:,:,layer)) = griddata(...
            vessel_centre_hi(layer_pts,1),...
            vessel_centre_hi(layer_pts,2),...
            vessel_widths(layer_pts), xi, yi, 'nearest'); %#ok

        ori_map_i = zeros(rows, cols);        
        ori_map_i(vessel_mask_layers(:,:,layer)) = griddata(...
            vessel_centre_hi(layer_pts,1),...
            vessel_centre_hi(layer_pts,2),...
            vessel_oris(layer_pts), xi, yi, 'nearest'); %#ok
        
        ori_map(vessel_mask_layers(:,:,layer)) =...
            ori_map_i(vessel_mask_layers(:,:,layer));
        ori_map_3d(:,:,layer) = ori_map_i;
    end
       
end


    
    

    
    
