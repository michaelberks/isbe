function [vessel_mask vessel_centre_mask width_map ori_map ori_map_3d] =...
    make_mask_from_contour(outer_edge, inner_edge, rows, cols)
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

if nargin == 4
    patch_size = [cols rows];
else
    patch_size = ceil(max(outer_edge))+50;    
end

vessel_mask = false(patch_size(2), patch_size(1));
      
num_pts = size(outer_edge,1);   
halfway = floor(num_pts/2);
for i_pt = 1:num_pts-1
    segment_xy = [...
        outer_edge(i_pt,:); inner_edge(i_pt,:);...
        inner_edge(i_pt+1,:); outer_edge(i_pt+1,:)];
    segment_mask = poly2mask(segment_xy(:,1), segment_xy(:,2),...
        patch_size(2), patch_size(1));
    vessel_mask = vessel_mask | segment_mask;
    
    if nargout > 2 
        if i_pt == halfway
            vessel_mask1 = vessel_mask;
            vessel_mask2 = false(patch_size(2), patch_size(1));
        elseif i_pt > halfway
            vessel_mask2 = vessel_mask2 | segment_mask;
        end
    end        
end

if nargout > 1
    vessel_centre_mask = false(patch_size(2), patch_size(1));
    vessel_centre = (inner_edge+outer_edge)/2;
    [vessel_centre, dists, dists_i] = spline_contour(vessel_centre, [], 1);
    
    x = repmat(1:patch_size(1), patch_size(2), 1);
    y = repmat((1:patch_size(2))', 1, patch_size(1));
    num_pts = size(vessel_centre,1);
    for i_pt = 1:num_pts
        pt_mask = ...
            ((x - vessel_centre(i_pt,1)).^2 + (y - vessel_centre(i_pt,2)).^2) < 1;
        vessel_centre_mask = vessel_centre_mask | pt_mask;
        
    end
end

halfway = floor(size(vessel_centre,1)/2);
if nargout > 2
    ori_map1 = exp(2i*rand(patch_size(2), patch_size(1))*pi);
    ori_map2 = exp(2i*rand(patch_size(2), patch_size(1))*pi);
    width_map = zeros(patch_size(2), patch_size(1));
    
    
    vessel_widths = sqrt(sum((inner_edge - outer_edge).^2,2));
    vessel_widths = interp1(dists, vessel_widths, dists_i);
    vessel_norms = compute_spline_normals(vessel_centre);
    vessel_oris = complex(vessel_norms(:,2), vessel_norms(:,1)).^2;
    
%     vessel_centre_idx = sub2ind([patch_size(2) patch_size(1)], ...
%         round(vessel_centre(:,2)), round(vessel_centre(:,1)));
%     seg_idx1 = vessel_mask1(vessel_centre_idx);    
    [yi1 xi1] = find(vessel_mask1);
    
%     width_map(vessel_mask1) = griddata(vessel_centre(seg_idx1,1), vessel_centre(seg_idx1,2), vessel_widths(seg_idx1),...
%         xi1, yi1, 'nearest');       
%     ori_map(vessel_mask1) = griddata(vessel_centre(seg_idx1,1), vessel_centre(seg_idx1,2), vessel_oris(seg_idx1),...
%         xi1, yi1, 'nearest');   
    
    width_map(vessel_mask1) = griddata(vessel_centre(1:halfway,1), vessel_centre(1:halfway,2), vessel_widths(1:halfway),...
        xi1, yi1, 'nearest');       
    ori_map1(vessel_mask1) = griddata(vessel_centre(1:halfway,1), vessel_centre(1:halfway,2), vessel_oris(1:halfway),...
        xi1, yi1, 'nearest');   
    
%     seg_idx2 = ~seg_idx1;
%     vessel_mask2 = vessel_mask & ~vessel_mask1;
    [yi2 xi2] = find(vessel_mask2);
    
%     width_map(vessel_mask2) = griddata(vessel_centre(seg_idx2,1), vessel_centre(seg_idx2,2), vessel_widths(seg_idx2),...
%         xi2, yi2, 'nearest');   
%     ori_map(vessel_mask2) = griddata(vessel_centre(seg_idx2,1), vessel_centre(seg_idx2,2), vessel_oris(seg_idx2),...
%         xi2, yi2, 'nearest');

    width_map(vessel_mask2) = griddata(vessel_centre(halfway+1:end,1), vessel_centre(halfway+1:end,2), vessel_widths(halfway+1:end),...
        xi2, yi2, 'nearest');   
    ori_map2(vessel_mask2) = griddata(vessel_centre(halfway+1:end,1), vessel_centre(halfway+1:end,2), vessel_oris(halfway+1:end),...
        xi2, yi2, 'nearest');
    
    vessel_mask12 = vessel_mask1 & ~vessel_mask2;
    vessel_mask21 = ~vessel_mask1 & vessel_mask2;
    ori_map1(vessel_mask21) = ori_map2(vessel_mask21);
    ori_map2(vessel_mask12) = ori_map1(vessel_mask12);
    
    ori_map = ori_map1;
    ori_map_3d = cat(3, ori_map1, ori_map2);
    
    
    
end
    
    

    
    
