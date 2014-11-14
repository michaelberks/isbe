function [apex_offset_map] = ...
    transform_apex_offset_preds(apex_class_pred, apex_offset_x_pred, apex_offset_y_pred,...
        vessel_centre, nrows, ncols, base_width, include_pts, separate_trees)
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

num_pts = length(vessel_centre.x);

if ~exist('separate_trees', 'var') || isempty(separate_trees)
    separate_trees = 1;
end
if ~separate_trees && size(apex_offset_x_pred, 2) > 1
    apex_offset_x_pred = mean(apex_offset_x_pred, 2);
    apex_offset_y_pred = mean(apex_offset_y_pred, 2);
end
num_preds = size(apex_offset_x_pred, 2);

%Now we can loop through all these points again and update the apex
%location map using the offset predictions and weightings
apex_offset_map = zeros(nrows, ncols);
for i_pt = 1:num_pts
    if include_pts(i_pt)
        vxc = vessel_centre.x(i_pt);
        vyc = vessel_centre.y(i_pt);
        ori_c = angle(vessel_centre.ori(i_pt))/2;
        scale = vessel_centre.width(i_pt) / base_width;

        %Get scale relative to base width a make rotation matrix
        rot = [cos(ori_c) -sin(ori_c); sin(ori_c) cos(ori_c)];

        for i_pred = 1:num_preds
            %Translate predicted offset from the frame coordinates into real
            %image coordinates
            apex_xy_t = [apex_offset_x_pred(i_pt,i_pred) apex_offset_y_pred(i_pt,i_pred)];
            apex_xy_p = scale*apex_xy_t*rot + [vxc vyc];                

            %Add a blurred vote to the image (width of blur depends on scale of
            %frame)
            [im_r, im_c, gxy] = add_gaussian_blur(apex_xy_p(1), apex_xy_p(2), 2*scale, nrows, ncols);
            if ~isempty(gxy)
                apex_offset_map(im_r, im_c) = apex_offset_map(im_r, im_c) + apex_class_pred(i_pt,1)*gxy;
            end
        end
    end
end          
    
    
