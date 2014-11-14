function [] = ...
    incremental_apex_offset_preds(apex_class_pred, apex_offset_x_pred, apex_offset_y_pred,...
        vessel_centre, nrows, ncols, include_pts, class_thresh, candidates_dir, im_name, varargin)
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
    0, ... % non-strict mode
    'separate_trees',       1,...
    'base_width',           20,...
    'exclsuion_zone',       20,...
    'fov_mask',             [],...
    'plot',              []);

num_pts = length(vessel_centre.x);

if ~args.separate_trees && size(apex_offset_x_pred, 2) > 1
    apex_offset_x_pred = mean(apex_offset_x_pred, 2);
    apex_offset_y_pred = mean(apex_offset_y_pred, 2);
end
num_preds = size(apex_offset_x_pred, 2);

num_thresh = length(class_thresh);
class_thresh = sort(class_thresh, 'descend');

if args.plot
    figure; axis equal ij; hold all;
end

%Now we can loop through all these points again and update the apex
%location map using the offset predictions and weightings
apex_offset_map = zeros(nrows, ncols);
for i_th = 1:num_thresh
    
    if ~exist([candidates_dir '\thresh' zerostr(i_th,2)], 'dir')
        mkdir([candidates_dir '\thresh' zerostr(i_th,2)]);
    end
    
    if i_th == 1
        include_pts_i = include_pts & (apex_class_pred >= class_thresh(i_th));
    else
        include_pts_i = include_pts & (apex_class_pred >= class_thresh(i_th)) & (apex_class_pred < class_thresh(i_th-1));
    end    
    
    if args.plot && any(include_pts_i)
        plot(vessel_centre.x(include_pts_i), vessel_centre.y(include_pts_i), '.');
    end
    
    for i_pt = 1:num_pts
        if include_pts_i(i_pt)
            vxc = vessel_centre.x(i_pt);
            vyc = vessel_centre.y(i_pt);
            ori_c = angle(vessel_centre.ori(i_pt))/2;
            scale = vessel_centre.width(i_pt) / args.base_width;

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
    
    [candidate_xy candidate_scores] = ...
        local_image_maxima(apex_offset_map, args.exclsuion_zone, args.fov_mask, 0); %#ok
    
    save([candidates_dir '\thresh' zerostr(i_th,2) '\' im_name '_candidates'],...
        'candidate_xy', 'candidate_scores');
    
end
if ~exist([candidates_dir '\apex_class_thresh.mat'], 'file')
    save([candidates_dir '\apex_class_thresh.mat'], 'class_thresh');
end
    
    
