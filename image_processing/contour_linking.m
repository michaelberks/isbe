function [groups_map group_vecs] = contour_linking(ori_map, ridge_map)
%CONTOUR_LINKING *Insert a one line summary here*
%   [group_map] = contour_linking(ori_map, ridge_map)
%
% Inputs:
%      ori_map - *Insert description of input variable here*
%
%      ridge_map - *Insert description of input variable here*
%
%
% Outputs:
%      groups_map - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 18-Nov-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

eta_c = 5;
eta_p = 0;
eta_o = 0.95;

%get size of the image
[rows cols] = size(ridge_map);

%pre-compue cos and sins of the orientations
g2_c = cos(ori_map);
g2_s = sin(ori_map);

%find all the initial ridge pixels (as x,y and indices)
[ridge_pixels_y ridge_pixels_x] = find(ridge_map);
[ridge_pixels_idx] = sub2ind([rows cols], ridge_pixels_y, ridge_pixels_x);

%Set up the listed of ugrouped pixels, at first, all are ungrouped
ungrouped = 1:length(ridge_pixels_idx);

%Set up map of groups
groups_map = zeros(rows, cols);
group_vecs = cell(0,1);
group_idx = 0;
select_new = 1;
while ~isempty(ungrouped)
    
    if select_new
        %Select a random point from the unrouped pixels
        rand_idx = randperm(length(ungrouped));
        next_g_idx = ungrouped(rand_idx(1));
        group_idx = group_idx + 1;
        group_vecs{group_idx} = [];
    end
    
    %Remove it from the grouped pixels
    ungrouped(ungrouped == next_g_idx) = [];
    x_g = ridge_pixels_x(next_g_idx(1));
    y_g = ridge_pixels_y(next_g_idx(1));
    
    %set the group map to the current group index
    groups_map(y_g, x_g) = group_idx;
    group_vecs{group_idx}(end+1,:) = [x_g y_g];
    
    %Get vectors and distances to all ungrouped pixels
    r = [x_g-ridge_pixels_x(ungrouped) y_g-ridge_pixels_y(ungrouped)];
    r_dist = sqrt(sum(r.^2,2));
    
    %Set select_new to true - it will be turned back to false if we find a
    %suitable pixel
    select_new = 1;
    
    %Check if we have any ungrouped pixels in local vicinty
    local = r_dist < eta_c;
    if any(local)
        
        %if so, get the orientation vectors of the grouped pixel and the
        %local ungrouped pixels
        u_idx = ungrouped(local);
        xg_c = g2_c(y_g, x_g);
        xg_s = g2_s(y_g, x_g);
        xu_c = g2_c(ridge_pixels_idx(u_idx));
        xu_s = g2_s(ridge_pixels_idx(u_idx));

        %compute the alignment of the orientations
        aligned = abs(xg_c*xu_c + xg_s*xu_s) > eta_o;
        
        %check if any are aligned
        if any(aligned)
            
            %if so, compute the alignment of the grouped orientation and
            %the vector to the ungrouped pixel
            r_dist = r_dist(aligned);
            r = r(aligned,:);
            u_idx = u_idx(aligned);
            r_aligned = (abs(xg_c*r(:,1) + xg_s*r(:,2)) ./ r_dist) > eta_p;
            
            %check if these aren't aligned
            if any(r_aligned)
                %if so we have candidates, so select the closest as the new
                %grouped pixel
                next_g_idx = u_idx(r_aligned);
                r_dist = r_dist(r_aligned);
                [dummy min_dist_idx] = min(r_dist);
                next_g_idx = next_g_idx(min_dist_idx);
                
                %Set select new to false
                select_new = 0;
            end
        end
    end
end

%         %Add any other pixels connected to the current group
%         [rgroup cgroup] = find(groups_map == group_idx);
%         groups_map(bwselect(ridge_map, cgroup, rgroup, 8)) = group_idx;
%         
%         %find all the initial ridge pixels (as x,y and indices)
%         [ridge_pixels_y ridge_pixels_x] = find(ridge_map);
%         [ridge_pixels_idx] = sub2ind([rows cols], ridge_pixels_y, ridge_pixels_x);
% 
%         %Set up the listed of ugrouped pixels, at first, all are ungrouped
%         ungrouped = 1:length(ridge_pixels_idx);
%         
