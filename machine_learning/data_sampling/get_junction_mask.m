function [junction_mask] = get_junction_mask(line_mask, ...
                                             centre_mask, ...
                                             output_type)
% Given masks of vessels and their centrelines, find the junctions and
% return a binary image indicating all pixels that can be considered as
% being at a junction.

% Get xy coords of points
[c_y c_x] = find(centre_mask);
num_pts = length(c_x);

% Distance transform of vessel mask - gives Euclidean distance of
% every point to background segment
v_thick = bwdist(~line_mask);

% Create a mask of size ws x ws, with ones on the border and zeros
% in the centre
ws = 5;
ws2 = floor(ws/2);
offs = -ws2:ws2;
count_mask = ones(ws); 
count_mask(2:ws-1,2:ws-1) = 0;

% Set up meshgrid and preallocate junction mask image
[rows cols] = size(line_mask);
[xx,yy] = meshgrid(1:cols, 1:rows);
junction_mask = false(rows, cols);

for ii = 1:num_pts
    p_rows = min(max(c_y(ii) + offs, 1),size(centre_mask,1));
    p_cols = min(max(c_x(ii) + offs, 1),size(centre_mask,2));
    
    patch = centre_mask(p_rows, p_cols);
    patch = bwselect(patch, ws2+1, ws2+1, 8) & count_mask;
    label = bwlabel(patch, 8);

    % if more than two vessel sprout from this point then add points
    % to the junction mask
    if max(label(:)) > 2
        if strcmp(output_type, 'junction_detection')
            % add a fat circle that lies just within the vessel
            radius = v_thick(c_y(ii), c_x(ii));
            circ = (xx - c_x(ii)).^2 + (yy - c_y(ii)).^2 < radius^2;
            circ = circ & line_mask;
            junction_mask = junction_mask | circ;
        else
            % just add the centreline point
            junction_mask(c_y(ii) + offs, c_x(ii) + offs) = true;
        end
    end
end
