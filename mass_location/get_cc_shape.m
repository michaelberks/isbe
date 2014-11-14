function [cc_shapes, cc_areas] = get_cc_shape(cc_list, seg_dir, n_pts, debug_mode)
%GET_CC_SHAPE *Insert a one line summary here*
%   [cc_shapes,cc_areas] = get_cc_shape(cc_list,n_pts)
%
% Inputs:
%      cc_list- *Insert description of input variable here*
%
%      n_pts- *Insert description of input variable here*
%
%
% Outputs:
%      cc_shapes- *Insert description of input variable here*
%
%      cc_areas- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 17-Nov-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%Set default arguments
if nargin < 3
    n_pts = 50;
end
if nargin < 4
    debug_mode = 0;
end

%Get number of input cc breast borders
ncc = length(cc_list);

%Pre-allocate output structures
cc_shapes = zeros(ncc, 6*n_pts);
cc_areas = zeros(ncc, 1);

%Go through input list, loading segmented breast border and converting to a
%standard shape
for ii = 1:ncc
    
    %Load segmentation
    segmentation = u_load([seg_dir cc_list(ii).name]);
    
    %Extract fields from saved segmentation
    breast_air = segmentation.breast_air;
    breast_border = segmentation.breast_border;
    
    r = segmentation.size(1);
    c = segmentation.size(2);
    
    %Work out if right or left breast, and flip border if necessary
    if ~isempty(strfind(cc_list(ii).name, 'R'));
        breast_border(:,1) = c - breast_border(:,1) + 1;
    end
    
    %throw away points at start and finish of breast air that lie on a horizontal edge
    go_on = true;
    while go_on
        if breast_border(breast_air(1),2) < breast_border(breast_air(2),2)
            go_on = false;
        else
            breast_air(1) = [];                 
        end
    end
    go_on = true;
    while go_on
        if breast_border(breast_air(end),2) > breast_border(breast_air(end-1),2)
            go_on = false;
        else
            breast_air(end) = [];                 
        end
    end
    
    %Extend the start and end of the breast air point to the chest wall 
    p_start = polyfit(breast_border(breast_air(1:5),1),breast_border(breast_air(1:5),2),1);
    p_end = polyfit(breast_border(breast_air(end-4:end),1),breast_border(breast_air(end-4:end),2),1);
    breast_air_start = [1 sum(p_start)];
    breast_air_end = [1 sum(p_end)];
    breast_air_pts = [breast_air_start; breast_border(breast_air,:); breast_air_end];
    
    %Get list of points inside polygon defined by breast air to get an
    %estimate of area and centroid
    xe = ceil(max(breast_air_pts(:,1)));
    ys = floor(min(breast_air_pts(:,2)));
    ye = ceil(max(breast_air_pts(:,2)));
    
    [bw] = poly2mask(breast_air_pts(:,1), breast_air_pts(:,2) - ys + 1, ye - ys + 1, xe);
    [yy xx] = find(bw);
    centre_x = mean(xx);
    centre_y = mean(yy);
    
    cc_areas(ii) = polyarea(breast_air_pts(:,1), breast_air_pts(:,2));
        
    clear xx yy bw
    centre_y = centre_y - 1 + ys;

    
    %Find the point on breast air closest to the intersection of a line
    %running normal from the chest wall through the centroid
    %The line by construction is horizontal, so just need to find the
    %breast air point with closest matching y to centre_y
    [dummy idx] = min((breast_air_pts(:,2)-centre_y).^2);
    
    %Split breast border into 3 segemnts - upper breast air, lower breast
    %air and chest wall and space 50 points equally along each
    segment1 = breast_air_pts(1:idx,:);
    segment2 = breast_air_pts(idx+1:end,:);
    
    temp_y = linspace(breast_air_pts(end,2), breast_air_pts(1,2), n_pts+2)';
    segment3 = [ones(n_pts,1) temp_y(2:end-1)];
    
    diff1 = diff(segment1);
    cum_dist1 = [0; cumsum(sqrt(sum(diff1.^2,2)))];
    segment1 = interp1(cum_dist1, segment1, linspace(0, cum_dist1(end), n_pts), 'linear');
    
    diff2 = diff(segment2);
    cum_dist2 = [0; cumsum(sqrt(sum(diff2.^2,2)))];
    segment2 = interp1(cum_dist2, segment2, linspace(0, cum_dist2(end), n_pts), 'linear');
    
    %Arrange segments in output structure
    cc_shapes(ii,:) = ...
        [segment1(:,1)' segment2(:,1)' segment3(:,1)' ...
         segment1(:,2)' segment2(:,2)' segment3(:,2)'];
    
    %If we're debugging, plot border and the segments we've formed
    if debug_mode
        figure; axis ij equal; hold on;
        plot([1 c c 1 1], [1 1 r r 1], 'b');
        plot(breast_border(:,1), breast_border(:,2), 'k.');
        plot(segment1(:,1), segment1(:,2), 'r');
        plot(segment2(:,1), segment2(:,2), 'g');
        plot(segment3(:,1), segment3(:,2), 'y');
        plot(segment1(:,1), segment1(:,2), 'rx', 'linewidth', 2);
        plot(segment2(:,1), segment2(:,2), 'gx', 'linewidth', 2);
        plot(segment3(:,1), segment3(:,2), 'yx', 'linewidth', 2);
        plot(centre_x, centre_y, 'r*', 'MarkerSize', 10);
        plot([1 breast_air_pts(idx,1)], [breast_air_pts(idx,2) breast_air_pts(idx,2)], 'c');
    end
    
end
