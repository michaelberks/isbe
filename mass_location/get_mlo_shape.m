function [ml_shapes, ml_areas, transforms, count_upper, count_lower] = get_mlo_shape(mlo_list, seg_dir, n_pts, debug_mode)
%GET_MLO_SHAPE *Insert a one line summary here*
%   [ml_shapes,ml_areas] = get_ml_shape(mlo_list,n_pts)
%
% Inputs:
%      mlo_list- *Insert description of input variable here*
%
%      n_pts- *Insert description of input variable here*
%
%
% Outputs:
%      ml_shapes- *Insert description of input variable here*
%
%      ml_areas- *Insert description of input variable here*
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

%Get number of input ml breast borders
nml = length(mlo_list);

%Pre-allocate output structures
ml_shapes = zeros(nml, 6*n_pts);
ml_areas = zeros(nml, 1);

count_upper = 0;
count_lower = 0;
%Go through input list, loading segmented breast border and converting to a
%standard shape
for ii = 1:nml
    
    %Load segmentation
    segmentation = u_load([seg_dir mlo_list(ii).name]);
    
    %Extract fields from saved segmentation
    breast_air = segmentation.breast_air;
    breast_border = segmentation.breast_border;
    seg_size = segmentation.size;
    
    %rows = segmentation.size(1);
    cols = segmentation.size(2);
    pectoral_edge = segmentation.pectoral_edge;
    
    %Work out if right or left breast, and flip border if necessary
    if ~isempty(strfind(mlo_list(ii).name, 'R'));
        breast_border(:,1) = cols - breast_border(:,1) + 1;
        pectoral_edge(:,1) = cols - pectoral_edge(:,1) + 1;
    end

    
    %throw away points at start and finish of breast air that lie on a horizontal edge
%     go_on = true;
%     while go_on
%         if breast_border(breast_air(1),2) < breast_border(breast_air(2),2)
%             go_on = false;
%         else
%             breast_air(1) = [];                 
%         end
%     end
%     go_on = true;
%     while go_on
%         if breast_border(breast_air(end),2) > breast_border(breast_air(end-1),2)
%             go_on = false;
%         else
%             breast_air(end) = [];                 
%         end
%     end
%     

    %Get mask of breast to compute centroid
    bw1 = poly2mask(breast_border(:,1), breast_border(:,2), seg_size(1), seg_size(2));
    bw2 = ~poly2mask([1; pectoral_edge(:,1)], [1; pectoral_edge(:,2)], seg_size(1), seg_size(2));
    
    [yy xx] = find(bw1 & bw2); clear bw*;
    breast_centre = [mean(xx) mean(yy)];

    %Compute theta from pectoral_edge
    theta = atan2(diff(pectoral_edge(:,2)), diff(pectoral_edge(:,1)));
    
    %Compute the rotation needed so the pectoral edge lies vertically in
    %the image
    rot_mat = [cos(pi/2-theta) sin(pi/2-theta); -sin(pi/2-theta) cos(pi/2-theta)];
    
    transforms(ii).rot = rot_mat; %#ok
    
    %Rotate the pectoral edge and compute the translation required to move
    %this to the y-axis
    pectoral_edge_rot = pectoral_edge * rot_mat;
    t = pectoral_edge_rot(1,:);
    transforms(ii).t = t; %#ok
    
    pectoral_edge_rot = pectoral_edge_rot - [t;t];
    
    %Rotate and translate the breast border
    breast_border_rot = breast_border * rot_mat - repmat(t, size(breast_border,1),1);
    
    xx = breast_border_rot(breast_air,1);
    yy = breast_border_rot(breast_air,2);

    xx = [xx(2); xx(1); xx; xx(end); xx(end-1);]; %#ok
    xt = xx(2:end-3) - xx(4:end-1);
    xtt = xx(1:end-4) - xx(5:end);
    yy = [yy(2); yy(1); yy; yy(end); yy(end-1);]; %#ok
    yt = yy(2:end-3) - yy(4:end-1);
    ytt = yy(1:end-4) - yy(5:end);

    xx([1:2, end-1:end],:) = [];
    yy([1:2, end-1:end],:) = [];
 
    k = (xt.* ytt + yt.*xtt) ./ ((xt.^2 + yt.^2).^(3/2));

    h = fspecial('gaussian', [30 1], 5);
    k_smooth = conv2(k,h, 'same');

    xyn = [-yt xt] ./ [sqrt(xt.^2 + yt.^2) sqrt(xt.^2 + yt.^2)];
    xd = [0; diff(xx)];
    yd = [0; diff(yy)];
    D = sum([xd yd] .* xyn, 2);
    D_smooth = conv2(D,h, 'same');
 
    halfway = round(length(D)/2);
    [min_up min_up_idx] = min(D_smooth(1:halfway));
    [min_do min_do_idx] = min(D_smooth(halfway+1:end));
    min_do_idx = min_do_idx + halfway;
    
    
    discard = false(length(breast_air),1);
    if min_up < 0
        discard(1:min_up_idx-1) = 1;
        count_upper = count_upper + 1;
    end
    if min_do < 0
        discard(min_do_idx+1:end) = 1;
        count_lower = count_lower + 1;
    end
    breast_air(discard) = [];
    
    if debug_mode
        figure; 
        subplot(2,2,[1 3]); axis equal ij; hold on;
        plot(breast_border_rot(:,1), breast_border_rot(:,2));
        plot([xx + xyn(:,1) xx - xyn(:,1)], [yy+xyn(:,2) yy-xyn(:,2)], 'g:');
        plot(pectoral_edge_rot(:,1), pectoral_edge_rot(:,2), 'r');
        if min_up < 0
            plot(xx(min_up_idx), yy(min_up_idx), 'ro');
        end
        if min_do < 0
            plot(xx(min_do_idx), yy(min_do_idx), 'ro');
        end

        subplot(2,2,2); plot(k); hold on;
        plot(k_smooth, 'r:');

        subplot(2,2,4); plot(D); hold on;
        plot(D_smooth, 'r:');
    end
    
    %Extend the start and end of the breast air points to the pectoral edge 
    %p_start = polyfit(breast_border_rot(breast_air(1:5),1),breast_border_rot(breast_air(1:5),2),1);
    %breast_air_start = [1 p_start(2)];
    
    p_end = polyfit(breast_border_rot(breast_air(end-4:end),1),breast_border_rot(breast_air(end-4:end),2),1);
    breast_air_end = [1 p_end(2)];
    
    %breast_air_pts = [breast_air_start; breast_border_rot(breast_air,:); breast_air_end];
    breast_air_pts = [breast_border_rot(breast_air,:); breast_air_end];
    %breast_air_pts = breast_border_rot(breast_air,:);
    
    %Find the point on breast air closest to the intersection of a line
    %running normal from the chest wall through the centroid
    %The line by construction is horizontal, so just need to find the
    %breast air point with closest matching y to centre_y
    [dummy idx] = min((breast_air_pts(:,2)-breast_centre(2)).^2);
    
    %Split breast border into 3 segemnts - upper breast air, lower breast
    %air and chest wall and space 50 points equally along each
    segment1 = breast_air_pts(1:idx,:);
    segment2 = breast_air_pts(idx+1:end,:);
    
    temp_y = linspace(breast_air_pts(end,2), breast_air_pts(1,2), n_pts+2)';
    segment3 = [zeros(n_pts,1) temp_y(2:end-1)];
    
    diff1 = diff(segment1);
    cum_dist1 = [0; cumsum(sqrt(sum(diff1.^2,2)))];
    segment1 = interp1(cum_dist1, segment1, linspace(0, cum_dist1(end), n_pts), 'linear');
    
    diff2 = diff(segment2);
    cum_dist2 = [0; cumsum(sqrt(sum(diff2.^2,2)))];
    segment2 = interp1(cum_dist2, segment2, linspace(0, cum_dist2(end), n_pts), 'linear');
    
    %Arrange segments in output structure
    ml_shapes(ii,:) = ...
        [segment1(:,1)' segment2(:,1)' segment3(:,1)' ...
         segment1(:,2)' segment2(:,2)' segment3(:,2)'];
     
    %Compute the area of the new shape
    ml_areas(ii) = polyarea(ml_shapes(ii,1:end/2), ml_shapes(ii,end/2+1:end));
    
    %If we're debugging, plot border and the segments we've formed
    if debug_mode
        figure; axis ij equal; hold on;
        %plot([1 c c 1 1], [1 1 r r 1], 'b');
        plot(breast_border_rot(:,1), breast_border_rot(:,2), 'k.');
        plot(segment1(:,1), segment1(:,2), 'r');
        plot(segment2(:,1), segment2(:,2), 'g');
        plot(segment3(:,1), segment3(:,2), 'y');
        plot(segment1(:,1), segment1(:,2), 'rx', 'linewidth', 2);
        plot(segment2(:,1), segment2(:,2), 'gx', 'linewidth', 2);
        plot(segment3(:,1), segment3(:,2), 'yx', 'linewidth', 2);
        plot(breast_centre(1), breast_centre(2), 'r*', 'MarkerSize', 10);
        plot([1 breast_air_pts(idx,1)], [breast_air_pts(idx,2) breast_air_pts(idx,2)], 'c');
    end
    
end
