function [offset, theta] = register_frame(frame1, frame2, offset_lim, theta_lim, plot_flag)
%REGISTER_FRAME *Insert a one line summary here*
%   [offset, theta] = register_frame(frame1, frame2, offset_lim, theta_lim)
%
% Inputs:
%      frame1 - *Insert description of input variable here*
%
%      frame2 - *Insert description of input variable here*
%
%      offset_lim - *Insert description of input variable here*
%
%      theta_lim - *Insert description of input variable here*
%
%
% Outputs:
%      offset - *Insert description of input variable here*
%
%      theta - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 23-Aug-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

%--------------------------------------------------
%Compute line points in first frame

%Apply Gaussian 2nd derivatives
[line_strength, line_orientation] = ...
    gaussian_clover_line(frame1, [4 8]);%gaussian_2nd_derivative_line2

%Apply non-maximal suppression to skeletonise the line strength map
line_nms = mb_non_maximal_supp(line_strength, line_orientation);

%Discard points from edges of map
line_nms([1:32 end-32:end], :) = 0;
line_nms(:, [1:32 end-32:end]) = 0;

%Apply hysterisis to select suitable lines from the NMS map
[line_mask] = hysterisis(line_nms);

%Extract (x,y) coordinates of the remaining lines
[y1 x1] = find(line_mask);

%--------------------------------------------------
%Repeat the process for the 2nd frame
%Apply Gaussian 2nd derivatives
[line_strength, line_orientation] = ...
    gaussian_clover_line(frame2, [4 8]);%gaussian_2nd_derivative_line2

%Apply non-maximal suppression to skeletonise the line strength map
line_nms = mb_non_maximal_supp(line_strength, line_orientation);

%Discard points from edges of map
line_nms([1:32 end-32:end], :) = 0;
line_nms(:, [1:32 end-32:end]) = 0;

%Apply hysterisis to select suitable lines from the NMS map
[line_mask] = hysterisis(line_nms);

%Extract (x,y) coordinates of the remaining lines
[y2 x2] = find(line_mask);

%--------------------------------------------------
%Generate offset and theta ranges
if nargin < 3 || isempty(offset_lim)
    offset_lim = 10; %Use default 10 if range not user specified
end
offset_sz = 2*offset_lim + 1;

if nargin < 4 || isempty(theta_lim)
    theta_lim = 10; %Use default 10 if range not user specified
end
theta_range = pi*(-theta_lim:theta_lim)/180;
theta_sz = length(theta_range);

%Compute the centre of the first frame
xc = size(frame1,2)/2;
yc = size(frame1,1)/2;

%--------------------------------------------------
%Loop through thetas finding rotation that highest offset max

max_count = 0;
for tt = 1:theta_sz
    theta = theta_range(tt);
    
    %Generate rotation matrix for this theta
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    
    %Rotate line points extracted from frame 1 - need to centre by (xc, yc)
    %first
    xyt = [x1-xc y1-yc]*R;

    %Set up hough matrix to count the votes for each offset
    offset_counts = zeros(offset_sz);
    
    %Loop through each line point from frame 1 and compute the offset to all
    %points in frame 2
    for ii = 1:length(x1);
        
        %Computes offset and decentre
        x_offset = round(x2 - xyt(ii,1) - xc);
        y_offset = round(y2 - xyt(ii,2) - yc);
        
        %Only count votes less than the offset limit
        keep = abs(x_offset) <= offset_lim & abs(y_offset) <= offset_lim;
        
        %Use sparse matrix trick to counts votes
        if any(keep)
            offset_counts = offset_counts + sparse(...
                y_offset(keep)+offset_lim+1,...
                x_offset(keep)+offset_lim+1,...
                1, offset_sz, offset_sz);
        end
    end

    %Workout the offset with maximum vote for this theta
    [max_count_theta max_idx_theta] = max(offset_counts(:));

    %If this vote count is larger than the existing maximum record the
    %offset index and theta
    if max_count_theta > max_count
        max_count = max_count_theta;
        max_idx = max_idx_theta;
        max_theta = theta;
        
    end
end

%Convert the maximum offset ID
[max_row max_col] = ind2sub([offset_sz offset_sz], max_idx);
offset = [max_col max_row] - offset_lim - 1;

%Return theta in degrees
theta = 180*max_theta/pi;

if plot_flag
    R = [cos(max_theta) -sin(max_theta); sin(max_theta) cos(max_theta)];
    xyt = [x1-xc y1-yc]*R;
    x_trans = xyt(:,1) + xc + offset(:,1);
    y_trans = xyt(:,2) + yc + offset(:,2);
            
    figure; imagesc(frame1); hold on; axis image; colormap(gray(256));
    plot(x1, y1, 'r.', 'markersize', 2);
    figure; imagesc(frame2); hold on; axis image; colormap(gray(256));
    plot(x2, y2, 'g.', 'markersize', 2);
    plot(x1, y1, 'rx', 'markersize', 2);
    plot(x_trans, y_trans, 'b.', 'markersize', 2);
            
end

