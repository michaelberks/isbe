function rad_scores = hough_circle(bw, centres, min_rad, max_rad)
%HOUGH_CIRCLE given a binary edge map compute hough transform to determine
%votes for circles in the image
%
% Inputs:
%      bw - a binary edge map of the image
%
%      centres - a 2xN array of possible (x,y) centres for the circles. If
%      not given defaults to all possible centres in the image;
%
%      min_rad - the minimum radius of circles that can be voted for -
%      default 5 pixels
%
%      max_rad - the maximum radius of circles - defaults to half the
%      smaller of the two dimensions of the image
%
% Outputs:
%
%      rad_scores - an array of votes where each row corresponds to a
%      possible centre and each column to a radius given that centre
%
%
%
% Example:
%
% Notes: This algorithm was written to find a circle within in a fairly
% constrained region. As a result it loops through all possible centres for
% each edge point, with votes cast for possible radii. However if you
% usually need to search over the whole image you may want to think about 
% restructuring the algorithm as this code may be very slow
%
% See also: HOUGH
%
% Created: 11-Oct-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%Get image size
[rows cols] = size(bw);

%Set defaults for missing arguments - can't be used to a full
%name/parameter search
if nargin < 2 
    %set centres to all possible (x,y) locations in image
    x = repmat(1:cols, rows, 1);
    y = repmat((1:rows)', 1, cols);
    centres = [x(:) y(:)];
    clear x y;
    
    %note we have to loop through the centres so this may make the
    %algorithm very slow
end
if nargin < 3
    %set default minimum radius to 5
    min_rad = 5;
end
if nargin < 4
    %set default maximum radius to half the smaller of the two image
    %dimensions
    max_rad = round(min([rows cols])/2);
end
    
%Find all edge pts in BW
[ey ex] = find(bw);

%Create table to store votes in
rad_scores = zeros(size(centres,1), max_rad - min_rad + 1); 

tolerance = 2;

%For each edge pt
for ii = 1:length(ex)

    %Loop through each radii selecting centres that match
    for r = min_rad:max_rad
        idx = abs( (centres(:,1) - ex(ii)).^2 + (centres(:,2) - ey(ii)).^2 - r.^2 ) < tolerance.^2;
        rad_scores(idx,r-min_rad+1) = rad_scores(idx,r-min_rad+1) + 1; 
    end

%     %Work out distance of point to all possible centres
%     dists = sqrt((ex(ii) - centres(:,1)).^2 + (ey(ii) - centres(:,2)).^2);
%     %For each centre
%     for jj = 1:size(centres, 1)
%         
%         %If the edge point is a valid distance from the centre...
%         if dists(jj) > min_rad && dists(jj) < max_rad
%             
%             %Vote for this radius/centre combination (note because circles
%             %lie on a pixel grid as opposed to a continuous space we have
%             %to round the computed radius - in this case we take both the
%             %floor and ceiling and assign a vote to both)
%             ru = ceil(dists(jj)) - min_rad + 1;
%             rd = floor(dists(jj)) - min_rad + 1;
%             rad_scores(jj, ru) = rad_scores(jj, ru)+1;
%         	rad_scores(jj, rd) = rad_scores(jj, rd)+1;
%         end
%     end
    
end

        
        