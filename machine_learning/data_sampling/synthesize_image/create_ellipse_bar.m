function [image_out, label, label_centre, label_orientation] = ...
    create_ellipse_bar(halfwidth, contrast, orientation, ...
                       row, col, centre_x, centre_y)
%
%CREATE_ELLIPSE_BAR create an image containing a ellipse bar (the intensity
%has ellipse distribution)
% USAGE:
%   [image_out, label, label_centre] = create_ellipse_spicule(halfwidth,
%   contrast, orientation, row, col)
%
% Inputs:
%      halfwidth - halfwidth of Gaussian profile at half its maximum height
%
%      contrast - maximum height of (scaled) Gaussian profile
%
%      orientation - orientation of bar in image in degrees
%
%      row - number of rows in image
%
%      col - number of columns in image
%
%
% Outputs:
%      image_out - image containing spicule
%
%      label - label of spicule (1) vs background (0) 
%
%      label_centre - the label of centre line (1) vs background (0)
%
% Example:
%
% Notes:  
%  
% Created: 04-February-2010
% Author: Zezhi Chen
% Email : zezhi.chen@manchester.ac.uk
% Phone : +44 (0)161 275 7669
% Copyright: (C) University of Manchester
% 

% create co-ordinates

[x, y] = meshgrid(1:col, 1:row);
a = sin(pi*orientation/180);
b = cos(pi*orientation/180);
c = -((a*centre_x) + (b*centre_y));
dx = a*x + b*y + c;
ex = halfwidth*halfwidth - dx.*dx;
ex(ex<=0) = 0;

image_out = (contrast/halfwidth)*sqrt(ex);

%get the ground truth label
if nargout > 1
    %get the ground truth label for the whole line
    label = (image_out > 0);
end
if nargout > 2
    %get the ground truth label for the centre line
    label_centre = (-.5 < dx) & (dx < .5);
end
if nargout > 3
    %get the ground truth label for the centre line
    label_orientation = nan(row, col);
    label_orientation(label) = mod(orientation,180);
end

return