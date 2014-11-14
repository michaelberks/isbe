function [image_out, label, label_centre, label_orientation] =...
    create_sin_bar(halfwidth, contrast, orientation, row, col, squash, centre_x, centre_y)
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
if nargin < 7
    centre_x = col/2;
    centre_y = row/2;
end
% create co-ordinates
[x, y] = meshgrid(1:col, 1:row);

a = sin(pi*orientation/180);
b = cos(pi*orientation/180);
c = -((a*centre_x) + (b*centre_y));
dx = a*x + b*y + c;

label = (dx > -halfwidth) & (dx < halfwidth);

ex(label) = 0.5 * (pi - pi*dx(label)/halfwidth);
base_wave(label) = sin(ex(label));

% for ii = 3:2:99; 
%     base_wave = base_wave + squash*sin(ii*ex)/ii; 
% end
% image_out = contrast*base_wave / max(base_wave(:));

image_out = zeros(size(label));
image_out(label) = contrast*(squash + (1-squash)*base_wave(label));

%get the ground truth label
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