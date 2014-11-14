function [image_out, label, label_centre, label_orientation] =...
    create_sin_curve(halfwidth, contrast, radius, tangent, squash, row, col, centre_x, centre_y)
%
%CREATE_SIN_CURVE create an image containing a curved linear structure with
% a cross-sectional profile that varies between sinusoidal and rectangular
% - depending on the value of squash
% USAGE:
%   [image_out, label, label_centre] = create_sin_curve(halfwidth,
%   contrast, radius, tangent, squash, row, col, centre_x, centre_y)
%
% Inputs:
%      halfwidth - halfwidth of cross sectional profile
%
%      contrast - maximum height of cross sectional profile
%
%      radius - radius of the curvature
%
%      tangent - orientation of tangent vector at bar centre
%
%      squash - controls whether the line has a sinusoidal profile (squash
%      = 0) or a rectangular profile (=1). Values inbetween [0, 1] generate
%      a linear interpolation between the two extremes
%
%      row/col - number of rows/cols in image
%
%      centre_x/y - centre of the line in the image
%
%
% Outputs:
%      image_out - image containing line
%
%      label - label of line (1) vs background (0) 
%
%      label_centre - the label of centre line (1) vs background (0)
%
%      label_orientation - orientation (in degrees 0-180) at ecah pixel
%      (only meaningful for line points)
%
% Example:
%
% Notes:
%
% Created: 04-February-2010
% Author: Michael Berks
% Email : michael.berks@manchester.ac.uk
% Phone : +44 (0)161 275 7669
% Copyright: (C) University of Manchester
% 
if nargin < 7
    centre_x = col/2;
    centre_y = row/2;
end

tangent = pi*tangent/180;

centre_x = radius*sin(tangent) + centre_x;
centre_y = radius*cos(tangent) + centre_y;

% create co-ordinates
[x, y] = meshgrid(1:col, 1:row);

dx = sqrt((x-centre_x).^2 + (y-centre_y).^2) - radius;

label = (dx > -halfwidth) & (dx < halfwidth);

ex = 0.5 * (pi - pi*dx/halfwidth);
base_wave = sin(ex);

% for ii = 3:2:99; 
%     base_wave = base_wave + squash*sin(ii*ex)/ii; 
% end
% image_out = contrast*base_wave / max(base_wave(:));
    
image_out = contrast*(squash*label + (1-squash)*base_wave);
image_out(~label) = 0;

if nargout > 2
    %get the ground truth label for the centre line
    label_centre = (-.5 < dx) & (dx < .5);
end
if nargout > 3
    %get the ground truth label for the centre line
    label_orientation = 180*atan((centre_x-x)./(centre_y-y)) / pi;
    idx = label_orientation < 0;
    label_orientation(idx) = label_orientation(idx) + 180;
end

return