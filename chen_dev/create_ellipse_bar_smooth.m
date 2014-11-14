function [image_out, label, label_centre] = create_ellipse_bar_smooth(halfwidth, contrast, orientation, row, col)
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

orientation = pi*orientation/180;
a = sin(orientation);
b = cos(orientation);
c = -((a*col/2) + (b*row/2));
dx = abs(a*x + b*y + c);
ex = halfwidth*halfwidth-dx.*dx;
ex(ex<=0) = 0;

image_out = (contrast/halfwidth)*sqrt(ex);
h = fspecial('average', 7);
image_out = imfilter(image_out, h);

% % get the ground truth label
label = (image_out > 0);

label_centre = (dx > -.5) & (dx < .5);


return