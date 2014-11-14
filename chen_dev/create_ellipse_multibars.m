function [image_out, label, label_centre] = create_ellipse_multibars(halfwidth, contrast, orientation, row, col, centre_x, centre_y)
% 
%CREATE_ELLIPSE_MULTIBARS create an image containing multiple ellipse bars (the intensity
%has ellipse distribution)
% USAGE:
%   [image_out, label, label_centre] = create_ellipse_multibars(halfwidth,
%   contrast, orientation,  row, col, centre_x, centre_y)
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
%      centre_x, centre_y - a point which the line through it.
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
% Notes: see also sample_testing_image_multibars.m
%  
% Created: 29-March-2010
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
c = -((a*centre_x) + (b*centre_y));
dx = abs(a*x + b*y + c);
ex = halfwidth*halfwidth-dx.*dx;
ex(ex<=0) = 0;

image_out = (contrast/halfwidth)*sqrt(ex);

% % get the ground truth label
label = (image_out > 0);

label_centre = (dx > -.5) & (dx < .5);


return