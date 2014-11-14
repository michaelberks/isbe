function [image_out, label, label_centre] = create_ellipse_spicule_mb(halfwidth, contrast, orientation, row, col)
%
%CREATE_ELLIPSE_SPICULE create an image containing a spicule
%
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

[x, y] = meshgrid(1:col, row:-1:1);

orientation = pi*orientation/180;
a = sin(orientation);
b = -cos(orientation);
c = -((a*col/2) + (b*row/2));
dx = abs(a*x + b*y + c);

t = tan(orientation);
if t < (-row/col)
    %intercept bottom
    sy = row;
    sx = -(b*sy + c)/a;
    ey = 1;
    ex = -(b*ey + c)/a;
elseif t < 0
    %intercept left
    sx = 1;
    sy = -(a*sx + c)/b;
    ex = col;
    ey = -(a*ex + c)/b;
elseif t < (row/col)
    %intercept right
    sx = col;
    sy = -(a*sx + c)/b;
    ex = 1;
    ey = -(a*ex + c)/b;
else
    %intercept top
    sy = 1;
    sx = -(b*sy + c)/a;
    ey = row;
    ex = -(b*ey + c)/a;
end
an = -b;
bn = a;
cn =  -(an*sx + bn*sy);
dxn = abs(an*x + bn*y + cn);
dxn = 1 - (dxn / abs(an*ex + bn*ey + cn));

halfwidth_mat = halfwidth*halfwidth*dxn.*dxn;
e_mat = halfwidth_mat - dx.*dx;
e_mat(e_mat<=0) = 0;
image_out = (contrast/halfwidth)*sqrt(e_mat);

%get the ground truth label
label = (image_out > 0);
label_centre = (dx < .5) & label;

return