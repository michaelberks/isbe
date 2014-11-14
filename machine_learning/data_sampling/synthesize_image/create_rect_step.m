function [image_out, label] = create_rect_step(halfwidth, contrast, orientation, row, col, centre_x, centre_y)
%CREATE_RECT_BAR *Insert a one line summary here*
%   [image_out,label] = create_rect_bar(halfwidth,contrast,orientation,row,col)
%
% Inputs:
%      halfwidth- halfwidth of the bar
%
%      contrast- maximum height of the bar
%
%      orientation- orientation of bar in image
%
%      row- number of rows in image
%
%      col- number of columns in image
%
%
% Outputs:
%      image_out- image containing rectangular bar
%
%      label- label of bar vs background
%
%
% Example:
%
% Notes:
%
% See also: CREATE_GAUSS_BAR
%
% Created: 28-Jan-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

x = repmat(1:col, row, 1);
y = repmat((1:row)', 1, col);
 
orientation = pi*orientation/180;
a = sin(orientation);
b = cos(orientation);
c = -((a*centre_x) + (b*centre_y));
dx = a*x + b*y + c;

label = (dx > -.5) & (dx < .5);
image_out = contrast*(dx < halfwidth);


