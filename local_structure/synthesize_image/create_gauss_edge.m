function [image_out, label, label_orientation] =...
    create_gauss_edge(halfwidth, contrast, orientation, row, col, centre_x, centre_y)
%CREATE_GAUSS_BAR create an image containing a Gaussian bar
%   [image_out,label] = create_gauss_bar(halfwidth,contrast,orientation,row,col)
%
% Inputs:
%      halfwidth- halfwidth of Gaussian profile at half its maximum height
%
%      contrast- maximum height of (scaled) Gaussian profile
%
%      orientation- orientation of bar in image
%
%      row- number of rows in image
%
%      col- number of columns in image
%
%
% Outputs:
%      image_out- image containing Gaussian bar
%
%      label- label of bar vs background (arbitrarily set to label pixels
%      as belonging to the bar +/- 1 std from the centre line
%
%
% Example:
%
% Notes:
%
% See also: CREATE_RECT_BAR
%
% Created: 28-Jan-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 
if nargin < 7
    centre_x = col/2;
    centre_y = row/2;
end

x = repmat(1:col, row, 1);
y = repmat((1:row)', 1, col);
 
a = sin(pi*orientation/180);
b = cos(pi*orientation/180);
c = -((a*centre_x) + (b*centre_y));
dx = abs(a*x + b*y + c);
        
sigma2 = (halfwidth^2) / log(2);

image_out = contrast*exp(-(dx.^2 / sigma2));

%get the ground truth label
if nargout > 1
    %get the ground truth label for the centre line
    label = dx < .5;
    %label = (1.5 < dx) & (dx < 4);
end
if nargout > 2
    %get the ground truth label for the centre line
    label_orientation = zeros(row, col);
    label_orientation(label) = mod(orientation,180);
end