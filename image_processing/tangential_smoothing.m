function [smoothed_image, orientations, z1, z2] = ...
    tangential_smoothing(image_in, orientations, normals)
%TANGENTIAL_SMOOTHING *Insert a one line summary here*
%   [smoothed_image] = tangential_smoothing(image_in, orientations)
%
% Inputs:
%      image_in - *Insert description of input variable here*
%
%      orientations - *Insert description of input variable here*
%
%
% Outputs:
%      smoothed_image - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 27-Jul-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('normals','var') || isempty(normals)
    normals = 0;
end
if ~exist('orientations','var') || isempty(orientations)
    [~, orientations] = gaussian_1st_derivative_gradient(image_in, 0.5);
    normals = 0;
end

%make meshgrid of (x, y) co-ordinates
[y_lim x_lim] = size(image_in);
xx = repmat(1:x_lim, y_lim, 1);
yy = repmat((1:y_lim)', 1, x_lim);

% Compute sin and cosine of normal orientations (i.e. orientations shifted
% through pi/2)
% If the orientations are already normals we don't need to rotate by pi/2;
if normals
	ori_c = cos(orientations + pi/2);
	ori_s = sin(orientations + pi/2);
else
	ori_c = cos(orientations);
	ori_s = sin(orientations);
end

%linearly interpolate values at normal coordinates
%Note the '-y' here because matrix indexing in 'y' direction is inverse to cartesian
z1 = interp2(image_in, xx + ori_c, yy - ori_s, '*linear');
z2 = interp2(image_in, xx - ori_c, yy + ori_s, '*linear');

%discard any point from image_in that is smaller than either of its two normals 
smoothed_image = (z1 + z2) / 2;
