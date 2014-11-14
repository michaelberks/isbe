function [normal_xy, normal_ppx, normal_ppy, dists] = compute_spline_normals(contour_xy)
%COMPUTE_SPLINE_NORMALS *Insert a one line summary here*
%   [normal_xy, normal_pp, dists] = compute_spline_normals(contour_xy)
%
% Inputs:
%      contour_xy - normal_xy(:,1),y coordinates of contour
%
%
% Outputs:
%      normal_xy - Unit normal vectors at each point on the contour
%
%      normal_ppx/y - Piece-wise quadratic to compute normal direction at
%       each contour point
%
%      dists - Distance between each contour point (of which the spline fit
%       is a function
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 21-Mar-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

%Compute dists
dists = cumsum([0; sum(diff(contour_xy).^2,2)]);

%Fit cubic spline to contour 
ppx = spline(dists, contour_xy(:,1));
ppy = spline(dists, contour_xy(:,2));

%Derive piecewise quadratics from the splines to compute the x and y
%gradients
normal_ppx = ppx;
normal_ppx.order = 3;
normal_ppx.coefs = ppx.coefs(:,1:3)*diag([3 2 1]);

normal_ppy = ppy;
normal_ppy.order = 3;
normal_ppy.coefs = ppy.coefs(:,1:3)*diag([3 2 1]);

%Evaluate the gradient functions and compute normals
yg = ppval(normal_ppy, dists);
xg = ppval(normal_ppx, dists);
normal_xy = [-yg ./ sqrt(xg.^2 + yg.^2) xg ./ sqrt(xg.^2 + yg.^2)];
