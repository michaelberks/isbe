function [spline_contour_xy, dists, dists_i] = spline_contour(contour_xy, num_pts, spacing, interp_method)
%SPLINE_CONTOUR *Insert a one line summary here*
%   [spline_contour_xy, dists] = spline_contour(contour_xy, num_pts)
%
% Inputs:
%      contour_xy - *Insert description of input variable here*
%
%      num_pts - *Insert description of input variable here*
%
%
% Outputs:
%      spline_contour_xy - *Insert description of input variable here*
%
%      dists - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 22-Mar-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

if nargin < 4
    interp_method = 'spline';
end

%Compute dists
dists = cumsum([0; sqrt(sum(diff(contour_xy).^2,2))]);

if nargin < 2
    num_pts = round(dists(end));
elseif isempty(num_pts)
    num_pts = round(dists(end) / spacing);
end

%Fit spline
dists_i = linspace(0, dists(end), num_pts);
spline_contour_xy = ...
    interp1(dists, contour_xy, dists_i, interp_method);
