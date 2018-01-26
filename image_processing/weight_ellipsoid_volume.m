function [wxyz, exyz] = weight_ellipsoid_volume(v, x, y, z, weights)
%WEIGHT_ELLIPSOID_VOLUME *Insert a one line summary here*
%   [wxyz] = weight_ellipsoid_volume(v, x, y, z)
%
% Inputs:
%      v - *Insert description of input variable here*
%
%      x - *Insert description of input variable here*
%
%      y - *Insert description of input variable here*
%
%      z - *Insert description of input variable here*
%
%
% Outputs:
%      wxyz - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 04-Aug-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
[exyz] = compute_ellipsoid_volume(v, x, y, z);
if ~exist('weights', 'var')
    weights = 1 / mean(exyz(exyz > 0));
end
wxyz = exp(-weights*abs(exyz));