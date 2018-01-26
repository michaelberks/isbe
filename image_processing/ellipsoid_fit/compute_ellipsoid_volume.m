function [exyz] = compute_ellipsoid_volume(v, x, y, z)
%COMPUTE_ELLIPSOID_VOLUME *Insert a one line summary here*
%   [ellipsoid] = compute_ellipsoid_volume(v, x, y, z)
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
%      ellipsoid - *Insert description of input variable here*
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
exyz = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
          2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
          2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z + v(10);
