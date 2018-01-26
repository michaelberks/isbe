function [e2_plane] = compute_ellipse_plane(v, x, y)
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

%Compute ax^2 + by^2 + 2cxy + 2dx + 2ey + f = 0 by substituting
%given z into:
% Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0 where
% v = [A B C ... J];
e2_plane = v(1)*x.^2 + v(2)*y.^2 + 2*v(3)*x.*y + 2*v(4)*x + 2*v(5)*y + v(6);



      


