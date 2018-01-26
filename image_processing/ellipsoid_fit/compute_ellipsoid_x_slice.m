function [v2] = compute_ellipsoid_x_slice(v, x)
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

%Compute terms for ax^2 + by^2 + 2cxy + 2dx + 2ey + f = 0 by substituting
%given z into:
% Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0 where
% v = [A B C ... J];

%a = B, b = C - shift dims y->x, z->x
a = v(2);
b = v(3);

%c = F, yz -> xy
c = v(6);

%d = Dx + H, e = Ex + I
d = v(4)*x + v(8);
e = v(5)*x + v(9);

%f = Ax^2 + 2Gx + J
f = v(1)*x*x + 2*v(7)*x + v(10);

v2 = [a b c d e f]';

% normalize to the more conventional form with constant term = -1
if abs( v2(end) ) > 1e-6
    v2 = -v2 / v2(end); 
else
    v2 = -sign( v2(end) ) * v2;
end



      


