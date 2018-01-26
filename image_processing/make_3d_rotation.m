function [R] = make_3d_rotation(theta, rot_axis)
%MAKE_3D_ROTATION *Insert a one line summary here*
%   [R] = make_3d_rotation(theta, u)
%
% Inputs:
%      theta - *Insert description of input variable here*
%
%      u - *Insert description of input variable here*
%
%
% Outputs:
%      R - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 10-Aug-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if nargin == 1
    v_x = [0 -theta(3) theta(2); theta(3) 0 -theta(1); -theta(2) theta(1) 0];
    R = expm(v_x);
    
else
    if norm(rot_axis) ~= 1
        rot_axis = rot_axis / norm(rot_axis);
    end
    c = cos(theta);
    s = sin(theta);
    ux = rot_axis(1);
    uy = rot_axis(2);
    uz = rot_axis(3);
    R = [...
        c + (ux^2)*(1 - c) ux*uy*(1 - c) - uz*s ux*uz*(1 - c) + uy*s;
        uy*ux*(1 - c) + uz*s c + (uy^2)*(1 - c) uy*uz*(1 - c) - ux*s;
        uz*ux*(1 - c) - uy*s uz*uy*(1 - c) + ux*s c + (uz^2)*(1 - c)];
end