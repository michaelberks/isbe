function [k_map] = complex_gaussian_curvature(theta, sigma)
%COMPLEX_GAUSSIAN_CURVATURE *Insert a one line summary here*
%   [k_map] = complex_gaussian_curvature(theta, sigma)
%
% Inputs:
%      theta - *Insert description of input variable here*
%
%      sigma - *Insert description of input variable here*
%
%
% Outputs:
%      k_map - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 03-May-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

[g, dg] = gaussian_filters_1d(sigma);

if isreal(theta)
    ctheta = complex(cos(2*theta), sin(2*theta));
else
    ctheta = theta;
    theta = angle(ctheta)/2;
end

dtheta_x = real( -.5i*conj(ctheta).*conv2(ctheta, g'*dg/sigma, 'same') );
dtheta_y = real( .5i*conj(ctheta).*conv2(ctheta, dg'*g/sigma, 'same') );

k_map = cos(theta).*dtheta_x + sin(theta).*dtheta_y;
