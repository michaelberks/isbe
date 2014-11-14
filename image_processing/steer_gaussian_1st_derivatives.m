function [g1d_responses_a] = steer_gaussian_1st_derivatives(g1d_responses, theta, num_angles)
%STEER_GAUSSIAN_2ND_DERIVATIVES compute the steered Gaussian 2nd derivative
%response for N equally spaced angles on the circle given raw separable
%Ixx, Iyy and Ixy responses
%   [g2d_responses_a] = steer_gaussian_2nd_derivatives(g1d_responses, num_angles)
%
% Inputs:
%      g1d_responses - 2nd deriv responses as computed by COMPUTE_GAUSSIAN_2ND_DERIVATIVES
%
%      num_angles - number of angles to steer the reposnes to at each scale
%
%
% Outputs:
%      g2d_responses_a - steered responses equally spaced around the circle
%
%
%
% Example:
%
% Notes:
%
% See also: COMPUTE_GAUSSIAN_2ND_DERIVATIVES
%
% Created: 01-Dec-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

%pre-allocate output arguments
if isempty(theta)
    [r c num_scales dim] = size(g1d_responses);
    g1d_responses_a = zeros(r, c, num_scales, num_angles);

    %Compute the angles equally spaced on the circle
    angles = linspace(0, pi, num_angles+1);
    angles(end) = [];

    for ii = 1:num_angles
        theta = -angles(ii);

        %I_th = Ixy*sin(2th) + Ixx*cos^2(th) + Iyy*sin^2(th); 
        g1d_responses_a(:,:,:,ii) = ...
            g1d_responses(:,:,:,1)*cos(theta) + g1d_responses(:,:,:,2)*sin(theta);    
    end
else
    g1d_responses_a = ...
            g1d_responses(:,:,:,1)*cos(-theta) + g1d_responses(:,:,:,2)*sin(-theta);
end
    