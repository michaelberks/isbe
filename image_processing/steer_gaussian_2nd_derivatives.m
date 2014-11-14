function [g2d_responses_a] = steer_gaussian_2nd_derivatives(g2d_responses, theta, num_angles)
%STEER_GAUSSIAN_2ND_DERIVATIVES compute the steered Gaussian 2nd derivative
%response for N equally spaced angles on the circle given raw separable
%Ixx, Iyy and Ixy responses
%   [g2d_responses_a] = steer_gaussian_2nd_derivatives(g2d_responses, num_angles)
%
% Inputs:
%      g2d_responses - 2nd deriv responses as computed by COMPUTE_GAUSSIAN_2ND_DERIVATIVES
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
    
    %Compute the angles equally spaced on the circle
    angles = linspace(0, pi, num_angles+1);
    angles(end) = [];
    
    if iscell(g2d_responses)
        num_scales = length(g2d_responses);
        g2d_responses_a = cell(num_scales, 1);
        
        for i_sc = 1:num_scales
            [r c dim] = size(g2d_responses{i_sc}); %#ok
            g2d_responses_a{i_sc} = zeros(r, c, num_angles);
            
            for i_th = 1:num_angles
                theta = angles(i_th);
                g2d_responses_a{i_sc}(:,:,i_th) = steer_func(g2d_responses{i_sc}, theta);    
            end
        end
        
    else 
        [r c num_scales dim] = size(g2d_responses); %#ok
        g2d_responses_a = zeros(r, c, num_scales, num_angles);

        for i_th = 1:num_angles
            theta = angles(i_th);
            g2d_responses_a(:,:,:,i_th) = steer_func(g2d_responses, theta);    
        end
    end
else
    g2d_responses_a = steer_func(g2d_responses, theta);
end

function steered_responses = steer_func(g2d_responses, theta)

    %I_th = Ixy*sin(2th) + Ixx*cos^2(th) + Iyy*sin^2(th); 
    cc = cos(theta)^2;
    ss = sin(theta)^2;
    s2 = sin(2*theta);
    
%     steered_responses = ...
%         g2d_responses(:,:,:,1)*s2 + ...
%         g2d_responses(:,:,:,2)*cc + ...
%         g2d_responses(:,:,:,3)*ss;  
    
    nd = ndims(g2d_responses);
    theta_vec = cat(nd, s2, cc, ss);
    steered_responses = sum(bsxfun(@times, g2d_responses, theta_vec), nd);