function [h2d_responses_a] = steer_hilbert_2nd_derivatives(h2d_responses, theta, num_angles)
%STEER_GAUSSIAN_2ND_DERIVATIVES compute the steered Gaussian 2nd derivative
%response for N equally spaced angles on the circle given raw separable
%Ixx, Iyy and Ixy responses
%   [h2d_responses_a] = steer_gaussian_2nd_derivatives(h2d_responses, num_angles)
%
% Inputs:
%      h2d_responses - 2nd deriv responses as computed by COMPUTE_GAUSSIAN_2ND_DERIVATIVES
%
%      num_angles - number of angles to steer the reposnes to at each scale
%
%
% Outputs:
%      h2d_responses_a - steered responses equally spaced around the circle
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
    angles = linspace(-pi/2, pi/2, num_angles+1);
    angles(end) = [];
    
    if iscell(h2d_responses)
        num_scales = length(h2d_responses);
        h2d_responses_a = cell(num_scales, 1);
        
        for i_sc = 1:num_scales
            [r c dim] = size(h2d_responses{i_sc}); %#ok
            h2d_responses_a{i_sc} = zeros(r, c, num_angles);
            
            for i_th = 1:num_angles
                theta = angles(i_th);
                h2d_responses_a{i_sc}(:,:,i_th) = steer_func(h2d_responses{i_sc}, theta);    
            end
        end
        
    else 
        [r c num_scales dim] = size(h2d_responses); %#ok
        h2d_responses_a = zeros(r, c, num_scales, num_angles);
        for ii = 1:num_angles
            theta = angles(ii);
            h2d_responses_a(:,:,:,ii) = steer_func(h2d_responses, theta);    
        end
    end
else
    h2d_responses_a = steer_func(h2d_responses, theta-pi/2);
end

function steered_responses = steer_func(h2d_responses, theta)

    ccc = cos(theta)^3;
    sss = sin(theta)^3;
    ccs3 = 3*cos(theta)^2*sin(theta);
    ssc3 = 3*cos(theta)*sin(theta)^2;
    
    %I_th = Ixy*sin(2th) + Ixx*cos^2(th) + Iyy*sin^2(th); 
%     steered_responses = ...
%         h2d_responses(:,:,:,1)*sss + ...
%         h2d_responses(:,:,:,2)*ssc3 + ...
%         h2d_responses(:,:,:,3)*ccs3 + ...
%         h2d_responses(:,:,:,4)*ccc;
    
    nd = ndims(h2d_responses);
    theta_vec = cat(nd, sss, ssc3, ccs3, ccc);
    steered_responses = sum(bsxfun(@times, h2d_responses, theta_vec), nd);
    

