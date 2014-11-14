function [line_strength, orientation, scale] = gaussian_haar_line(im, scales)
%GAUSSIAN_2ND_DERIVATIVE_LINE *Insert a one line summary here*
%   [line_strength, orientation, scale] = gaussian_2nd_derivative_line(im, scales)
%
% Inputs:
%      im - *Insert description of input variable here*
%
%      scales - *Insert description of input variable here*
%
%
% Outputs:
%      line_strength - *Insert description of input variable here*
%
%      orientation - *Insert description of input variable here*
%
%      scale - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 01-Dec-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

%pre-allocate output arguments
im = double(im);
[r c] = size(im);
line_strength = zeros(r, c);
orientation = zeros(r, c);
scale = zeros(r, c);

padsize = 3*max(scales);
im = padarray(im, [padsize padsize], 'replicate');

for ii = 1:length(scales)
    
    %Filter the image
    [haar_responses] = compute_haar_responses(im, scales(ii));
    
    %Extract output and remove padding
    Ixy = -haar_responses(padsize+(1:r), padsize+(1:c),1);
    Ixx_minus_yy = haar_responses(padsize+(1:r), padsize+(1:c),2);
    
    theta_a = atan(2*Ixy ./ Ixx_minus_yy) / 2;
    theta_b = theta_a + pi/2;
    wo_theta_a = filter_output(Ixx, Iyy, Ixy, theta_a);
    wo_theta_b = filter_output(Ixx, Iyy, Ixy, theta_b);
    
    %if wo_theta_a at this scale so far is biggest, swap it into the main
    %line_strength, orientation and scale matrices - note the theta are
    %normal to the line direction, so for theta_a we swap in theta_b and
    %vice-versa
    swap_idx = abs(wo_theta_a) > abs(line_strength);
    line_strength(swap_idx) = wo_theta_a(swap_idx);
    orientation(swap_idx) = theta_b(swap_idx);
    scale(swap_idx) = scales(ii);
    
    %Repeat for theta_b
    swap_idx = abs(wo_theta_b) > abs(line_strength);
    line_strength(swap_idx) = wo_theta_b(swap_idx);
    orientation(swap_idx) = theta_a(swap_idx);
    scale(swap_idx) = scales(ii);
    
end

%wrap orientations back into [-pi pi] range
orientation(orientation > pi) = orientation(orientation > pi) - 2*pi;

function wo_theta = filter_output(Ixx, Iyy, Ixy, theta)
%Compute the output of the filter for an arbitrary angle
cc = cos(theta).^2;
ss = sin(theta).^2;
s2 = sin(2*theta);
wo_theta = Ixx.*cc + Iyy.*ss + Ixy.*s2;
