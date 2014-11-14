function [line_strength, orientation, scale] = gaussian_2nd_derivative_line2(im, scales)
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
[r c] = size(im);
line_strength = zeros(r, c);
orientation = zeros(r, c);
scale = zeros(r, c);

im = double(im);

for ii = 1:length(scales)
    
    %Make 2nd order directional filters
    [g g1 g2] = make_g_filt(scales(ii));
    
    %Filter the image
    Ixx = conv2(g', g2, im, 'same');
    Iyy = conv2(g2', g, im, 'same');
    Ixy = conv2(-g1', g1, im, 'same');
    
    theta_a = atan(2*Ixy ./ (Ixx - Iyy)) / 2;
    
    theta_b = theta_a + pi/2;
    I_theta_a = filter_output(Ixx, Iyy, Ixy, theta_a);
    I_theta_b = filter_output(Ixx, Iyy, Ixy, theta_b);
    
    %if wo_theta_a at this scale so far is biggest, swap it into the main
    %line_strength, orientation and scale matrices - note the theta are
    %normal to the line direction, so for theta_a we swap in theta_b and
    %vice-versa
    swap_idx = abs(I_theta_a) > abs(line_strength);
    line_strength(swap_idx) = I_theta_a(swap_idx);
    orientation(swap_idx) = theta_b(swap_idx);
    scale(swap_idx) = scales(ii);
    
    %Repeat for theta_b
    swap_idx = abs(I_theta_b) > abs(line_strength);
    line_strength(swap_idx) = I_theta_b(swap_idx);
    orientation(swap_idx) = theta_a(swap_idx);
    scale(swap_idx) = scales(ii);
    
end

%wrap orientations back into [-pi pi] range
orientation(orientation > pi) = orientation(orientation > pi) - 2*pi;

function [g g1 g2] = make_g_filt(sigma)
%aux function to make the gaussian 2nd derivative filters

width = round(6*sigma);
ssq = sigma^2;

x = -width:width;
g = exp(-(x.*x)/(2*ssq)) / sqrt(2*pi*ssq);%
g1 = -x .* g / ssq;
g2 = (x.*x/ssq - 1) .* g / ssq;

function I_theta = filter_output(Ixx, Iyy, Ixy, theta)
%Compute the output of the filter for an arbitrary angle
cc = cos(theta).^2;
ss = sin(theta).^2;
s2 = sin(2*theta);
I_theta = Ixx.*cc + Iyy.*ss + Ixy.*s2;
