function [g2d_responses] = compute_gaussian_2nd_derivatives(im, scales)
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
g2d_responses = zeros(r, c, length(scales), 3);

padsize = round(5*max(scales));
im = padarray(double(im), [padsize padsize], 'replicate');

for ii = 1:length(scales)
    
    %Make 2nd order directional filters
    [g,dg,ddg] = gaussian_filters_1d(scales(ii));
    
    %Filter the image
    % Filter the image with separable filters and remove padding
    Ixy = conv2(dg',dg,im,'same');
	Ixx = conv2(g',ddg,im,'same');
	Iyy = conv2(ddg',g,im,'same');
    
    %Filter the image
    g2d_responses(:,:,ii,1) = -Ixy(padsize+(1:r), padsize+(1:c));
    g2d_responses(:,:,ii,2) = Ixx(padsize+(1:r), padsize+(1:c));
    g2d_responses(:,:,ii,3) = Iyy(padsize+(1:r), padsize+(1:c));
    
end

% %pre-allocate output arguments
% [r c] = size(im);
% g2d_responses = zeros(r, c, length(scales), num_angles);
% im = double(im);
% 
% for ii = 1:length(scales)
%     
%     %Make 2nd order directional filters
%     filters = make_g2d_filt(scales(ii), num_angles);
%     
%     for jj = 1:num_angles
%     
%         %Filter the image
%         g2d_responses(:,:,ii,jj) = imfilter(im, filters(:,:,jj), 'same', 'replicate');
%     end
%     
% end
% 
% 
% function [filters] = make_g2d_filt(sigma, num_angles)
% %aux function to make the gaussian 2nd derivative filters
% 
% width = round(6*sigma);
% ssq = sigma^2;
% 
% [x0 y0] = meshgrid(-width:width, -width:width);
% 
% filters = zeros(2*width+1, 2*width+1, num_angles);
% for ii = 1:num_angles
%     theta = (ii-1)*pi/num_angles;
%     xyt = [cos(theta) -sin(theta); sin(theta) cos(theta)] * [x0(:) y0(:)]';
% 
%     xt = reshape(xyt(1,:), 2*width+1, 2*width+1);
%     yt = reshape(xyt(2,:), 2*width+1, 2*width+1);
% 
%     filters(:,:,ii) = ...
%         exp(-(xt.*xt)/(2*ssq)) .* exp(-(yt.*yt)/(2*ssq)) .* ((xt.*xt)/ssq - 1) * sqrt(2) / sqrt(pi*ssq*ssq);
% end
