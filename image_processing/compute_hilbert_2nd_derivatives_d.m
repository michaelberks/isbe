function [g2dh_responses] = compute_hilbert_2nd_derivatives_d(im, init_scale, num_scales)
% COMPUTE_GAUSSIAN_2ND_DERIVATIVES_D Applies Gaussian 2nd derivative filters to an
% image at multiple scales. Multi-scale processing is achieved by successively smoothing
% and down-sampling the original image. See COMPUTE_GAUSSIAN_2ND_DERIVATIVES 
% for a version in which the filters are instead increased
%   [line_strength, orientation, scale] = gaussian_2nd_derivative_line(im, init_scale, num_scales)
%
% Inputs:
%      im - Input image
%
%      init_scale - Sigma of Gaussian 2nd derivative applied at high
%      frequency (i.e. the original image)
%
%      num_scales - number of scales at which to apply the filters
%
%      num_angles - number of oriented filters to aplly at each scale,
%      defaults to 3
%
%
% Outputs:
%      g2d_responses - Responses of the 2nd derivative filters. Returned as
%      a num_scales by 1 cell, where each cell contains the responses of
%      the oriented 2nd derivative filters at that scale and is thus an
%      r by c by num_angles array 
%
%
% Example:
%
% Notes:
%
% See also: COMPUTE_GAUSSIAN_2ND_DERIVATIVES
%
% Created: 16-Feb-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

%Make separable Hilbert filters at initial scale
filters = gaussian_filters_2d_hilbert_sep(init_scale);

%pre-allocate output arguments
g2dh_responses = cell(num_scales,1);
im = double(im);
padsize = round(5*init_scale);

% x = -1.75:0.5:1.75;
% f = (sin(pi*x) .* sin(pi*x/2) + eps) ./ ((pi^2 * x.^2 / 2) + eps);
% f = f / sum(f);
f = [1 5 8 5 1] / 20;

for i_scale = 1:num_scales
    %pre-allocate array for filter responses
    [r c] = size(im);
    g2dh_responses{i_scale} = zeros(r, c, 3);
    
    % Filter the image with separable filters and remove padding
    im = padarray(im, [padsize padsize], 'replicate');

    H2a = conv2(filters(4, :)', filters(1, :), im, 'same');
    H2b = conv2(filters(3, :)', filters(2, :), im, 'same');
    H2c = conv2(filters(2, :)', filters(3, :), im, 'same');
    H2d = conv2(filters(1, :)', filters(4, :), im, 'same');
    
    %Filter the image
    g2dh_responses{i_scale}(:,:,1) = H2a(padsize+(1:r), padsize+(1:c));
    g2dh_responses{i_scale}(:,:,2) = H2b(padsize+(1:r), padsize+(1:c));
    g2dh_responses{i_scale}(:,:,3) = H2c(padsize+(1:r), padsize+(1:c));
    g2dh_responses{i_scale}(:,:,4) = H2d(padsize+(1:r), padsize+(1:c));
    
    %Smooth and downsample the image
    im = imresize(im(padsize+(1:r), padsize+(1:c)), 0.5, 'lanczos2');
    %im = conv2(im(padsize+(1:r), padsize+(1:c)), f, 'same');
    %im = conv2(im(:,1:2:end), f', 'same');
    %im = im(1:2:end,:);
end

