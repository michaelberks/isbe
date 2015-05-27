function [g2d_responses] = compute_gaussian_2nd_derivatives_d(im, init_scale, num_scales)
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

%Make 2nd order directional filters at inti_scaleal scale
[g,dg,ddg] = gaussian_filters_1d(init_scale);

%pre-allocate output arguments
g2d_responses = cell(num_scales,1);
im = double(im);
padsize = round(5*init_scale);

% x = -1.75:0.5:1.75;
% f = (sin(pi*x) .* sin(pi*x/2) + eps) ./ ((pi^2 * x.^2 / 2) + eps);
% f = f / sum(f);
f = [1 5 8 5 1] / 20;

for i_scale = 1:num_scales
    %pre-allocate array for filter responses
    [r c] = size(im);
    g2d_responses{i_scale} = zeros(r, c, 3);
    
    % Filter the image with separable filters and remove padding
    im = padarray(im, [padsize padsize], 'replicate');
    Ixy = conv2(dg',dg,im,'same');
	Ixx = conv2(g',ddg,im,'same');
	Iyy = conv2(ddg',g,im,'same');
    
    %Filter the image
    g2d_responses{i_scale}(:,:,1) = -Ixy(padsize+(1:r), padsize+(1:c));
    g2d_responses{i_scale}(:,:,2) = Ixx(padsize+(1:r), padsize+(1:c));
    g2d_responses{i_scale}(:,:,3) = Iyy(padsize+(1:r), padsize+(1:c));
    
    %Smooth and downsample the image
    im = imresize(im(padsize+(1:r), padsize+(1:c)), 0.5, 'lanczos2');   
    %im = conv2(im(padsize+(1:r), padsize+(1:c)), f, 'same');
    %im = conv2(im(:,1:2:end), f', 'same');
    %im = im(1:2:end,:);
end

%--------------------------------------------------------------------------
%OLD Function
function old_func(im, init_scale, num_scales)

if nargin < 4
    num_angles = 3;
end

%pre-allocate output arguments
g2d_responses = cell(num_scales,1);
im = double(im);

%Make 2nd order directional filters
filters = make_g2d_filt(init_scale, num_angles);

for i_scale = 1:num_scales
    %pre-allocate array for filter responses
    [r c] = size(im);
    g2d_responses{i_scale} = zeros(r, c, num_angles);
    
    %Filter the image by each oriented filter
    for jj = 1:num_angles
        g2d_responses{i_scale}(:,:,jj) = imfilter(im, filters(:,:,jj), 'same', 'replicate');
    end
    
    %Smooth and downsample the image
    %im = imfilter(im, fspecial('Guassian'), 'same', 'replicate');
    im = imresize(im, 0.5, 'lanczos2');
end


function [filters] = make_g2d_filt(sigma, num_angles)
%aux function to make the gaussian 2nd derivative filters

width = round(6*sigma);
ssq = sigma^2;

[x0 y0] = meshgrid(-width:width, -width:width);

filters = zeros(2*width+1, 2*width+1, num_angles);
for i_scale = 1:num_angles
    theta = (i_scale-1)*pi/num_angles;
    xyt = [cos(theta) -sin(theta); sin(theta) cos(theta)] * [x0(:) y0(:)]';

    xt = reshape(xyt(1,:), 2*width+1, 2*width+1);
    yt = reshape(xyt(2,:), 2*width+1, 2*width+1);

    filters(:,:,i_scale) = ...
        exp(-(xt.*xt)/(2*ssq)) .* exp(-(yt.*yt)/(2*ssq)) .* ((xt.*xt)/ssq - 1) * sqrt(2) / sqrt(pi*ssq*ssq);
end
