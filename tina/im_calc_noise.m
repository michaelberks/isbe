function [variance, variance_x, variance_y] = im_calc_noise(input_im, roi, mask, plot_output)
%IM_CALC_NOISE estimate noise on an image by analysing variance of 2nd
%deriavtives
%   [variance] = im_calc_noise(input_im, mask, roi)
%
% Inputs:
%      input_im - *Insert description of input variable here*
%
%      roi - *Insert description of input variable here*
%
%
% Outputs:
%      variance - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 29-Jun-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

n_planes = size(input_im,3);
if ~exist('roi', 'var') || isempty(roi)
    roi = repmat([1 size(input_im,1) 1 size(input_im,2)], n_planes, 1);
end

if ~exist('plot_output', 'var')
    plot_output = false;
end

if ~exist('mask', 'var') || isempty(mask)
    mask = true(size(input_im));
end

variance = zeros(n_planes,1);
variance_x = zeros(n_planes,1);
variance_y = zeros(n_planes,1);
for i_p = 1:n_planes
    
    %Create a mask of 2nd derivative magnitudes at a smoother scale -
    %remove these from image so they don't get counted in the noise
    [g2mag g2ori] = gaussian_2nd_derivative_line(input_im(:,:,i_p), [1 2 4]);
    g2nms = mb_non_maximal_supp(abs(g2mag), g2ori);
    g2mask = ~hysterisis(g2nms, [], 0.95);

    %Compute 2nd derivatives in x and y
    dxx = imfilter(input_im(:,:,i_p), [0.5 -1 0.5], 'symmetric');
    dyy = imfilter(input_im(:,:,i_p), [0.5; -1; 0.5], 'symmetric');

    %Discard zero values adjacent to other zero to remove large flat areas from
    %image
    zero_counts_x = conv2(double(~dxx), [1 1 1], 'same');
    valid_x = zero_counts_x < 3;

    zero_counts_y = conv2(double(~dyy), [1 1 1], 'same');
    valid_y = zero_counts_y < 3;

    %Estimate variance from histogram of second derivatives
    [variance_x(i_p,:)] = imf_sigma(dxx, valid_x & mask(:,:,i_p) & g2mask, roi(i_p,:));
    [variance_y(i_p,:)] = imf_sigma(dyy, valid_y & mask(:,:,i_p) & g2mask, roi(i_p,:));
    variance(i_p,:) = (variance_x(i_p,:) + variance_y(i_p,:) ) /2; 

    if plot_output

        [hist_x] = initialise_histogram(dxx(valid_x), [], [], 50);
        [hist_y] = initialise_histogram(dyy(valid_y), [], [], 50);

        figure; 
        subplot(2,3,1); imgray(input_im(:,:,i_p));
        subplot(2,3,2); imgray(dxx);
        subplot(2,3,3); imgray(dyy);
        subplot(2,3,5); bar(hist_x.xbins, hist_x.xcounts);
        title(['Var_x = ' num2str(variance_x(i_p,:),3)]);
        subplot(2,3,6); bar(hist_y.xbins, hist_y.xcounts);
        title(['Var_y = ' num2str(variance_y(i_p,:),3)]);

    end
end