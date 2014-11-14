function [gabor_responses] = compute_gabor_responses_d(im, init_scale, num_scales, n_angles)
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

% pre-allocate output arguments
% Note: responses will be complex

%Make gabor filters at initial scale
filters = gabor_filters(n_angles, init_scale);

%pre-allocate output arguments
gabor_responses = cell(num_scales,1);
im = double(im);
padsize = round(5*init_scale);

for i_scale = 1:num_scales
    %pre-allocate array for filter responses
    [r c] = size(im);
    gabor_responses{i_scale} = zeros(r, c, 3);
    
    % Filter the image with separable filters and remove padding
    im = padarray(im, [padsize padsize], 'replicate');
    for itheta = 1:n_angles
        Ifilt = conv2(im, filters(:,:,itheta), 'same');
        gabor_responses{i_scale}(:,:,itheta) = Ifilt(padsize+(1:r), padsize+(1:c));
    end
    
    %Smooth and downsample the image
    im = imresize(im(padsize+(1:r), padsize+(1:c)), 0.5, 'lanczos2');
end
