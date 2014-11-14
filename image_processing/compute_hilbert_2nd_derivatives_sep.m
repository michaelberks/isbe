function [g2dh_responses] = compute_hilbert_2nd_derivatives_sep(im, scales)
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
g2dh_responses = zeros(r, c, length(scales), 4);

padsize = round(5*max(scales));
im = padarray(double(im), [padsize padsize], 'replicate');

for isigma = 1:length(scales)
    % Create Hilbert transform filters of second order derivative
    filters = gaussian_filters_2d_hilbert_sep(scales(isigma));
    
    H2a = conv2(filters(4, :)', filters(1, :), im, 'same');
    H2b = conv2(filters(3, :)', filters(2, :), im, 'same');
    H2c = conv2(filters(2, :)', filters(3, :), im, 'same');
    H2d = conv2(filters(1, :)', filters(4, :), im, 'same');
    
    %Filter the image
    g2dh_responses(:,:,isigma,1) = H2a(padsize+(1:r), padsize+(1:c));
    g2dh_responses(:,:,isigma,2) = H2b(padsize+(1:r), padsize+(1:c));
    g2dh_responses(:,:,isigma,3) = H2c(padsize+(1:r), padsize+(1:c));
    g2dh_responses(:,:,isigma,4) = H2d(padsize+(1:r), padsize+(1:c));
    
end
