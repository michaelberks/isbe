function [g_responses] = compute_gaussian_responses(im, scales)
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
g_responses = zeros(r, c, length(scales));

padsize = round(5*max(scales));
im = padarray(double(im), [padsize padsize], 'replicate');

for ii = 1:length(scales)
    
    %Make 2nd order directional filters
    [g] = gaussian_filters_1d(scales(ii));
    g = g./sum(g);
    
    %Filter the image
    % Filter the image with separable filters and remove padding
    Ig = conv2(g',g,im,'same');
    
    %Filter the image
    g_responses(:,:,ii) = Ig(padsize+(1:r), padsize+(1:c));
    
end