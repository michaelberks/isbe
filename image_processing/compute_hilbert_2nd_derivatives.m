function [g2dh_responses] = compute_hilbert_2nd_derivatives(im, scales)
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

n_angles = 4;

%pre-allocate output arguments
[r c] = size(im);
g2dh_responses = zeros(r, c, length(scales), n_angles);

padsize = round(5*max(scales));
im = padarray(double(im), [padsize padsize], 'replicate');

for isigma = 1:length(scales)
    % Create Hilbert transform filters of second order derivative
    % Slow but unavoidable :(
    h_filters = gaussian_filters_2d_hilbert(n_angles, scales(isigma));
    
    for itheta = 1:n_angles
        Ifilt = conv2(im, h_filters(:,:,itheta), 'same');
        g2dh_responses(:,:,isigma,itheta) = Ifilt(padsize+(1:r), padsize+(1:c));
    end
end
