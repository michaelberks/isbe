function [gabor_responses] = compute_gabor_responses(im, scales, n_angles)
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
[r c] = size(im);
gabor_responses = zeros(r, c, length(scales), n_angles);

padsize = round(5*max(scales));
im = padarray(double(im), [padsize padsize], 'replicate');

for isigma = 1:length(scales)
    % Create Gabor filters
    % Slow but unavoidable :(
    filters = gabor_filters(n_angles, scales(isigma));
    
    for itheta = 1:n_angles
        Ifilt = conv2(im, filters(:,:,itheta), 'same');
        gabor_responses(:,:,isigma,itheta) = Ifilt(padsize+(1:r), padsize+(1:c));
    end
end
