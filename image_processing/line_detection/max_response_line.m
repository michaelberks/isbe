function [max_response, max_band, max_scale] = max_response_line(responses)
%GAUSSIAN_2ND_DERIVATIVE_LINE *Insert a one line summary here*
%   [max_response, max_band, max_scale] = gaussian_2nd_derivative_line(im, scales)
%
% Inputs:
%      im - *Insert description of input variable here*
%
%      scales - *Insert description of input variable here*
%
%
% Outputs:
%      max_response - *Insert description of input variable here*
%
%      max_band - *Insert description of input variable here*
%
%      max_scale - *Insert description of input variable here*
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

%Get number of scales and image size from input responses
if iscell(responses)
    [r c num_bands] = size(responses{1});
    num_scales = length(responses);
else
    [r c num_scales num_bands] = size(responses);
end

%pre-allocate output arguments 
max_response = zeros(r, c);
max_band = zeros(r, c);
max_scale = zeros(r, c);

for i_scale = 1:num_scales
    
    for i_band = 1:num_bands
        if iscell(responses)
            %Interpolate the responses back to the full pixel grid
            scaling = 2^(i_scale-1);
            offset = (i_scale-1) / (2^i_scale);
            rr = offset + (1:r)' / scaling;
            cc = offset + (1:c) / scaling;
            band_response = interp2(responses{i_scale}(:,:,i_band), cc, rr, 'cubic');
        else
            band_response = responses(:,:,i_scale,i_band);
        end
        
        swap_idx = abs(band_response) > abs(max_response);
        max_response(swap_idx) = band_response(swap_idx);
        max_band(swap_idx) = i_band;
        max_scale(swap_idx) = i_scale;
        
    
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------