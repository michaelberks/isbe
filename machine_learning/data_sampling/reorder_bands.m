function [responses n_bands] = reorder_bands(responses, do_max)
%REORDER_BANDS reorders filter response sub-bands so that the maximum
%response is in the first index of the band dimension
%   [data] = reorder_bands(data, shift_or_max)
%
% Inputs:
%      responses - input data, 4d array, [n_samples, n_pixels, n_levels, n_bands]
%
%      do_max - 1 takes only maximum response, 0 circ shifts the data
%      so max is along first index of bands dim
%
%
% Outputs:
%      responses - reordered input data
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 01-Oct-2012
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
[n_rows, n_cols, n_levels, n_bands] = size(responses);

%Reshape responses so the bands form columns of a matrix
responses = reshape(responses, [], n_bands);

%Work which band has maximum absolute response for each sample
[d, max_band] = max(abs(responses), [], 2); clear d; %#ok

%Loop through each band rotating/selecting the max as appropriate
for i_filter = 1:n_bands

    %Get which samples were maximal in this band
    band_idx = max_band==i_filter;

    if do_max
        %Note taking max in this way select maximum absolute response
        %whilst preserving sign
        responses(band_idx,1) = responses(band_idx,i_filter);
    else
        %We must be rotating (we'll assume the user hasn't being silly
        %enough to select do_max and rotate!) so shift such that
        %i_filter occupies the 1st column
        responses(band_idx,:) = circshift(responses(band_idx,:), [0 1-i_filter]);
    end
end
if do_max
    %Throw away the remaining bands
    responses(:,2:n_bands) = [];
    n_bands = 1;
end
%Reshape responses
responses = reshape(responses, n_rows, n_cols, n_levels, n_bands);
