function [pyramid_coeffs] = mb_get_pyramid_coefficients(aPyramid, aLocations)
%
% [pyramid_coeffs] = cjr_get_pyramid_coefficients(aPyramid, aLocations)
%
% CJR_GET_PYRAMID_COEFFICIENTS
%		- get vector(s) of pyramid coeffcients at specified location(s)
%
% Given a pyramid and spatial (row, col) location(s) in that pyramid
% (i.e. position(s) in the image corresponding to the pyramid), this
% function returns the pyramid coefficients at the location for each level
%	of the pyramid as a row vector (or a matrix where each row is the
% coefficients, in the case where multiple locations are specified).
%
%	'aPyramid'
%		- the pyramid.
%
% 'aLocations'
%		- an Nx2 matrix where the first column specifies the row indices
%		and the second column specifies the column indices of the
%		locations in the pyramid we want the coeffcients for. N can be 1
%		if we only need the coefficients for one location.
%
% Return Values:
%
%	CJR_GET_PYRAMID_COEFFICIENTS returns pyramid_coeffs, a matrix of
% the pyramid coefficients where each row corresponds to the coefficients
% at a particular pyramid location.

% work out the number of pyramid levels and orientations
[num_orientations, num_levels] = mb_pyr_info(aPyramid);

% pre-allocate the pyramid_coeffs matrix
pyramid_coeffs = zeros(size(aLocations, 1), (num_orientations * num_levels)+2); % +2 because we also have 1 low- and 1 high-pass band

% loop over the pyramid levels
subband_counter = 1; % this will count across the column of pyramid_coeffs
for level = 1 : num_levels+2 % +2 because we didn't count the high- and low-pass bands 
	for orientation = 1 : num_orientations
		% compute the location for each level
		[r, c] = mb_get_indices_into_pyramid('Function', 'UpPyramid', 'Row', aLocations(:,1), 'Col', aLocations(:,2), 'Level', level);
		
		% place the coefficient for the current level and orientation for each location
		% to be sampled into our data matrix, pyramid_coeffs
		%pyramid_coeffs(:, subband_counter) = diag(aPyramid{level, orientation}(r, c)); % returns a col vector
        idx = sub2ind(size(aPyramid{level,orientation}), r, c);
        pyramid_coeffs(:, subband_counter) = aPyramid{level, orientation}(idx); % returns a col vector
		% diag is used here because we want to be able to cope with multiple indices (r and c can be vectors)
		% -- have a play with a small sample image interactively to see why diag is used
		
		% increment the suband (i.e. column or dimension) counter
		subband_counter = subband_counter + 1;
		
		% break out if doing the very top or very bottom
		if any(level == [1 num_levels+2])
				break;
		end
	end
end
