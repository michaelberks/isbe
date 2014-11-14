function [num_orientations, num_levels] = mb_pyr_info(in_arg)
%
% Get the number of orientations and levels in a pyramid.
% 
% [num_orientations, num_levels] = cjr_pyr_info(pind)
% (for a pyramid in the Simoncelli form)
%
% OR
%
% [num_orientations, num_levels] = cjr_pyr_info(pyramid)
% (for a pyramid in CJR's form)
%
% Get the number of orientations and levels
% in the pyramid described by either the indices in
% pind (as returned by buildSFpyr.m) or in the pyramid stored in the cell array.
% We do not include the high- and low-pass residuals in the number of levels.

if iscell(in_arg)
	num_orientations = size(in_arg, 2);
	num_levels = size(in_arg, 1) - 2; % as we don't count the high- and low-pass bands
else
		% this works by finding all the levels with the same number of
		% pixels as the first band (the high-pass residual), counting them
		% and then subtracting one (i.e. not counting the high-pass residual)
		num_orientations = length(find(in_arg(:,1) == in_arg(1,1)))-1;
		
		% this works by identifying the first different number of rows in
		% the levels and counting them. Try the algorithm with a short
		% list like [1 1 2 2 3 3] to see how it works.
		if nargout == 2
				a = in_arg(:,1);
				b = [0; a];
				a = [a; 0];
				c = a - b;
				c(end) = 0;
				num_levels = length(a(c~=0)) - 1;
		end
end
