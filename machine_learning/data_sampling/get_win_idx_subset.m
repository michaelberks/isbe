function [subset_idx] = get_win_idx_subset(old_win_size, new_win_size)
%GET_WIN_IDX_SUBSET get indices for a new window size that is a subset of
%a larger window
%   [subset_idx] = get_win_idx_subset(old_win_size, new_win_size)
%
% Inputs:
%      old_win_size - 
%
%      new_win_size - 
%
%
% Outputs:
%      subset_idx - indices into old window for new window size;
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

if new_win_size > old_win_size
    error('New window size must be less than or equal to old window size');
    
elseif ~rem(new_win_size, 2)
    error('New window size must be odd');
    
elseif new_win_size == old_win_size
    subset_idx = true(1, old_win_size^2);
    
else
    mask = false(old_win_size);
    centre_pixel = (old_win_size + 1) / 2;
    
    new_half_size = (new_win_size-1) / 2;

    mask(...
        centre_pixel + (-new_half_size:new_half_size),...
        centre_pixel + (-new_half_size:new_half_size) ) = 1;
    
    subset_idx = mask(:)';
end