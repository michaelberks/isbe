function [input_image] = mb_place_window(input_image, sampled_window, row, col)
%MB_PLACE_WINDOW *Insert a one line summary here*
%   [input_image] = mb_place_window(input_image,sampled_window,row,col)
%
% Inputs:
%      input_image- *Insert description of input variable here*
%
%      sampled_window- *Insert description of input variable here*
%
%      row- *Insert description of input variable here*
%
%      col- *Insert description of input variable here*
%
%
% Outputs:
%      input_image- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 02-Sep-2008
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

window_size = size(sampled_window,1);
window_overlap = (window_size - 1) / 2;

[im_rows im_cols] = size(input_image);

top_overlap = max(0, window_overlap - row + 1);
bottom_overlap = max(0, window_overlap - (im_rows - row));
left_overlap = max(0, window_overlap - col + 1);
right_overlap = max(0, window_overlap - (im_cols - col));

from_row = row - window_overlap + top_overlap;
from_col = col - window_overlap + left_overlap;
to_row = row + window_overlap - bottom_overlap;
to_col = col + window_overlap - right_overlap;

input_image(from_row : to_row, from_col : to_col) = ...
    sampled_window(1+top_overlap:end-bottom_overlap, 1+left_overlap:end-right_overlap);

