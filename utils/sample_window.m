function [sampled_window, centre_value] = sample_window(input_matrix, window_size, row, col, missing_val)
%
% This function samples a square window of size window_size around (row,col) in input_matrix.
% Elements of sampled_window which are off the edge of input_matrix are set to NaN's.
%
% This function assumes that window_size is odd (i.e. there will be a central element in
% sampled_window).
%
% Also return the centre_value of the window
input_class = class(input_matrix);

overlap = (window_size-1)/2;

sampled_window = NaN(window_size, window_size);

if nargin > 4
    sampled_window(:) = missing_val;
end

[rows cols] = size(input_matrix);
top_overlap = max(0, overlap - row + 1);
bottom_overlap = max(0, overlap - (rows - row));
left_overlap = max(0, overlap - col + 1);
right_overlap = max(0, overlap - (cols - col));

from_row = row - floor(overlap - top_overlap);
from_col = col - floor(overlap - left_overlap);
to_row = row + ceil(overlap - bottom_overlap);
to_col = col + ceil(overlap - right_overlap);
sampled_window(1+floor(top_overlap):end-ceil(bottom_overlap), 1+floor(left_overlap):end-ceil(right_overlap)) = ...
    input_matrix(from_row : to_row, from_col : to_col);

eval(['sampled_window = ' input_class '(sampled_window);']);
if nargout > 1
    centre_value = input_matrix(row, col);
end

