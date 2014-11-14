function [dec_col] = bin_row_to_dec(bin_array)
%BIN_ROW_TO_DEC *Insert a one line summary here*
%   [dec_col] = bin_row_to_dec(bin_array)
%
% Inputs:
%      bin_array - *Insert description of input variable here*
%
%
% Outputs:
%      dec_col - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 25-Feb-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
bin_array = logical(bin_array);
[nr nc] = size(bin_array);
dec_col = zeros(nr,1);
for cc = 1:nc
    dec_col = dec_col + (2^(cc-1))*bin_array(:,cc);
end
    
