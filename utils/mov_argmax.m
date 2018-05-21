function [argmax_idx, argmax_unique] = mov_argmax(x, win_sz)
%MOV_ARGMAX *Insert a one line summary here*
%   [argmax_idx, argmax_unique] = mov_argmax(x, win_sz)
%
% Inputs:
%      x - *Insert description of input variable here*
%
%      win_sz - *Insert description of input variable here*
%
%
% Outputs:
%      argmax_idx - *Insert description of input variable here*
%
%      argmax_unique - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 26-Apr-2018
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
off_lo = floor(win_sz / 2);
if rem(win_sz, 2)
    off_hi = off_lo;
else
    off_hi = off_lo - 1;
end

n_x = length(x);
argmax_idx = zeros(size(x));
for i_x = 1:n_x
    s_idx = max(1, i_x - off_lo);
    e_idx = min(i_x + off_hi, n_x);
    
    [~, max_idx] = max(x(s_idx:e_idx));
    argmax_idx(i_x) = s_idx + max_idx - 1;
end
argmax_unique = unique(argmax_idx);
