function [local_max_idx] = local_argmax(x, win_sz, min_val)
%LOCAL_ARGMAX *Insert a one line summary here*
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
if ~exist('win_sz', 'var') 
    win_sz = [];
end
if ~exist('min_val', 'var') || isempty(min_val)
    min_val = -inf;
end
x = x(:);
n_x = length(x);

%Get list of all points bigger than their neighbours and bigger than the
%min val
gl = [false; x(2:end) > x(1:end-1)];
gr = [x(2:end) < x(1:end-1); false];
local_max_idx = find(gl & gr & (x > min_val));

n_max = length(local_max_idx);

%If we've got no local max, or we don't need to exclude in a window, we can
%return
if isempty(win_sz) || ~n_max
    return;
end
discard_idx = false(n_max,1);

off_lo = floor(win_sz / 2);
if rem(win_sz, 2)
    off_hi = off_lo;
else
    off_hi = off_lo - 1;
end

for i_m = 1:n_max
    i_x = local_max_idx(i_m);
    s_idx = max(1, i_x - off_lo);
    e_idx = min(i_x + off_hi, n_x);
    
    discard_idx(i_m) = any(x(s_idx:e_idx) > x(i_x));
end
local_max_idx(discard_idx)= [];








