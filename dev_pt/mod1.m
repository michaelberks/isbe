function n_out = mod1(n_in,k)
% MOD1 Does the same as mod except it copies to the range 1..n instead of
% 0..n-1
%
% e.g. mod1(5,3) = 2
%   1 2 3 4 5 6 7 8 ...
%   1 2 3 1 2 3 1 2 ...
n_out = mod(n_in-1,k)+1;
