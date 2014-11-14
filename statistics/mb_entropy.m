function [e] = mb_entropy(p, dim)
%MB_ENTROPY compute entropy of the vector of values in p
%   [e] = mb_entropy(p)
%
% Inputs:
%      p - vector of values to compute entropy
%      dim - computes entropy along specific dimension of p
%
% Outputs:
%      e - entropy of p, defined as -sum(p.*log2(p))
%
%
% Example:
%
% Notes: values of p should be non-negative or a complex result will be
% returned
%
% See also:
%
% Created: 29-Oct-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

%make p sum to 1 (along specified dimension if given)
if nargin < 2
    p = p / sum(p(:));
else
    p = bsxfun(@rdivide, p, sum(p, dim));
end
%remove any nan's introduced by dividing by zero
p(isnan(p)) = 0;

%compute logs, then set values where p equals zero to zero
lp = log2(p);
lp(~p) = 0;

%compute entropy (along specified dimension if given)
if nargin < 2
    e = -sum(p .* lp);
else
    e = -sum(p .* lp, dim);
end

