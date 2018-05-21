function [sub1, sub2, sub3, sub4, sub5] = find_nd(array_in)
%FIND_ND *Insert a one line summary here*
%   [sub1, sub2, sub3, sub4, sub5] = find_nd(array_in)
%
% Inputs:
%      array_in - *Insert description of input variable here*
%
%
% Outputs:
%      sub1 - *Insert description of input variable here*
%
%      sub2 - *Insert description of input variable here*
%
%      sub3 - *Insert description of input variable here*
%
%      sub4 - *Insert description of input variable here*
%
%      sub5 - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 24-Apr-2018
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
N = ndims(array_in);

if nargout > N
    error('More output subscript dimensions (%d) than dimenison in array (%d)',...
        nargout, N);
end
switch N
    case 2
        [sub1, sub2] = find(array_in);
        
    case 3
        idx = find(array_in);
        [sub1, sub2, sub3] = ind2sub(size(array_in), idx);
        
    case 4
        idx = find(array_in);
        [sub1, sub2, sub3, sub4] = ind2sub(size(array_in), idx);
        
    case 5
        idx = find(array_in);
        [sub1, sub2, sub3, sub4, sub5] = ind2sub(size(array_in), idx);
        
    otherwise
        error('find_nd is only implemented for arrays of 5 dimensions or fewer');
end
        
        
