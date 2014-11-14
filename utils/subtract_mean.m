function [output_mat] = subtract_mean(input_mat, dim)
%SUBTRACT_MEAN *Insert a one line summary here*
%   [output_mat] = subtract_mean(input_mat, dim)
%
% Inputs:
%      input_mat - *Insert description of input variable here*
%
%      dim - dimension upon which to subtract mean
%
%
% Outputs:
%      output_mat - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 25-Mar-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

if nargin < 2
    dim = 1;
end
output_mat = bsxfun(@minus, input_mat, mean(input_mat,dim));
