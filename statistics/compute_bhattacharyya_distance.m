function [D_b] = compute_bhattacharyya_distance(p, q)
%COMPUTE_BHATTACHARYYA_DISTANCE *Insert a one line summary here*
%   [D_b] = compute_bhattacharyya_distance(p, q)
%
% Inputs:
%      p - *Insert description of input variable here*
%
%      q - *Insert description of input variable here*
%
%
% Outputs:
%      D_b - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 24-Mar-2014
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
p = p(:) / sum(p(:));
q = q(:) / sum(q(:));

D_b = -log( sum(sqrt(p.*q)) );
