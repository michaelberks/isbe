function [connectivity_map] = make_connectivity_map(prob_image, cx, cy, num_c_pts)
%MAKE_CONNECTIVITY_MAP *Insert a one line summary here*
%   [connectivity_map] = make_connectivity_map(prob_image, num_c_pts)
%
% Inputs:
%      prob_image - *Insert description of input variable here*
%
%      num_c_pts - *Insert description of input variable here*
%
%
% Outputs:
%      connectivity_map - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 12-Dec-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

if ~exist('num_c_pts', 'var') || isempty(num_c_pts)
    num_c_pts = 20;
end

c_pts = linspace(0, 1, num_c_pts);

connectivity_map = zeros(size(prob_image));

for i_c = c_pts
    prob_mask = prob_image >= i_c;
    
    connectivity_mask = bwselect(prob_mask, cx, cy);
    
    connectivity_map(connectivity_mask) = i_c;
end
