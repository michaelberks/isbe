function [contour_probs] = compute_contour_probability(contour_xy, im_prob, im_ori, n_probs, lambda)
%COMPUTE_CONTOUR_PROBABILITY *Insert a one line summary here*
%   [] = compute_contour_probability()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 28-Jun-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
n_pts = size(contour_xy,1);

contour_probs = zeros(n_pts, length(lambda));
contour_probs(1,:) = 1;

for i_pt = 2:n_pts
    r_i = round(contour_xy(i_pt,2));
    c_i = round(contour_xy(i_pt,1));
    
    p_i = im_prob(r_i, c_i);
    ori_i = im_ori(r_i, c_i);
    
    contour_probs(i_pt,:) = contour_probs(i_pt-1,:) .* (1 - (1 - p_i).*lambda);
    
    if contour_probs(i_pt) < 1e-6
        break;
    end
end

contour_probs = interp1(1:n_pts, contour_probs, linspace(1, n_pts, n_probs), 'linear');