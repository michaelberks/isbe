function [distal_idx] = evaluate_distal_candidates(distal_xy, candidates_xy, dist_thresh)
%SELECT_DISTAL_CANDIDATES *Insert a one line summary here*
%   [distal_idx] = select_distal_candidates(candidates_xy, max_iterations)
%
% Inputs:
%      candidates_xy - *Insert description of input variable here*
%
%      max_iterations - *Insert description of input variable here*
%
%
% Outputs:
%      distal_idx - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 12-Nov-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

if ~exist('dist_thresh', 'var') || isempty(dist_thresh)
    dist_thresh = 20;
end

num_candidates = size(candidates_xy, 1);
distal_idx = false(num_candidates, 1);

num_distal = size(distal_xy,1);

if num_distal < 2
    return;
end

for i_can = 1:num_candidates
    
    [sorted_d, sort_idx] = sort(distal_xy(:,1)-candidates_xy(i_can,1));
    
    i1 = find(sorted_d < 0, 1, 'last');
    
    if i1 == num_distal
        i2 = i1 - 1;
    elseif isempty(i1)
        i1 = 1;
        i2 = 2;
    else
        i2 = i1 + 1;
    end
    
    n1_xy = distal_xy(sort_idx(i1),:);
    n2_xy = distal_xy(sort_idx(i2),:);
    
    n = n2_xy - n1_xy;
    n = n / sqrt(sum(n.^2));
    p_to_a = n1_xy - candidates_xy(i_can,:);
    
    d = p_to_a - dot(p_to_a, n)*n;
    d = sqrt(sum(d.^2));

    distal_idx(i_can) = d < dist_thresh;

end