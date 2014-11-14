function [reintegrate discard] = revaluate_distal_candidates(distal_xy, candidates_xy, max_angle)
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

if ~exist('max_angle', 'var') || isempty(max_angle)
    max_angle = pi/4;
end

num_candidates = size(candidates_xy, 1);
num_distal = size(distal_xy,1);

reintegrate = false(num_candidates, 1);
discard = false(num_distal, 1);

for i_can = 1:num_candidates

    tan_vectors = abs(atan(...
        (distal_xy(:,2) - candidates_xy(i_can,2)) ./...
        (distal_xy(:,1) - candidates_xy(i_can,1)) ));
    
    violated = tan_vectors > max_angle;
    
    if sum(~violated & ~discard) < 5
        display('Not enough candidates to process revaluation function');
        continue;
    end
    
    rho1 = fit_local_line(distal_xy(~discard,:), candidates_xy(i_can,1));
    rho2 = fit_local_line([distal_xy(~violated & ~discard,:); candidates_xy(i_can,:)], candidates_xy(i_can,1));
    
    if rho2 < rho1
        reintegrate(i_can) = 1;
        discard(violated) = 1;
    end
        

end

function rho = fit_local_line(distal_xy, candidate_x)

%Find nearest 5 candidates in x-direction
[~, sort_idx] = sort(abs(distal_xy(:,1)-candidate_x));
%rho = abs(corr(distal_xy(sort_idx(1:5),1), distal_xy(sort_idx(1:5),2)));
pp = polyfit(distal_xy(sort_idx(1:5),1), distal_xy(sort_idx(1:5),2), 1);
rho = sum((polyval(pp, distal_xy(sort_idx(1:5),1)) - distal_xy(sort_idx(1:5),2)).^2);
