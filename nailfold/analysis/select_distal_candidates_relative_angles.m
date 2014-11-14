function [distal_idx] = select_distal_candidates_relative_angles(candidates_xy, candidate_scores, max_iterations, max_angle, do_plot)
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

if ~exist('candidate_scores', 'var')
    candidate_scores = [];
end
if ~exist('max_iterations', 'var') || isempty(max_iterations)
    max_iterations = 1;
end
if ~exist('max_angle', 'var') || isempty(max_angle)
    max_angle = pi/4;
end
if ~exist('do_plot', 'var') || isempty(do_plot)
    do_plot = 0;
end

num_candidates = size(candidates_xy, 1);

violations_mat = true(num_candidates, num_candidates);
tan_mat = zeros(num_candidates, num_candidates);

for i_can = 1:num_candidates

    leave_out_idx = setdiff(1:num_candidates, i_can);
    
    pp = polyfit(candidates_xy(leave_out_idx,1), candidates_xy(leave_out_idx,2), 1);
    rotate_angle = atan(pp(1));
    
    %rot_mat = [cos(rotate_angle) sin(rotate_angle); -sin(rotate_angle) cos(rotate_angle)];
    %candidates_xyr = candidates_xy*rot_mat;

    tan_vectors = abs(atan(...
        (candidates_xy(leave_out_idx,2) - candidates_xy(i_can,2)) ./...
        (candidates_xy(leave_out_idx,1) - candidates_xy(i_can,1)) ) - rotate_angle );

    violations_mat(leave_out_idx, i_can) = tan_vectors < max_angle;
    tan_mat(leave_out_idx,i_can) = tan_vectors;

end

happy_list = 1:num_candidates;
violations_i = violations_mat;
tan_mat_i = tan_mat;

if ~isempty(candidate_scores)
    candidate_scores_i = candidate_scores;
end


for iter = 1:max_iterations

    if do_plot
        figure; hold all; axis equal ij; title(['Iteration: ' num2str(iter)]);
        plot(candidates_xy(happy_list,1), candidates_xy(happy_list,2), 'kx');
    end
    
    while any(~violations_i(:))
        %vio_count = sum(violations_i);
        %[~, worst_offender] = min(vio_count);
        
        if ~isempty(candidate_scores)
            %violations_score = ~violations_i * candidate_scores_i;
            violations_score =  -candidate_scores_i;
            violations_score(all(violations_i)) = -inf;
        else
            violations_score = sum(tan_mat_i);
        end
        
        [~, worst_offender] = max(violations_score);

        violations_i(worst_offender,:) = [];
        violations_i(:,worst_offender) = [];

        tan_mat_i(worst_offender,:) = [];
        tan_mat_i(:,worst_offender) = [];
        
        if ~isempty(candidate_scores)
            candidate_scores_i(worst_offender) = [];
        end
        
        if do_plot
            plot(candidates_xy(happy_list(worst_offender),1), candidates_xy(happy_list(worst_offender),2), 'o');
        end
        happy_list(worst_offender) = [];
    end
    
    readmit_list = find(all(violations_mat(happy_list,:)));
    if isempty(setdiff(readmit_list, happy_list))
        break;
    else
        happy_list = readmit_list;
        violations_i = violations_mat(happy_list, happy_list);
        tan_mat_i = tan_mat(happy_list, happy_list);
        if ~isempty(candidate_scores)
            candidate_scores_i = candidate_scores(happy_list);
        end
    end
end

if do_plot
    display(['Stopped after ' num2str(iter) ' iterations']);
    figure; hold all; axis equal ij;
    plot(candidates_xy(:,1), candidates_xy(:,2), 'kx');
    plot(candidates_xy(happy_list,1), candidates_xy(happy_list,2), 'ro');
end

distal_idx = happy_list;

