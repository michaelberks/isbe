function [selected_distal selected_non_distal intermediate_selections candidate_displacements candidate_class_probs candidate_class] = ...
    select_vessels_from_candidates_old(candidate_xy, candidate_scores, rot_struc, varargin)
%SELECT_VESSELS_FROM_CANDIDATES *Insert a one line summary here*
%   [] = select_vessels_from_candidates(varargin)
%
% SELECT_VESSELS_FROM_CANDIDATES uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 20-Nov-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, 0, ...
    'max_candidates', 15,...
    'poly_n', 5,...
    'initial_thresh', 0.3,...
    'class_map', [],...
    'weak_vessel_thresh', 0.3,...
    'strong_vessel_thresh', 0.8,...
    'upper_ydist', -140,...
    'lower_ydist', 90,...
    'angle_discard_thresh', pi/2.5,...
    'angle_keep_thresh', 40,...
    'do_distal_sub', 0, ...
    'do_fill_gaps', 1,...
    'do_final_cull', 1,...
    'vessel_centre', []);
clear varargin;

%Make sure candidates are in descending order of score
if ~issorted(candidate_scores(end:-1:1))
    [candidate_scores cs_i] = sort(candidate_scores, 'descend');
    candidate_xy = candidate_xy(cs_i,:);
end
    
num_candidates = length(candidate_scores);
if nargout > 2
    intermediate_selections = false(num_candidates, 4);
end



%Transform candidates 
candidate_xyr = bsxfun(@plus, bsxfun(@minus, candidate_xy,...
    [rot_struc.ncols rot_struc.nrows]/2)*rot_struc.rot_mat',...
    [rot_struc.ncolsr rot_struc.nrowsr]/2);

do_testing = 1;
if do_testing
   xx = 1:8:rot_struc.ncols;
   yy = (1:8:rot_struc.nrows)';
   grid_x = repmat(xx, length(yy), 1);
   grid_y = repmat(yy, 1, length(xx));
   
   [location_distribution] = build_2d_kernel_distribution(...
       candidate_xy(candidate_scores > args.initial_thresh,:),...
       [grid_x(:) grid_y(:)],...
       candidate_scores(candidate_scores > args.initial_thresh,:), 0);
   location_distribution.D_f = reshape(location_distribution.D_f, size(grid_x));
   %location_distribution.D_f = imresize(location_distribution.D_f, [rot_struc.nrowsr rot_struc.ncolsr]);
   
   [~, y_max] = max(location_distribution.D_f);
   candidate_polyfit = interp1(xx, yy(y_max), candidate_xy(:,1), 'linear');
   candidate_displacements = candidate_xy(:,2) - candidate_polyfit;
   
%    figure; imgray(location_distribution.D_f);
%    plot(xx/8, y_max, 'g');
%    for i_can = 1:num_candidates
%        text(candidate_xyr(i_can,1)/8, candidate_xyr(i_can,2)/8, num2str(candidate_scores(i_can),3), 'color', 'r');
%    end
       
else
 
    %Fit polynomial to transformed centres
    [pp, s, mu] = polyfit(rot_struc.centres_r(:,1), rot_struc.centres_r(:,2), args.poly_n); 

    %Compute y-coord for each x-coord of the transformed candidates, and take
    %the difference between the fitted and actual y values
    candidate_polyfit = polyval(pp, candidate_xyr(:,1), s, mu);
    candidate_displacements = candidate_xyr(:,2) - candidate_polyfit;

    %Use initial selected candidates to determine compute the mean displacement
    num_to_select = min(sum(candidate_scores > args.initial_thresh,1), args.max_candidates);
    offset = median(candidate_displacements(1:num_to_select)); 
    candidate_displacements = candidate_displacements - offset;
end

%Get initial selection of distal row    
if ~isstruct(args.class_map)
    kept = (candidate_scores > args.weak_vessel_thresh) & ...
        ((candidate_displacements > args.upper_ydist) | (candidate_scores > args.strong_vessel_thresh) ) &...
        (candidate_displacements < args.lower_ydist );
    candidate_class_probs = candidate_scores;
    candidate_class = [];
else
    candidate_class = interp2(args.class_map.x, args.class_map.y, args.class_map.post_class,...
        candidate_scores, candidate_displacements, 'linear', 0); 
    candidate_class_probs = interp2(args.class_map.x, args.class_map.y, args.class_map.post_probs,...
        candidate_scores, candidate_displacements, 'linear', 0);
    kept = candidate_class == 1;        
end

if nargout > 2
    %Save this selection in the first column of our output
    intermediate_selections(:,1) = kept;
end

%Now see if there are some strong candidates we should really add back in,
%potentially in place of other weaker candidates included because of
%mis-alignment of the image mid line
if args.do_distal_sub
    potential = ~kept & (candidate_scores > args.strong_vessel_thresh);
    [reintegrate discard] = revaluate_distal_candidates(candidate_xy(kept,:),...
        candidate_xy(potential,:));
    kept(kept) = ~discard;
    kept(potential) = reintegrate;

    if nargout > 2
        %Save this selection in the first column of our output
        intermediate_selections(:,2) = kept;
    end
end

%Next up, see if we can add any weak vessels that appear to lie very
%closely on the distal row
if args.do_fill_gaps
    potential = ~kept & (candidate_scores > args.weak_vessel_thresh);
    [reclaimed] = evaluate_distal_candidates(candidate_xy(kept,:), ...
        candidate_xy(potential,:), args.angle_keep_thresh);
    kept(potential) = reclaimed;

    if nargout > 2
        %Save this selection in the first column of our output
        intermediate_selections(:,3) = kept;
    end
end

%Finally, we discard any vessels that lie on top of one another
if args.do_final_cull
    fell_at_the_last = false(num_candidates,1);
    if any(kept)
        %[distal_idx] = select_distal_candidates_relative_angles(candidate_xy(kept,:), candidate_scores(kept,:), 10, args.angle_discard_thresh, 0);
        [distal_idx] = select_distal_candidates(candidate_xy(kept,:), candidate_class_probs(kept,:), 10, args.angle_discard_thresh, 0, 0);

        distal_idx_tf = true(sum(kept),1);
        distal_idx_tf(distal_idx) = 0;

        fell_at_the_last(kept) = distal_idx_tf;
        kept  = kept & ~fell_at_the_last;
    end 
    if nargout > 2
        %Save this selection in the first column of our output
        intermediate_selections(:,4) = kept;
    end
end

%Select any remaining strong vessel candidates below the lower limit as
%non-distal vessels
if ~isstruct(args.class_map)
    selected_non_distal = ~kept & (candidate_displacements >= args.lower_ydist) & (candidate_scores > 0.8);
else
    selected_non_distal = candidate_class == 2;
end
selected_distal = kept;

% if ~isempty(args.vessel_centre)
%     vessel_centre_xyr = bsxfun(@plus, bsxfun(@minus, [args.vessel_centre.x, args.vessel_centre.y],...
%         [rot_struc.ncols rot_struc.nrows]/2)*rot_struc.rot_mat',...
%         [rot_struc.ncolsr rot_struc.nrowsr]/2);
%     
%     vessel_centre_polyfit = polyval(pp, vessel_centre_xyr(:,1), s, mu);
%     vessel_centre_displacements = vessel_centre_xyr(:,2) - vessel_centre_polyfit;
%     valid_centres = (vessel_centre_displacements > (offset + args.upper_ydist)) &...
%         (vessel_centre_displacements < (offset + 2*args.lower_ydist));
%     
%     vessel_centre_sum(1) = sum(args.vessel_centre.prob(valid_centres));
%     vessel_centre_sum(2) = rot_struc.centres_r(end,1) - rot_struc.centres_r(1,1);
% else
%     vessel_centre_sum = [];
% end
    
    
    
    
