function [out_stats] = analyse_incremental_rocs(varargin)
%ANALYSE_DETECTION_RESULTS *Insert a one line summary here*
%   [] = make_detection_results_struc(varargin)
%
% ANALYSE_DETECTION_RESULTS uses the U_PACKARGS interface function
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
% Created: 31-Oct-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
args = u_packargs(varargin, 0,...
    {...
    'results_name',...
    'candidates_dir'    },... %mandatory arguments
    'apex_gt_dir', [nailfoldroot 'data/rsa_study/test/apex_gt/'],...
    'results_dir', [nailfoldroot 'data/rsa_study/test/results/'],...
    'selected_images', [],...
    'use_only_gradeable', 1,...
    'min_num_markers', 1,...
    'max_missing_markers', inf);

detections = [];
gt_list = [];
candidates_list = [];
load([args.results_dir args.results_name '.mat']);

num_images = length(detections);
if isempty(args.selected_images)
    selected_images = 1:num_images;
else
    selected_images = args.selected_images(:)';
end

%% ------------------------------------------------------------------------

total_distal = 0;
total_non_distal = 0;
total_undefined = 0;

detected_distal = 0;
detected_non_distal = 0;
detected_undefined = 0;

total_candidates = 0;
correct_candidates = 0;
    
total_images = 0;
total_images_with_valid_apexes = 0;

for i_im = selected_images


    %Load in the GT structure
    load([args.apex_gt_dir gt_list(i_im).name]);

    if args.use_only_gradeable && ~gradeable
        continue;
    end
    total_images = total_images + 1;

    valid_apexes = ...
        (num_apex_markers >= args.min_num_markers) &...
        (num_im_markers - num_apex_markers <= args.max_missing_markers);

    detected_distal = detected_distal + sum(detections{i_im,2}(is_distal & valid_apexes));
    detected_non_distal = detected_non_distal + sum(detections{i_im,2}(is_non_distal & valid_apexes));
    detected_undefined = detected_undefined + sum(detections{i_im,2}(is_undefined & valid_apexes));

    total_distal = total_distal + sum(is_distal & valid_apexes);
    total_non_distal = total_non_distal + sum(is_non_distal & valid_apexes);
    total_undefined = total_undefined + sum(is_undefined & valid_apexes);

    total_candidates = total_candidates + length(detections{i_im,1});
    correct_candidates = correct_candidates + sum(detections{i_im,1});

    if any(valid_apexes)
        total_images_with_valid_apexes = total_images_with_valid_apexes + 1;
    end
end
    
display('*-------------------------------------------------------------*');  
display('*-------------- Apex detections overall summary --------------*');  
display(['% distal apices detected: ' num2str(detected_distal / total_distal, 3)...
    ' (' num2str(detected_distal) ' of ' num2str(total_distal) ' in ' num2str(total_images) ' images)']);
display(['% non-distal apices detected : ' num2str(detected_non_distal / total_non_distal, 3)...
    ' (' num2str(detected_non_distal) ' of ' num2str(total_non_distal) ' in ' num2str(total_images) ' images)']);
display(['% undefined apices detected: ' num2str(detected_undefined / total_undefined, 3)...
    ' (' num2str(detected_undefined) ' of ' num2str(total_undefined) ' in ' num2str(total_images) ' images)']);
display(['% Candidates that were apices: ' num2str(correct_candidates / total_candidates, 3)...
    ' (' num2str(correct_candidates) ' of ' num2str(total_candidates) ' in ' num2str(total_images) ' images)']);
display(['# Images with valid apices: ' num2str(total_images_with_valid_apexes)]);
    


all_candidate_pos = false(total_candidates,1);
all_candidate_neg = false(total_candidates,1);
all_candidate_scores = zeros(total_candidates,1);

curr_idx = 0;

for i_im = selected_images;

    load([args.apex_gt_dir gt_list(i_im).name]);

    if args.use_only_gradeable && ~gradeable
        continue;
    end

    load([args.candidates_dir candidates_list(i_im).name], 'candidate_scores');

    if isempty(candidate_scores)
        continue;
    end
    
    ranks = 1:length(candidate_scores);         
    idx = curr_idx + ranks;
    curr_idx = idx(end);

    valid_apexes = ...
        (num_apex_markers >= args.min_num_markers) &...
        (num_im_markers - num_apex_markers <= args.max_missing_markers) &...
        is_distal;

    included_detections = detections{i_im,1};
    included_detections(included_detections) = ...
        valid_apexes(detections{i_im,3}(included_detections));


    all_candidate_pos(idx) = included_detections;
    all_candidate_neg(idx) = ~detections{i_im,1};
    all_candidate_scores(idx) = candidate_scores;

end

%Compute the standard ROC
[roc_pts_orig, auc_orig] = ...
    calculate_roc_curve(all_candidate_scores, all_candidate_pos);

%Compute the ROC using modified exclusion criteria
[roc_pts, auc, ~, fp_counts, ~, op_pts] = ...
        calculate_roc_curve_exclusions(all_candidate_scores,...
            all_candidate_pos, all_candidate_neg, 1e3);
 
roc_pts_adjusted = roc_pts;
roc_pts_adjusted(:,2) = roc_pts_adjusted(:,2)*detected_distal / total_distal;
auc_adjusted = trapz(roc_pts_adjusted(:,1), roc_pts_adjusted(:,2));        

[~, equal_error_pt] = min(abs(sum(roc_pts, 2)-1));
        out_stats.equal_error_thresh = op_pts(equal_error_pt);

th90i = find(roc_pts(:,1) < 0.1, 1, 'last');
out_stats.thresh90 = op_pts(th90i);

th95i = find(roc_pts(:,1) < 0.05, 1, 'last');
out_stats.thresh95 = op_pts(th95i);

th99i = find(roc_pts(:,1) < 0.01, 1, 'last');
out_stats.thresh99 = op_pts(th99i);

if args.plot
    figure; hold all;
    legend_txt = cell(3,1);

    plot(roc_pts_orig(:,1), roc_pts_orig(:,2), '-', 'linewidth', 3);
    legend_txt{1} = ['Original, A_z = ' num2str(auc_orig, 3)   ...' \pm ' num2str(auc_se,3)
                ' (' num2str(total_candidates) ' candidates)']; % in ' num2str(total_images) ' images



    plot(roc_pts(:,1), roc_pts(:,2), '-', 'linewidth', 3);
    legend_txt{2} = ['Modified, A_z = ' num2str(auc, 3)  ... ' \pm ' num2str(auc_se,3)
                ' (' num2str(total_candidates) ' candidates)']; % in ' num2str(total_images) ' images

    plot(roc_pts_adjusted(:,1), roc_pts_adjusted(:,2), '-', 'linewidth', 3); 
    legend_txt{3} = ['Modified adjusted, A_z = ' num2str(auc_adjusted, 3)  ... '\pm' num2str(auc_se,2)
        ' (' num2str(total_candidates) ' candidates)']; % in ' num2str(total_images_i) ' images

    title('ROC curves of model score');
    xlabel('1 - specificity', 'fontsize', 14);
    ylabel('Sensitivity', 'fontsize', 14);
    set(gca, 'fontsize', 14);
    legend(legend_txt, 'location', 'southeast', 'fontsize', 14);
    axis([0 1 0 1]);
end
        
out_stats.roc_pts_orig = roc_pts_orig;
out_stats.roc_pts = roc_pts;
out_stats.roc_pts_adjusted = roc_pts_adjusted;
out_stats.fp_counts = fp_counts;

out_stats.total_distal = total_distal;
out_stats.total_non_distal = total_non_distal;
out_stats.total_undefined = total_undefined;

out_stats.detected_distal = detected_distal;
out_stats.detected_non_distal = detected_non_distal;
out_stats.detected_undefined = detected_undefined;

out_stats.total_candidates = total_candidates;
out_stats.correct_candidates = correct_candidates;
    
out_stats.total_images = total_images;
out_stats.total_images_with_valid_apexes = total_images_with_valid_apexes;        
        
        




