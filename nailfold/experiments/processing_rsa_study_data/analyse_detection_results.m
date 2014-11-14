function [im_by_im_counts] = analyse_detection_results(varargin)
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
    'max_missing_markers', inf,...
    'overall_summary', 1,...
    'summary_by_grade', 1,...
    'summary_by_shape', 1,...
    'summary_by_size', 1,...
    'compute_rocs', 1,...
    'calculate_frocs', 1,...
    'analysis_by_width', 1,...
    'fixed_counting', 20,...
    'pct_counting', 0.5,...
    'perfect_counting', 1,...
    'final_selection', 1,...
    'intermediate_selection', 0,...
    'width_bins', 5:5:200,...
    'save_fig_path', []);

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

if args.overall_summary
    total_distal = 0;
    total_non_distal = 0;
    total_undefined = 0;
    %total_confirmed = 0;

    detected_distal = 0;
    detected_non_distal = 0;
    detected_undefined = 0;
    %detected_confirmed = 0;

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
    
end
%% ------------------------------------------------------------------------
% Analyse detections as a function of image type
if args.summary_by_grade
    
    image_grade = cell(num_images,1); 
    for i_im = selected_images
        load([args.apex_gt_dir gt_list(i_im).name], 'majority_grade');
        image_grade{i_im} = majority_grade;
    end

    [image_grade_idx image_grade_labels] = grp2idx(image_grade);
    image_grade_counts = zeros(7,1);

    %
    detection_by_grade = zeros(1,7);
    vessels_by_grade = zeros(1,7);

    for i_im = selected_images
        
        load([args.apex_gt_dir gt_list(i_im).name]);
        
        if args.use_only_gradeable && ~gradeable
            continue;
        end
        i_gr = image_grade_idx(i_im);
        
        image_grade_counts(i_gr) = image_grade_counts(i_gr) + 1;
        
        valid_apexes = ...
            (num_apex_markers >= args.min_num_markers) &...
            (num_im_markers - num_apex_markers <= args.max_missing_markers);
        
        
        detection_by_grade(i_gr) = detection_by_grade(i_gr) + ...
            sum(detections{i_im,2}(is_distal & valid_apexes));

       vessels_by_grade(i_gr) = vessels_by_grade(i_gr) + ...
           sum((is_distal & valid_apexes));
    end
   
    display('*-------------------------------------------------------------*'); 
    display('*-------------- Apex detections by image grade ---------------*');
    for i_gr = 1:7
        display([image_grade_labels{i_gr} ': % apices detected: ' num2str(detection_by_grade(i_gr) / vessels_by_grade(i_gr), 3)...
            ' (' num2str(detection_by_grade(i_gr)) ' of ' num2str(vessels_by_grade(i_gr)) ' in ' num2str(image_grade_counts(i_gr)) ' images)']);
    end
end

%% ------------------------------------------------------------------------
% Analyse detections by vessel shape
if args.summary_by_shape
    %--------------------------------------------------------------------------

shape_labels = {'Non-specific', 'Normal', 'Angiogenic','Meandering'};

    detections_by_shape = zeros(1,4);
    counts_by_shape = zeros(1,4);

    for i_im = selected_images
        
        load([args.apex_gt_dir gt_list(i_im).name]);
        
        if args.use_only_gradeable && ~gradeable
            continue;
        end
        
        valid_apexes = ...
            (num_apex_markers >= args.min_num_markers) &...
            (num_im_markers - num_apex_markers <= args.max_missing_markers) &...
            is_distal;
        
        detected = detections{i_im,2}(valid_apexes);
        distal_shape = apex_shape(valid_apexes);
        for i_ap = 1:length(detected)
            i_sh = find(strcmp(distal_shape{i_ap}, shape_labels)); 
            detections_by_shape(i_sh) = detections_by_shape(i_sh) + detected(i_ap);
            counts_by_shape(i_sh) = counts_by_shape(i_sh) + 1;
        end
    end

    display('*-------------------------------------------------------------*'); 
    display('*-------------- Apex detections by vessel shape --------------*');
    for i_sh = 1:4
        display([shape_labels{i_sh} ': % apices detected by correlation or shape: '...
            num2str(detections_by_shape(i_sh) / counts_by_shape(i_sh), 3)...
            ' (' num2str(detections_by_shape(i_sh)) ' of ' num2str(counts_by_shape(i_sh)) ')']);
    end
    

end

%% ------------------------------------------------------------------------
% Analyse detections by vessel size
if args.summary_by_size
    size_labels = {'Normal', 'Enlarged', 'Giant', 'Irregular', 'Undefined'};

    detections_by_size = zeros(1,5);
    counts_by_size = zeros(1,5);

    for i_im = selected_images
        
        load([args.apex_gt_dir gt_list(i_im).name]);
        
        if args.use_only_gradeable && ~gradeable
            continue;
        end
        
        valid_apexes = ...
            (num_apex_markers >= args.min_num_markers) &...
            (num_im_markers - num_apex_markers <= args.max_missing_markers) &...
            is_distal;
        
        detected = detections{i_im,2}(valid_apexes);
        distal_size = apex_size(valid_apexes);
        
        for i_ap = 1:length(detected)
            i_sz = find(strcmp(distal_size{i_ap}, size_labels));
            detections_by_size(i_sz) = detections_by_size(i_sz) + detected(i_ap);
            counts_by_size(i_sz) = counts_by_size(i_sz) + 1;
        end
    end
    
    display('*-------------------------------------------------------------*'); 
    display('*-------------- Apex detections by vessel size ---------------*');
    for i_sz = 1:5
        display([size_labels{i_sz} ': % apices detected by correlation or shape: '...
            num2str(detections_by_size(i_sz) / counts_by_size(i_sz), 3)...
            ' (' num2str(detections_by_size(i_sz)) ' of ' num2str(counts_by_size(i_sz)) ')']);
    end  
end
%% ------------------------------------------------------------------------

if args.compute_rocs
    if ~args.overall_summary || ~args.summary_by_grade
        display('You must also select to compute the overall and summary by image grades to do the ROC analysis');        
    else
        
        all_candidate_pos = false(total_candidates,1);
        all_candidate_neg = false(total_candidates,1);
        all_candidate_scores = zeros(total_candidates,1);
        all_candidate_grades = zeros(total_candidates,1);
        all_candidate_ranks = zeros(total_candidates,1);

        curr_idx = 0;

        for i_im = selected_images;

            load([args.apex_gt_dir gt_list(i_im).name]);
            
            if args.use_only_gradeable && ~gradeable
                continue;
            end
            
            load([args.candidates_dir candidates_list(i_im).name], 'candidate_scores');
            
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
            all_candidate_grades(idx) = image_grade_idx(i_im);
            all_candidate_ranks(idx) = ranks(end:-1:1);

        end

        legend_txt = cell(8,4);
        image_grade_labels{2} = 'Ungradeable (Q)';
        image_grade_labels{7} = 'Ungradeable (C)';

        fig_h = zeros(4,1);
        ax_h = zeros(4,1);
        for i_f = 1:4
            fig_h(i_f) = figure;%
            ax_h(i_f) = gca;  axis equal; axis([0 1 0 1]); hold all;
        end
        
        im_by_im_counts.roc_pts = zeros(1000,2,8);
        
        for i_gr = 1:7
            
            if args.use_only_gradeable && ismember(i_gr, [2 7]);
                continue;
            end

            grade_idx = all_candidate_grades == i_gr;
            
            if any(grade_idx)
                total_candidates_i = sum(grade_idx);
                total_images_i = image_grade_counts(i_gr);
                
                %Compute the standard ROC for each image grade
                [roc_pts, auc] = ...
                    calculate_roc_curve(all_candidate_scores(grade_idx), all_candidate_pos(grade_idx));
                
                plot(ax_h(1), roc_pts(:,1), roc_pts(:,2), '-', 'linewidth', 2); 
                legend_txt{i_gr,1} = [image_grade_labels{i_gr} ', A_z = ' num2str(auc, 3) ... '\pm' num2str(auc_se,2)
                    ' (' num2str(total_candidates_i) ' candidates)']; %in ' num2str(total_images_i) ' images
                
                %Compute the ROC using modified exclusion criteria
                [roc_pts, auc, ~, fp_counts, ~, op_pts] = ...
                    calculate_roc_curve_exclusions(all_candidate_scores(grade_idx),...
                        all_candidate_pos(grade_idx),all_candidate_neg(grade_idx), 1e3);
            
                plot(ax_h(2), roc_pts(:,1), roc_pts(:,2), '-', 'linewidth', 2); 
                legend_txt{i_gr,2} = [image_grade_labels{i_gr} ', A_z = ' num2str(auc, 3)  ... '\pm' num2str(auc_se,2)
                    ' (' num2str(total_candidates_i) ' candidates)']; % in ' num2str(total_images_i) ' images
                
                [~, equal_error_pt] = min(abs(sum(roc_pts, 2)-1));
                equal_error_thresh = op_pts(equal_error_pt);
                
                display(['EER threshold for ' image_grade_labels{i_gr} ': ' num2str(equal_error_thresh)]);
                
                
                roc_pts_adjusted = roc_pts;
                roc_pts_adjusted(:,2) = roc_pts_adjusted(:,2)*detection_by_grade(i_gr) / vessels_by_grade(i_gr);
                auc_adjusted = trapz(roc_pts_adjusted(:,1), roc_pts_adjusted(:,2));
                
                plot(ax_h(3), roc_pts_adjusted(:,1), roc_pts_adjusted(:,2), '-', 'linewidth', 2); 
                legend_txt{i_gr,3} = [image_grade_labels{i_gr} ', A_z = ' num2str(auc_adjusted, 3)  ... '\pm' num2str(auc_se,2)
                    ' (' num2str(total_candidates_i) ' candidates)']; % in ' num2str(total_images_i) ' images
                
                im_by_im_counts.roc_pts(:,:,i_gr) = roc_pts_adjusted;
                
%                 %Compute the ROC using modified exclusion criteria and
%                 %rankings instead of score (equivalent to asking, what if
%                 %kept the first n candidates)
%                 [roc_pts, auc, ~, ~, auc_se] = ...
%                     calculate_roc_curve_exclusions(all_candidate_ranks(grade_idx),...
%                         all_candidate_pos(grade_idx),all_candidate_neg(grade_idx));
            
                plot(ax_h(4), fp_counts / total_images_i, roc_pts(:,2), '-', 'linewidth', 2); 
                legend_txt{i_gr,4} = [image_grade_labels{i_gr} ' (' num2str(total_candidates_i) ' candidates in ' num2str(total_images_i) ' images)'];
                
            else
                
                for i_f = 1:4
                    plot(ax_h(i_f), -1, -1, '-', 'linewidth', 2); 
                    legend_txt{i_gr,i_f} = [image_grade_labels{i_gr} ', no data' ];
                end
                
            end

        end
        
        if args.use_only_gradeable
            legend_txt([2 7],:) = [];
        end
        
        %Compute the standard ROC
        [roc_pts, auc] = ...
            calculate_roc_curve(all_candidate_scores, all_candidate_pos);
        plot(ax_h(1), roc_pts(:,1), roc_pts(:,2), 'k--', 'linewidth', 3);
        legend_txt{end,1} = ['All grades, A_z = ' num2str(auc, 3)   ...' \pm ' num2str(auc_se,3)
                    ' (' num2str(total_candidates) ' candidates)']; % in ' num2str(total_images) ' images
        
        %Compute the ROC using modified exclusion criteria
        [roc_pts, auc, ~, fp_counts, ~, op_pts] = ...
                calculate_roc_curve_exclusions(all_candidate_scores,...
                    all_candidate_pos, all_candidate_neg, 1e3);

        plot(ax_h(2), roc_pts(:,1), roc_pts(:,2), 'k--', 'linewidth', 3);
        legend_txt{end,2} = ['All grades, A_z = ' num2str(auc, 3)  ... ' \pm ' num2str(auc_se,3)
                    ' (' num2str(total_candidates) ' candidates)']; % in ' num2str(total_images) ' images
                
        [~, equal_error_pt] = min(abs(sum(roc_pts, 2)-1));
                equal_error_thresh = op_pts(equal_error_pt);
                
        th90i = find(roc_pts(:,1) < 0.1, 1, 'last');
        thresh90 = op_pts(th90i);
        
        th95i = find(roc_pts(:,1) < 0.05, 1, 'last');
        thresh95 = op_pts(th95i);
        
        th99i = find(roc_pts(:,1) < 0.01, 1, 'last');
        thresh99 = op_pts(th99i);
                
        display(['EER threshold for all grades: ' num2str(equal_error_thresh)]);
        display(['90% specificity threshold for all grades: ' num2str(thresh90)]);
        display(['95% specificity threshold for all grades: ' num2str(thresh95)]);
        display(['99% specificity threshold for all grades: ' num2str(thresh99)]);
        
        roc_pts_adjusted = roc_pts;
        roc_pts_adjusted(:,2) = roc_pts_adjusted(:,2)*detected_distal / total_distal;
        auc_adjusted = trapz(roc_pts_adjusted(:,1), roc_pts_adjusted(:,2));
        im_by_im_counts.roc_pts(:,:,8) = roc_pts_adjusted;
        
        plot(ax_h(3), roc_pts_adjusted(:,1), roc_pts_adjusted(:,2), 'k--', 'linewidth', 3); 
        legend_txt{end,3} = ['All grades, A_z = ' num2str(auc_adjusted, 3)  ... '\pm' num2str(auc_se,2)
            ' (' num2str(total_candidates_i) ' candidates)']; % in ' num2str(total_images_i) ' images
                
        %Compute the ROC using modified exclusion criteria
%         [roc_pts, auc, ~, ~, auc_se] = ...
%                 calculate_roc_curve_exclusions(all_candidate_ranks,...
%                     all_candidate_pos, all_candidate_neg);

        plot(ax_h(4), fp_counts/total_images, roc_pts(:,2), 'k--', 'linewidth', 3);
        legend_txt{end,4} = ['All grades (' num2str(total_candidates) ' candidates in ' num2str(total_images) ' images)'];
            
        for i_f = 1:4
            title(ax_h(i_f), 'ROC curves of model score');
            xlabel(ax_h(i_f), '1 - specificity', 'fontsize', 14);
            ylabel(ax_h(i_f), 'Sensitivity', 'fontsize', 14);
            set(ax_h, 'fontsize', 14);
            legend(ax_h(i_f), legend_txt(:,i_f), 'location', 'southeast', 'fontsize', 14);
            
            if ~isempty(args.save_fig_path)
                figure(fig_h(i_f));
                set(fig_h(i_f), 'windowstyle', 'normal');
                exportfig([args.save_fig_path 'roc' num2str(i_f) '.png']);
            end
            
        end
        
    end
end

%% ------------------------------------------------------------------------
% 4) Analyse the widths of detected/non detected vessels
if args.analysis_by_width

    if ~args.overall_summary
        display('You must also select to compute the overall summary to do the analysis by width');        
    else
    
        total_undetected = total_distal - detected_distal;
        detected_widths = zeros(detected_distal,1);
        undetected_widths = zeros(total_undetected, 1);

        curr_idx_d = 0;
        curr_idx_m = 0;
        for i_im = selected_images

            load([args.apex_gt_dir gt_list(i_im).name]);

            if args.use_only_gradeable && ~gradeable
                continue;
            end

            valid_apexes = ...
                (num_apex_markers >= args.min_num_markers) &...
                (num_im_markers - num_apex_markers <= args.max_missing_markers) &...
                is_distal;

            detected = detections{i_im,2}(valid_apexes);
            distal_width = apex_widths(valid_apexes);

            if any(detected)
                idx = curr_idx_d + (1:sum(detected));
                curr_idx_d = idx(end);
                detected_widths(idx) = distal_width(detected);
            end
            if any(~detected)
                idx = curr_idx_m + (1:sum(~detected));
                curr_idx_m = idx(end);
                undetected_widths(idx) = distal_width(~detected);
            end
        end
        %
        detected_bin_counts = hist(detected_widths, args.width_bins);
        undetected_bin_counts = hist(undetected_widths, args.width_bins);

        figure; bar(args.width_bins, ...
            [detected_bin_counts/sum(detected_bin_counts);...
            undetected_bin_counts/sum(undetected_bin_counts);...
            (detected_bin_counts+undetected_bin_counts) / sum(detected_bin_counts + undetected_bin_counts)]');
        legend({'Missed', 'Detected', 'All'});
        title('Distribution of widths for detected and missed apices');
        
        figure; hold on;
        plot(args.width_bins, undetected_bin_counts/sum(undetected_bin_counts), 'r', 'linewidth', 2);
        plot(args.width_bins, detected_bin_counts/sum(detected_bin_counts), 'b', 'linewidth', 2);
        plot(args.width_bins, (detected_bin_counts+undetected_bin_counts) / sum(detected_bin_counts + undetected_bin_counts), 'g', 'linewidth', 2);
        legend({'Missed', 'Detected', 'All'});
        title('Distribution of widths for detected and missed apices');

        figure; plot(args.width_bins,  detected_bin_counts./ (detected_bin_counts+undetected_bin_counts));
        title('Detection rate as a function of vessel width');
        ylabel('% of vessels detected');
        xlabel('Vessel width');

        
        num_bins = length(args.width_bins);
        detection_ratios = zeros(num_bins,1);

        for i_width = 1:num_bins
            width_i = args.width_bins(i_width);
            detected_i = sum(detected_widths < width_i);
            undetected_i = sum(undetected_widths < width_i);
            detection_ratios(i_width) = detected_i / (detected_i + undetected_i);    
        end
        figure; hold on;
        plot(args.width_bins, detection_ratios, '-');
        title('Detection rate of all vessels up to a given width');
        ylabel('% of vessels detected');
        xlabel('Vessel width');

    end
end
%% ------------------------------------------------------------------------
% How well would we do if we fixed the number of candidates allowed in each
% image
if args.fixed_counting
    fixed_distal_detections = 0;
    fixed_non_distal_detections = 0;
    fixed_undefined_detections = 0;  
    false_positives = 0;
    
    for i_im = selected_images
        
        
        %Load in the GT structure
        load([args.apex_gt_dir gt_list(i_im).name]);
        
        if args.use_only_gradeable && ~gradeable
            continue;
        end
            
        valid_apexes = ...
                (num_apex_markers >= args.min_num_markers) &...
                (num_im_markers - num_apex_markers <= args.max_missing_markers);
        
        num_valid = min(args.fixed_counting, length(detections{i_im,1}));
        
        detected_idx = detections{i_im,3}(1:num_valid);
        true_detected_idx = ...
            unique(detected_idx(detected_idx > 0));
        
        included_distal_detections = ...
            valid_apexes(true_detected_idx)  & ...
            is_distal(true_detected_idx);
            
        included_non_distal_detections = ...
            valid_apexes(true_detected_idx)  & ...
            is_non_distal(true_detected_idx);
        
        included_undefined_detections = ...
            valid_apexes(true_detected_idx)  & ...
            is_undefined(true_detected_idx);
        
        fixed_distal_detections = fixed_distal_detections + ...
            sum(included_distal_detections);
        
        fixed_non_distal_detections = fixed_non_distal_detections + ...
            sum(included_non_distal_detections);
        
        fixed_undefined_detections = fixed_undefined_detections + ...
            sum(included_undefined_detections);
        
        false_positives = false_positives + sum(~detected_idx);
        
        if num_valid == 1 && ~detected_idx && any(valid_apexes)
            display(['Image ' num2str(i_im)]);
        end

    end
    
    display('*-------------------------------------------------------------*');  
    display(['*------------- If we only use the first ' num2str(args.fixed_counting) ' of candidates.... -------------*']);
    display(['% distal apices detected: ' num2str(fixed_distal_detections / total_distal, 3)]);
    display(['% non-distal apices detected: ' num2str(fixed_non_distal_detections / total_non_distal, 3)]);
    display(['% undefined apices detected: ' num2str(fixed_undefined_detections / total_non_distal, 3)]);
    display(['False positives: ' num2str(false_positives) ' in ' num2str(total_images) ' image']);
    
end
%% ------------------------------------------------------------------------
% How well would we do if we kept only a fixed percentage of the candidates
% in each image
if args.pct_counting
    fixed_distal_detections = 0;
    fixed_non_distal_detections = 0;
    fixed_undefined_detections = 0;  
    false_positives = 0;
    
    for i_im = selected_images
        
        
        %Load in the GT structure
        load([args.apex_gt_dir gt_list(i_im).name]);
        
        if args.use_only_gradeable && ~gradeable
            continue;
        end
            
        valid_apexes = ...
                (num_apex_markers >= args.min_num_markers) &...
                (num_im_markers - num_apex_markers <= args.max_missing_markers);
        
        num_valid = round(length(detections{i_im,1})*args.pct_counting);
        
        detected_idx = detections{i_im,3}(1:num_valid);
        true_detected_idx = ...
            unique(detected_idx(detected_idx > 0));
        
        included_distal_detections = ...
            valid_apexes(true_detected_idx)  & ...
            is_distal(true_detected_idx);
            
        included_non_distal_detections = ...
            valid_apexes(true_detected_idx)  & ...
            is_non_distal(true_detected_idx);
        
        included_undefined_detections = ...
            valid_apexes(true_detected_idx)  & ...
            is_undefined(true_detected_idx);
        
        fixed_distal_detections = fixed_distal_detections + ...
            sum(included_distal_detections);
        
        fixed_non_distal_detections = fixed_non_distal_detections + ...
            sum(included_non_distal_detections);
        
        fixed_undefined_detections = fixed_undefined_detections + ...
            sum(included_undefined_detections);
        
        false_positives = false_positives + sum(~detected_idx);

    end
    
    display('*-------------------------------------------------------------*');  
    display(['*------------- If we only use the first ' num2str(100*args.pct_counting) '% of candidates.... -------------*']);  
    display(['% distal apices detected: ' num2str(fixed_distal_detections / total_distal, 3)]);
    display(['% non-distal apices detected: ' num2str(fixed_non_distal_detections / total_non_distal, 3)]);
    display(['% undefined apices detected: ' num2str(fixed_undefined_detections / total_non_distal, 3)]);
    display(['False positives: ' num2str(false_positives) ' in ' num2str(total_images) ' image']);
    
end
%% ------------------------------------------------------------------------
% How well would we do if we knew how many vessels we were supposed to be
% looking for?
if args.perfect_counting
    fixed_distal_detections = 0;
    fixed_non_distal_detections = 0;
    fixed_undefined_detections = 0;  
    false_positives = 0;
    
    for i_im = selected_images
        
        
        %Load in the GT structure
        load([args.apex_gt_dir gt_list(i_im).name]);
        
        if args.use_only_gradeable && ~gradeable
            continue;
        end
            
        valid_apexes = ...
                (num_apex_markers >= args.min_num_markers) &...
                (num_im_markers - num_apex_markers <= args.max_missing_markers);
        
        num_valid = min(length(valid_apexes), length(detections{i_im,1}));%sum(valid_apexes)
        
        detected_idx = detections{i_im,3}(1:num_valid);
        true_detected_idx = ...
            unique(detected_idx(detected_idx > 0));
        
        included_distal_detections = ...
            valid_apexes(true_detected_idx)  & ...
            is_distal(true_detected_idx);
            
        included_non_distal_detections = ...
            valid_apexes(true_detected_idx)  & ...
            is_non_distal(true_detected_idx);
        
        included_undefined_detections = ...
            valid_apexes(true_detected_idx)  & ...
            is_undefined(true_detected_idx);
        
        fixed_distal_detections = fixed_distal_detections + ...
            sum(included_distal_detections);
        
        fixed_non_distal_detections = fixed_non_distal_detections + ...
            sum(included_non_distal_detections);
        
        fixed_undefined_detections = fixed_undefined_detections + ...
            sum(included_undefined_detections);
        
        false_positives = false_positives + sum(~detected_idx);

    end
    
    display('*-------------------------------------------------------------*');  
    display('*----------------- If we use perfect counting ----------------*');  
    display(['% distal apices detected: ' num2str(fixed_distal_detections / total_distal, 3)]);
    display(['% non-distal apices detected: ' num2str(fixed_non_distal_detections / total_non_distal, 3)]);
    display(['% undefined apices detected: ' num2str(fixed_undefined_detections / total_non_distal, 3)]);
    display(['False positives: ' num2str(false_positives) ' in ' num2str(total_images) ' image']);
    
end
%%
if args.final_selection
    
    im_by_im_counts.all_distal = zeros(num_images, 1);
    im_by_im_counts.all_non_distal = zeros(num_images, 1);
    im_by_im_counts.all_undefined = zeros(num_images, 1);
    im_by_im_counts.valid_distal = zeros(num_images, 1);
    im_by_im_counts.valid_non_distal = zeros(num_images, 1);
    im_by_im_counts.valid_undefined = zeros(num_images, 1);
    im_by_im_counts.detected_distal = zeros(num_images, 1);
    im_by_im_counts.detected_non_distal = zeros(num_images, 1);
    im_by_im_counts.false_positives = zeros(num_images, 1);
    im_by_im_counts.bad_mosaic = false(num_images, 1);
    im_by_im_counts.ungradeable = false(num_images, 1);

    final_distal_total = 0;
    correct_distal_detections = 0;
    non_distal_as_distal = 0;
    undefined_as_distal = 0;
    
    distal_as_non_distal = 0;
    correct_non_distal_detections = 0;
    undefined_as_non_distal = 0;
    
    false_positives = 0;
    missed_distal = 0;
    
    
    
    for i_im = selected_images   
        
        %Load in the GT structure
        load([args.apex_gt_dir gt_list(i_im).name]);
        
        load([args.candidates_dir candidates_list(i_im).name], 'kept', 'non_distal', 'intermediate_selections');
        
        valid_apexes = ...
                (num_apex_markers >= args.min_num_markers) &...
                (num_im_markers - num_apex_markers <= args.max_missing_markers);
            
        im_by_im_counts.all_distal(i_im) = sum(is_distal);
        im_by_im_counts.all_non_distal(i_im) = sum(is_non_distal);
        im_by_im_counts.all_undefined(i_im) = sum(is_undefined);
        im_by_im_counts.valid_distal(i_im) = sum(is_distal & valid_apexes);
        im_by_im_counts.valid_non_distal(i_im) = sum(is_non_distal & valid_apexes);
        im_by_im_counts.valid_undefined(i_im) = sum(is_undefined & valid_apexes);
        im_by_im_counts.ungradeable(i_im) = ~gradeable;
        
        im_by_im_counts.detected_distal(i_im) = sum(kept); %#ok
        im_by_im_counts.detected_non_distal(i_im) = sum(non_distal);
        im_by_im_counts.false_positives(i_im) = sum(~detections{i_im,3}(kept));
        
        if isempty(intermediate_selections) %#ok
            im_by_im_counts.bad_mosaic(i_im) = 1;
            continue;
        end
            
        %-----------------------------------------------
        distal_idx = detections{i_im,3}(kept);
        true_distal_idx = ...
            unique(distal_idx(distal_idx > 0));

        included_distal_detections = ...
            valid_apexes(true_distal_idx)  & ...
            is_distal(true_distal_idx);

        included_non_distal_detections = ...
            valid_apexes(true_distal_idx)  & ...
            is_non_distal(true_distal_idx);

        included_undefined_detections = ...
            valid_apexes(true_distal_idx)  & ...
            is_undefined(true_distal_idx);

        correct_distal_detections = correct_distal_detections + ...
            sum(included_distal_detections);

        non_distal_as_distal = non_distal_as_distal + ...
            sum(included_non_distal_detections);

        undefined_as_distal = undefined_as_distal + ...
            sum(included_undefined_detections);

        false_positives = false_positives + im_by_im_counts.false_positives(i_im);
        %----------------------------------------------------------
        
        non_distal_idx = detections{i_im,3}(non_distal);
        true_non_distal_idx = ...
            unique(non_distal_idx(non_distal_idx > 0));

        included_distal_detections = ...
            valid_apexes(true_non_distal_idx)  & ...
            is_distal(true_non_distal_idx);

        included_non_distal_detections = ...
            valid_apexes(true_non_distal_idx)  & ...
            is_non_distal(true_non_distal_idx);

        included_undefined_detections = ...
            valid_apexes(true_non_distal_idx)  & ...
            is_undefined(true_non_distal_idx);

        distal_as_non_distal = distal_as_non_distal + ...
            sum(included_distal_detections);

        correct_non_distal_detections = correct_non_distal_detections + ...
            sum(included_non_distal_detections);

        undefined_as_non_distal = undefined_as_non_distal + ...
            sum(included_undefined_detections);
        
        %--------------------------------------------------------------
        
        missed_distal = missed_distal + ...
            length(setdiff(find(valid_apexes & is_distal), true_distal_idx));
        final_distal_total = final_distal_total + sum(valid_apexes & is_distal);
        
        %---------------------------------------------------------------

    end
    
    num_ungradeable_images = sum(~im_by_im_counts.bad_mosaic & im_by_im_counts.ungradeable);
    num_bad_images = sum(im_by_im_counts.bad_mosaic);
    
    display('*-------------------------------------------------------------*');  
    display('*------- Summary of the final selection of candidates --------*');  
    display(['Correctly detected distal apexes: ' num2str(correct_distal_detections) ...
        ' out of ' num2str(final_distal_total) ' (' num2str(100*correct_distal_detections/final_distal_total) '%)']);
    
    display(['Non-distal vessels marked as distal: ' num2str(non_distal_as_distal)]);
    display(['Undefined vessels marked as distal:' num2str(undefined_as_distal)]);
    
    display(['Distal vessels incorrectly marked as non-distal: ' num2str(distal_as_non_distal)]);
    display(['Correctly marked non-distal vessels: ' num2str(correct_non_distal_detections)]);
    display(['Undefined vessels marked as non-distal: ' num2str(undefined_as_non_distal)]);
    
    display(['False positives (vessels marked as distal but not marked by any observer in any category) '...
        num2str(false_positives)]);
    
    display(['Distal vessels not detected ' num2str(missed_distal)]);
    
    display(['Num images with badly registered mosaics ' num2str(num_bad_images)]);
    display(['Num ungradeable images ' num2str(num_ungradeable_images)]);
end

%%
if args.intermediate_selection
    
    for i_step = 1:4
        final_distal_total = 0;
        correct_distal_detections = 0;
        non_distal_as_distal = 0;
        undefined_as_distal = 0;

        false_positives = 0;
        missed_distal = 0;

        for i_im = selected_images   

            %Load in the GT structure
            load([args.apex_gt_dir gt_list(i_im).name]);
            if args.use_only_gradeable && ~gradeable
                continue;
            end

            load([args.candidates_dir candidates_list(i_im).name], 'kept', 'intermediate_selections');

            if isempty(intermediate_selections)
                continue;
            end
            
            if i_step < 4
                kept = intermediate_selections(:,i_step); 
            end
            
            valid_apexes = ...
                    (num_apex_markers >= args.min_num_markers) &...
                    (num_im_markers - num_apex_markers <= args.max_missing_markers);

            %-----------------------------------------------
            distal_idx = detections{i_im,3}(kept);
            true_distal_idx = ...
                unique(distal_idx(distal_idx > 0));

            included_distal_detections = ...
                valid_apexes(true_distal_idx)  & ...
                is_distal(true_distal_idx);

            included_non_distal_detections = ...
                valid_apexes(true_distal_idx)  & ...
                is_non_distal(true_distal_idx);

            included_undefined_detections = ...
                valid_apexes(true_distal_idx)  & ...
                is_undefined(true_distal_idx);

            correct_distal_detections = correct_distal_detections + ...
                sum(included_distal_detections);

            non_distal_as_distal = non_distal_as_distal + ...
                sum(included_non_distal_detections);

            undefined_as_distal = undefined_as_distal + ...
                sum(included_undefined_detections);

            false_positives = false_positives + sum(~distal_idx);
            %--------------------------------------------------------------

            missed_distal = missed_distal + ...
                length(setdiff(find(valid_apexes & is_distal), true_distal_idx));
            final_distal_total = final_distal_total + sum(valid_apexes & is_distal);

            %---------------------------------------------------------------

        end

        display('*-------------------------------------------------------------*');  
        display(['*------- Summary of the intermediate selection step ' num2str(i_step) ' --------*']);  
        display(['Correctly detected distal apexes: ' num2str(correct_distal_detections) ...
            ' out of ' num2str(final_distal_total) ' (' num2str(100*correct_distal_detections/final_distal_total) '%)']);

        display(['Non-distal vessels marked as distal: ' num2str(non_distal_as_distal)]);
        display(['Undefined vessels marked as distal:' num2str(undefined_as_distal)]);

        display(['False positives (vessels marked as distal but not marked by any observer in any category) '...
            num2str(false_positives)]);

        display(['Distal vessels not detected ' num2str(missed_distal)]);
    end
end
            
        
        
        
        
        




