function [auto_stats, sequence_names, subject_summary] = analyse_wellcome_apex_measures(varargin)
%ANALYSE_EXTRACTED_APEX_MEASURES *Insert a one line summary here*
%   [] = extract_vessel_centres()
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
% Created: 18-Jun-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
args = u_packargs(varargin,... % the user's input
    0, ... % non-strict mode
    'subject_summary',  [], ...
    'auto_stats',       [], ...
    'sequence_names',   [], ...
    'study_dir',        'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\',...
    'capillary_dir',    'C:\isbe\nailfold\data\wellcome_study\capillary_data\',...
    'selected_features', [],...
    'do_xls',           0,...
    'do_auto_stats',    0,...
    'do_image_plots',   0, ...
    'do_people_plots',  0,...
    'fig_dir',          [],...
    'save_dir',         'C:\isbe\nailfold\data\wellcome_study\results\',...
    'um_per_pix',       1.6077,...
    'xls_wide_format',  1, ...
    'feature_display_names',    [], ...
    'xls_filename_subjects',     'subject_auto_stats.xls');

if isempty(args.subject_summary)    
    %Get subject details from summary file
    fid = fopen([args.study_dir 'subject_summary.txt']);
    frewind(fid);
    s = textscan(fid, '%s'); 
    fclose(fid);
    s = s{1};
    subject_summary = reshape(s, 3, [])';
    [~,s_idx] = sort(subject_summary(:,1));
    subject_summary = subject_summary(s_idx,:);    
else
    subject_summary = args.subject_summary;
end

if isempty(args.sequence_names)  
    num_subs = size(subject_summary,1);
    sequence_names = cell(num_subs, 10);
    
    for i_sub = 1:num_subs
        session_dirs = dir([args.study_dir subject_summary{i_sub,1} '\2015*']);

        num_sessions = length(session_dirs);
        if ~num_sessions
            display(['Subject ' num2str(i_sub) ' has no sessions available.']);
            continue;
        else
            if num_sessions > 1
                display(['Subject ' num2str(i_sub) ' has ' num2str(num_sessions) ' sessions']);
            end
        end
        for i_ses = 1:num_sessions
            curr_seq = 1;
            session_dir = [args.study_dir subject_summary{i_sub,1} '\' session_dirs(i_ses).name '\'];
            for hand = 'LR'
                for digit = '12345'
                    sequence_dirs = ...
                        dir([session_dir hand digit '*']);

                    num_sequences = length(sequence_dirs);
                    if ~num_sequences
                        display(['Subject ' num2str(i_sub) ' has no sequences for ' hand digit]);
                    else
                        sequence_dir = sequence_dirs(end).name;
                        if num_sequences > 1
                            display(['Subject ' num2str(i_sub) ' has ' num2str(num_sequences) ' sequences for ' hand digit '. Using the latest:' sequence_dir]);
                        end
                        if ~exist([session_dir sequence_dir '\sequence_frames_data.dat'], 'file')
                            display(['Sequence ' sequence_dir ' for subject ' num2str(i_sub) ' has no sequence frames data file']);
                        else
                            sequence_names{i_sub,curr_seq} = [session_dir sequence_dir];
                        end
                    end
                    curr_seq = curr_seq+1;
                end
            end
        end
    end
else
    sequence_names = args.sequence_names;
    [num_subs, num_digits] = size(sequence_names);
    if num_digits ~= 10
        error('Sequence names should have ten columns - 1 for each finger and thumb!');
    end
end

if isempty(args.selected_features)
    selected_features = {...
        'num_distal_vessels',...
        'num_nondistal_vessels',...
        'num_giant_vessels',...
        'num_enlarged_vessels',...
        'median_apex_width',...
        'mean_mean_width',...
        'mean_weighted_width',...
        'median_weighted_width',...
        'mean_median_width',...
        'mean_std_width',...
        'max_mean_width',...
        'mean_max_width',...
        'mean_min_width',...
        'max_max_width',...
        'std_mean_width',...
        'mean_orientation_entropy',...
        'median_orientation_entropy',...
        'std_orientation_entropy',...
        'mean_connected_orientation_dispersion',...
        'mean_weighted_orientation_dispersion',...
        'dispersion_base_orientation',...
        'dispersion_connected_orientation',...
        'dispersion_weighted_orientation',...
        'mean_mean_flow',...
        'mean_median_flow',...
        'total_vessel_prob',...
        'mean_vessel_prob',...
          'total_scores',...
           'mean_scores',...
         'score_density',...
        'vessel_density',...
        'mean_inter_capillary_distance',...
        'std_inter_capillary_distance',...
        'median_inter_capillary_distance'};
else
    selected_features = args.selected_features;
end

if isempty(args.auto_stats)
    auto_stats = [];
    args.do_auto_stats = 1;
else
    auto_stats = args.auto_stats;
    args = rmfield(args, 'auto_stats');
end

if args.do_auto_stats

    auto_stats.present = false(num_subs,10);
    auto_stats.gradeable = false(num_subs,10);

    auto_stats.num_distal_vessels = zeros(num_subs,10);
    auto_stats.num_nondistal_vessels = zeros(num_subs,10);
    auto_stats.num_giant_vessels = zeros(num_subs,10);
    auto_stats.num_enlarged_vessels = zeros(num_subs,10);

    auto_stats.median_apex_width = nan(num_subs,10);
    auto_stats.mean_mean_width = nan(num_subs,10);
    auto_stats.mean_weighted_width = nan(num_subs,10);
    auto_stats.median_weighted_width = nan(num_subs,10);
    auto_stats.mean_median_width = nan(num_subs,10);
    auto_stats.mean_std_width = nan(num_subs,10);
    auto_stats.max_mean_width = nan(num_subs,10);
    auto_stats.mean_max_width = nan(num_subs,10);
    auto_stats.max_max_width = nan(num_subs,10);
    auto_stats.mean_min_width = nan(num_subs,10);
    auto_stats.std_mean_width = nan(num_subs,10); 
    
    auto_stats.mean_mean_flow = nan(num_subs,10);
    auto_stats.mean_median_flow = nan(num_subs,10);

    auto_stats.total_vessel_prob = zeros(num_subs,10);
    auto_stats.mean_vessel_prob = zeros(num_subs,10);

    auto_stats.total_scores = zeros(num_subs,10);
    auto_stats.mean_scores = zeros(num_subs,10);

    auto_stats.mean_orientation_entropy = nan(num_subs,10);
    auto_stats.median_orientation_entropy = nan(num_subs,10);
    auto_stats.std_orientation_entropy = nan(num_subs,10);
    
    auto_stats.mean_connected_orientation_dispersion = nan(num_subs,10);
    auto_stats.mean_weighted_orientation_dispersion = nan(num_subs,10);
    
    auto_stats.dispersion_base_orientation = nan(num_subs,10);
    auto_stats.dispersion_connected_orientation = nan(num_subs,10);
    auto_stats.dispersion_weighted_orientation = nan(num_subs,10);

    auto_stats.vessel_density = nan(num_subs,10);

    auto_stats.mean_inter_capillary_distance = nan(num_subs,10);
    auto_stats.std_inter_capillary_distance = nan(num_subs,10);
    auto_stats.median_inter_capillary_distance = nan(num_subs,10);
        
    auto_stats.score_density = nan(num_subs,10);

    for i_sub = 1:num_subs
        for i_digit = 1:10

            auto_stats.present(i_sub, i_digit) = ~isempty(sequence_names{i_sub, i_digit});
            if ~auto_stats.present(i_sub, i_digit)
                continue;
            end
            
            %Get the modified sequence name
            seq_name = sequence_names{i_sub, i_digit};
            dividers = seq_name == '\';
            pos = find(dividers, 3, 'last');
            seq_name(dividers) = '_';
            seq_name = seq_name(pos(1)+1:end);
            capillary_name = [seq_name '_capillary_data.mat'];
        
            %Load apex measures and extract data
            load([args.capillary_dir capillary_name], 'apex_measures');

            auto_stats.gradeable(i_sub, i_digit) = ~isempty(apex_measures.distal);
            
            if ~auto_stats.gradeable(i_sub, i_digit)
                continue;
            end
            
            auto_stats.num_distal_vessels(i_sub, i_digit) = sum(~isnan(apex_measures.distal.mean_weighted_width));
            auto_stats.num_nondistal_vessels(i_sub, i_digit) = sum(~isnan(apex_measures.nondistal.mean_weighted_width));

            if auto_stats.num_distal_vessels(i_sub, i_digit)
                auto_stats.median_apex_width(i_sub, i_digit) = naNmedian(apex_measures.distal.width_at_apex) * args.um_per_pix;
                auto_stats.mean_mean_width(i_sub, i_digit) = naNmean(apex_measures.distal.mean_width) * args.um_per_pix;
                auto_stats.mean_weighted_width(i_sub, i_digit) = naNmean(apex_measures.distal.mean_weighted_width) * args.um_per_pix;
                auto_stats.median_weighted_width(i_sub, i_digit) = naNmedian(apex_measures.distal.mean_weighted_width) * args.um_per_pix;
                auto_stats.mean_median_width(i_sub, i_digit) = naNmean(apex_measures.distal.median_width) * args.um_per_pix;
                auto_stats.mean_std_width(i_sub, i_digit) = naNmean(apex_measures.distal.std_width) * args.um_per_pix;
                auto_stats.max_mean_width(i_sub, i_digit) = nanmax(apex_measures.distal.mean_weighted_width) * args.um_per_pix;          
                auto_stats.mean_max_width(i_sub, i_digit) = naNmean(apex_measures.distal.max_width) * args.um_per_pix;
                auto_stats.max_max_width(i_sub, i_digit) = nanmax(apex_measures.distal.max_width) * args.um_per_pix;
                auto_stats.mean_min_width(i_sub, i_digit) = naNmean(apex_measures.distal.min_width) * args.um_per_pix;
                auto_stats.std_mean_width(i_sub, i_digit) = naNstd(apex_measures.distal.mean_weighted_width) * args.um_per_pix;

                auto_stats.num_giant_vessels(i_sub, i_digit) = sum((apex_measures.distal.mean_weighted_width* args.um_per_pix) > 50);
                auto_stats.num_enlarged_vessels(i_sub, i_digit) = sum((apex_measures.distal.mean_weighted_width* args.um_per_pix) > 20);

%             if auto_stats.num_nondistal_vessels(i_sub, i_digit)
%                 auto_stats.num_giant_vessels(i_sub, i_digit) = auto_stats.num_giant_vessels(i_sub, i_digit) +...
%                     sum((apex_measures.nondistal.mean_weighted_width* args.um_per_pix) > 50);
%                 auto_stats.num_enlarged_vessels(i_sub, i_digit) = auto_stats.num_enlarged_vessels(i_sub, i_digit) +...
%                     sum((apex_measures.nondistal.mean_weighted_width* args.um_per_pix) > 20);
%             end
                auto_stats.num_enlarged_vessels(i_sub, i_digit) = auto_stats.num_enlarged_vessels(i_sub, i_digit) -...
                    auto_stats.num_giant_vessels(i_sub, i_digit);

                ori_entropy = mb_entropy(apex_measures.distal.orientation_hist,2);
                auto_stats.mean_orientation_entropy(i_sub, i_digit) = naNmean(ori_entropy);
                auto_stats.median_orientation_entropy(i_sub, i_digit) = naNmedian(ori_entropy);
                auto_stats.std_orientation_entropy(i_sub, i_digit) = naNstd(ori_entropy);

                auto_stats.mean_connected_orientation_dispersion(i_sub, i_digit) = naNmean(abs(apex_measures.distal.connected_orientation));
                auto_stats.mean_weighted_orientation_dispersion(i_sub, i_digit) = naNmean(abs(apex_measures.distal.weighted_orientation));

                if auto_stats.num_distal_vessels(i_sub, i_digit) > 1
                    auto_stats.dispersion_base_orientation(i_sub, i_digit) = abs(naNmean(exp(1i*angle(apex_measures.distal.base_orientation))));
                    auto_stats.dispersion_connected_orientation(i_sub, i_digit) = abs(naNmean(exp(1i*angle(apex_measures.distal.connected_orientation))));
                    auto_stats.dispersion_weighted_orientation(i_sub, i_digit) = abs(naNmean(exp(1i*angle(apex_measures.distal.weighted_orientation))));
                end
                
                auto_stats.mean_mean_flow(i_sub, i_digit) = naNmean(apex_measures.distal.mean_flow) * args.um_per_pix;
                auto_stats.mean_median_flow(i_sub, i_digit) = naNmean(apex_measures.distal.median_flow) * args.um_per_pix;
                
                auto_stats.total_vessel_prob(i_sub, i_digit) = naNsum(apex_measures.distal.total_prob);
                auto_stats.mean_vessel_prob(i_sub, i_digit) = naNmean(apex_measures.distal.total_prob);

                auto_stats.total_scores(i_sub, i_digit) = naNsum(apex_measures.distal.candidate_scores);
                auto_stats.mean_scores(i_sub, i_digit) = naNmean(apex_measures.distal.candidate_scores);
            
            
                %Compute density measures
                if auto_stats.num_distal_vessels(i_sub, i_digit) > 2
                    apex_xy = sortrows(apex_measures.distal.apex_xy(~isnan(apex_measures.distal.width_at_apex),:));
                    d2 = sqrt(sum((apex_xy(1,:)-apex_xy(end,:)).^2)) * args.um_per_pix / 1000;
                    auto_stats.vessel_density(i_sub, i_digit) = auto_stats.num_distal_vessels(i_sub, i_digit) / d2;
                    auto_stats.score_density(i_sub, i_digit) = auto_stats.total_scores(i_sub, i_digit) / d2;
                    
                    inter_d = sqrt(sum(diff(apex_xy).^2,2)) * args.um_per_pix / 1000;
                    auto_stats.mean_inter_capillary_distance(i_sub, i_digit) = mean(inter_d);
                    auto_stats.std_inter_capillary_distance(i_sub, i_digit) = std(inter_d);
                    auto_stats.median_inter_capillary_distance(i_sub, i_digit) = median(inter_d);
                    
                                                                       
                end
            else
                auto_stats.gradeable(i_sub, i_digit) = 0;
            end
        end
    end
    
    if ~isempty(args.save_dir)
        create_folder(args.save_dir);
        save([args.save_dir '\auto_stats.mat'], 'auto_stats', 'sequence_names', 'subject_summary');
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
if args.do_xls
       
    if args.xls_wide_format
        xls_data = cell(num_subs+1, 10*(length(selected_features)+2) + 3);
        xls_data{1,1} = 'Subject ID';
        xls_data{1,2} = 'Dominant hand';
        xls_data{1,3} = 'Subject group';
        
        xls_data(2:end,1:3) = subject_summary;
        curr_col = 4;
        for i_digit = 1:10
            xls_data{1,curr_col} = 'Present';
            xls_data(2:end,curr_col) = num2cell(auto_stats.present(:,i_digit));
            curr_col = curr_col + 1;
            
            xls_data{1,curr_col} = 'Gradeable';
            xls_data(2:end,curr_col) = num2cell(auto_stats.gradeable(:,i_digit));
            curr_col = curr_col + 1; 
            
            for i_f = 1:length(selected_features);

                feature = selected_features{i_f};
                
                if isempty(args.feature_display_names)
                    feature_str = feature;
                    feature_str(feature_str == '_') = ' ';
                    feature_str(1) = feature_str(1) - 32;
                else
                    feature_str = args.feature_display_names{i_f};
                end

                col_data = auto_stats.(feature)(:,i_digit);
                col_data(~auto_stats.present(:,i_digit)) = -1;
                col_data(~auto_stats.gradeable(:,i_digit)) = -2;
                col_data(isnan(col_data)) = -3;
                
                xls_data{1, curr_col} = feature_str;
                xls_data(2:end, curr_col) = num2cell(col_data);
                curr_col = curr_col + 1;                
            end
        end
        xlswrite([args.save_dir args.xls_filename_subjects], xls_data, 1, 'A2');
    else
        
    end
        
        
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
if args.do_image_plots
    subject_group = str2num(cell2mat(subject_summary(:,3))); %#ok
    subject_group = repmat(subject_group, 1, 10);
    
    for i_f = 1:length(selected_features);

        feature = selected_features{i_f};
        
        if isempty(args.feature_display_names)
            feature_str = feature;
            feature_str(feature_str == '_') = ' ';
            feature_str(1) = feature_str(1) - 32;
        else
            feature_str = args.feature_display_names{i_f};
        end
        
        dist_ss = ...
            auto_stats.(feature)(subject_group >= 2  & auto_stats.gradeable);
        dist_ss(isnan(dist_ss) | ~(dist_ss < inf)) = [];
        dist_pr = ...
            auto_stats.(feature)(subject_group == 1  & auto_stats.gradeable);
        dist_pr(isnan(dist_pr) | ~(dist_pr < inf)) = [];
        dist_hc = ...
            auto_stats.(feature)(subject_group == 0  & auto_stats.gradeable);
        dist_hc(isnan(dist_hc) | ~(dist_hc < inf)) = [];

        p_ss_hc = ranksum(dist_ss, dist_hc);
        p_ss_pr = ranksum(dist_ss, dist_pr);
        p_hc_pr = ranksum(dist_hc, dist_pr);

        if p_ss_hc < 0.0001
            p_ss_hc_str = '< 0.0001';
        else
            p_ss_hc_str = ['= ' num2str(p_ss_hc,2)];
        end
        if p_ss_pr < 0.0001
            p_ss_pr_str = '< 0.0001';
        else
            p_ss_pr_str = ['= ' num2str(p_ss_pr,2)];
        end
        if p_hc_pr < 0.0001
            p_hc_pr_str = '< 0.0001';
        else
            p_hc_pr_str = ['= ' num2str(p_hc_pr,2)];
        end

        display(['Feature: ' feature]);
        display(['SS vs HC: p = ' num2str(p_ss_hc)]);
        display(['SS vs PR: p = ' num2str(p_ss_pr)]);
        display(['PR vs HC: p = ' num2str(p_hc_pr)]);
        display('');

        min_val = min([dist_ss; dist_pr; dist_hc]);
        max_val = max([dist_ss; dist_pr; dist_hc]);
        grid_pts = linspace(min_val, max_val, 100);

        kdist_ss = build_1d_kernel_distribution(dist_ss, grid_pts, 0);
        kdist_pr = build_1d_kernel_distribution(dist_pr, grid_pts, 0);
        kdist_hc = build_1d_kernel_distribution(dist_hc, grid_pts, 0);

        if ~isempty(args.fig_dir)
            figure;
            hold all;
            title({['Feature: \bf' feature_str ', \rm'];...
                ['SS vs HC: p ' p_ss_hc_str ...
                ', SS vs PR: p ' p_ss_pr_str ...
                ', PR vs HC: p ' p_hc_pr_str]},...
                'fontsize', 18);
            set(gca, 'fontSize', 18);
            plot(kdist_ss.x, kdist_ss.D_f, 'linewidth', 2);
            plot(kdist_pr.x, kdist_pr.D_f, 'linewidth', 2);
            plot(kdist_hc.x, kdist_hc.D_f, 'linewidth', 2);
            legend({'SSc', 'PR', 'HC'});
            saveas(gcf, [args.fig_dir feature '_pdf.fig']);
            exportfig([args.fig_dir feature '_pdf.png']);
        else
            figure; 
            subplot(12,2,1:2);
            title({['Feature: \bf' feature_str ', \rm images as independent samples'];...
                ['SS vs HC: p = ' num2str(p_ss_hc,3) ...
                ', SS vs PR: p = ' num2str(p_ss_pr,3) ...
                ', PR vs HC: p = ' num2str(p_hc_pr,3)]});
            axis off;
            subplot(12,2,3:2:23); hold all;
            title('Kernel estimated PDF');
            plot(kdist_ss.x, kdist_ss.D_f, 'linewidth', 2);
            plot(kdist_pr.x, kdist_pr.D_f, 'linewidth', 2);
            plot(kdist_hc.x, kdist_hc.D_f, 'linewidth', 2);
            legend({'SSc', 'PR', 'HC'});

            subplot(12,2,4:2:24); hold all;
            title('Kernel estimated CDF');

            plot(prctile(dist_ss, 0:100), (0:100)/100, 'linewidth', 2);
            plot(prctile(dist_pr, 0:100), (0:100)/100, 'linewidth', 2);
            plot(prctile(dist_hc, 0:100), (0:100)/100, 'linewidth', 2);
            set(gca, 'ylim', [0 1]);
        end
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
if args.do_people_plots
       
    ss_idx = strcmpi(subject_summary(:,3),'2') | strcmpi(subject_summary(:,3),'3');
    pr_idx = strcmpi(subject_summary(:,3),'1');
    hc_idx = strcmpi(subject_summary(:,3),'0');
    n_ss = sum(ss_idx);
    n_pr = sum(pr_idx);
    n_hc = sum(hc_idx);
    for i_f = 1:length(selected_features);

        feature = selected_features{i_f};
        if isempty(args.feature_display_names)
            feature_str = feature;
            feature_str(feature_str == '_') = ' ';
            feature_str(1) = feature_str(1) - 32;
        else
            feature_str = args.feature_display_names{i_f};
        end

        display('***********************');
        display(['Feature: ' feature]);

        dist_f = auto_stats.(feature);
        valid_ims = sum(auto_stats.present & auto_stats.gradeable...
            & ~isnan(dist_f) & (dist_f < inf), 2);

        max_f = max(dist_f, [], 2);
        mean_f = naNsum(dist_f, 2) ./  valid_ims;

        p_ss_hc_max = ranksum(max_f(ss_idx), max_f(hc_idx));
        p_ss_pr_max = ranksum(max_f(ss_idx), max_f(pr_idx));
        p_hc_pr_max = ranksum(max_f(hc_idx), max_f(pr_idx));

        p_ss_hc_avg = ranksum(mean_f(ss_idx), mean_f(hc_idx));
        p_ss_pr_avg = ranksum(mean_f(ss_idx), mean_f(pr_idx));
        p_hc_pr_avg = ranksum(mean_f(hc_idx), mean_f(pr_idx));

%         display(['SS vs HC: p = ' num2str(p_ss_hc_max)]);
%         display(['SS vs PR: p = ' num2str(p_ss_pr_max)]);
%         display(['PR vs HC: p = ' num2str(p_hc_pr_max)]);
%         display('');

%         display(['SS vs HC: p = ' num2str(p_ss_hc_avg)]);
%         display(['SS vs PR: p = ' num2str(p_ss_pr_avg)]);
%         display(['PR vs HC: p = ' num2str(p_hc_pr_avg)]);
%         display('');

        hc_mean = mean(mean_f(hc_idx));
        hc_sse = 1.96*std(mean_f(hc_idx),1) / sqrt(n_hc);
        pr_mean = mean(mean_f(pr_idx));
        pr_sse = 1.96*std(mean_f(pr_idx),1) / sqrt(n_pr);
        ss_mean = mean(mean_f(ss_idx));
        ss_sse = 1.96*std(mean_f(ss_idx),1) / sqrt(n_ss);
        display('Means');
        display(['HC: ' num2str(hc_mean,3) ' (' num2str(hc_mean-hc_sse,3) ', ' num2str(hc_mean+hc_sse,3) ')']);
        display(['PR: ' num2str(pr_mean,3) ' (' num2str(pr_mean-pr_sse,3) ', ' num2str(pr_mean+pr_sse,3) ')']);
        display(['SS: ' num2str(ss_mean,3) ' (' num2str(ss_mean-ss_sse,3) ', ' num2str(ss_mean+ss_sse,3) ')']);
        
        
        if hc_mean < ss_mean
            [~, auc, ~, ~, auc_se] = calculate_roc_curve(mean_f, ss_idx);
        else
            [~, auc, ~, ~, auc_se] = calculate_roc_curve(mean_f, ~ss_idx);
        end
            
        display(['ROC A_z: ' num2str(auc,2) ' (' num2str(auc-1.96*auc_se,2) ', ' num2str(auc+1.96*auc_se,2) ')']);
        continue;

        min_val = min(max_f);
        max_val = max(max_f);
        grid_pts = linspace(min_val, max_val, 100);

        kdist_ss = build_1d_kernel_distribution(max_f(ss_idx), grid_pts, 0);
        kdist_pr = build_1d_kernel_distribution(max_f(pr_idx), grid_pts, 0);
        kdist_hc = build_1d_kernel_distribution(max_f(hc_idx), grid_pts, 0);

        figure; 
        subplot(10,2,1:2);
        title({['Feature: \bf' feature_str];...
            ['\bf Max over images: \rm SS vs HC: p = ' num2str(p_ss_hc_max,3) ...
            ', SS vs PR: p = ' num2str(p_ss_pr_max,3) ...
            ', PR vs HC: p = ' num2str(p_hc_pr_max,3)];...
            ['\bf Mean over images: \rm SS vs HC: p = ' num2str(p_ss_hc_avg,3) ...
            ', SS vs PR: p = ' num2str(p_ss_pr_avg,3) ...
            ', PR vs HC: p = ' num2str(p_hc_pr_avg,3)]});
        axis off;

        subplot(10,2,3:2:9); hold all;
        title('Kernel estimated PDF');
        plot(kdist_ss.x, kdist_ss.D_a, 'linewidth', 2);
        plot(kdist_pr.x, kdist_pr.D_a, 'linewidth', 2);
        plot(kdist_hc.x, kdist_hc.D_a, 'linewidth', 2);
        legend({'SSc', 'PR', 'HC'});
        ylabel('Max of all images per subject');

        subplot(10,2,4:2:10); hold all;
        title('Kernel estimated CDF');
        plot(kdist_ss.x, cumsum(kdist_ss.D_a), 'linewidth', 2);
        plot(kdist_pr.x, cumsum(kdist_pr.D_a), 'linewidth', 2);
        plot(kdist_hc.x, cumsum(kdist_hc.D_a), 'linewidth', 2);
        set(gca, 'ylim', [0 1]);       

        %------------------------------------------------------
        min_val = min(mean_f);
        max_val = max(mean_f);
        grid_pts = linspace(min_val, max_val, 100);

        kdist_ss = build_1d_kernel_distribution(mean_f(ss_idx), grid_pts, 0);
        kdist_pr = build_1d_kernel_distribution(mean_f(pr_idx), grid_pts, 0);
        kdist_hc = build_1d_kernel_distribution(mean_f(hc_idx), grid_pts, 0);

        subplot(10,2,13:2:19); hold all;
        plot(kdist_ss.x, kdist_ss.D_a, 'linewidth', 2);
        plot(kdist_pr.x, kdist_pr.D_a, 'linewidth', 2);
        plot(kdist_hc.x, kdist_hc.D_a, 'linewidth', 2);
        ylabel('Mean of all images per subject');

        subplot(10,2,14:2:20); hold all;
        plot(kdist_ss.x, cumsum(kdist_ss.D_a), 'linewidth', 2);
        plot(kdist_pr.x, cumsum(kdist_pr.D_a), 'linewidth', 2);
        plot(kdist_hc.x, cumsum(kdist_hc.D_a), 'linewidth', 2);
        set(gca, 'ylim', [0 1]);
        
        if ~isempty(args.fig_dir)
            exportfig([args.fig_dir 'subject_' feature '_pdf.png']);
        end

    end

end
