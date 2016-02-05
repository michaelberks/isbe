function [co_occurrence detected_widths missed_widths a_b_widths a_c_widths b_c_widths] = cxx_results_cooc_fun(varargin)
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

args = u_packargs(varargin,... % the user's input
    0, ... % non-strict mode
    {'image_names'},         ...
    'data_dir',             [nailfoldroot 'data/rsa_study/master_set/'],...
    'prob_dir',             'rf_classification/296655',...
    'cluster_dir',          'apex_clusters_merged',...
    'candidates_dir',       'vessel_hogs',...
    'metrics_dir',          [],...'apex_maps\set12g_half_296655\miccai_maxima\apex_metrics'
    'width_feature',        'width_at_apex',...
    'width_fudge',          0,...
    'do_distal',            1,...
    'do_nondistal',         1,...
    'use_post_processing',  1,...
    'plot', 0);

%Build up the co-occurence matrix for the observers
load('C:\isbe\nailfold\data\rsa_study\data_lists\image_id_data.mat');
load('C:\isbe\nailfold\data\rsa_study\data_lists\miccai_lists.mat', 'miccai_selection');


prob_dir = [args.data_dir 'predictions/detection/' args.prob_dir '/'];
cluster_dir = [args.data_dir '/' args.cluster_dir '/'];
candidates_dir = [args.data_dir '/' args.candidates_dir '/'];

if ~isempty(args.metrics_dir)
    metrics_dir = [args.data_dir '/' args.metrics_dir '/'];
end



im_names = args.image_names;

num_images = length(im_names);

g = gaussian_filters_1d(2);
g = g / sum(g); 
    

co_occurrence = zeros(3,3,3, num_images);
a_b_widths = cell(num_images,1);
a_c_widths = cell(num_images,1);
b_c_widths = cell(num_images,1);
detected_widths = zeros(0,1);
missed_widths = zeros(0,1);

for i_im = 1:num_images
        
    im_name = im_names{i_im};
    load([cluster_dir im_name '_apex_clusters.mat'], 'vessels');

    %load in selected candidates
    candidates_data = load([candidates_dir im_name '_apex_candidates.txt'], 'candidate_xy');
    candidate_xy = candidates_data(:,1:2) + 1;
    %selected_distal = candidates_data(:,11) > -1;%0;
    %selected_non_distal = candidates_data(:,12) > 10;%0;
    if args.use_post_processing
        %selected_distal = candidates_data(:,11)>0;
        %selected_non_distal = candidates_data(:,12)>0;
        load([candidates_dir 'mat/' im_name '_caps.mat']);
    else   
        selected_distal = candidates_data(:,9) == 1;
        selected_non_distal = candidates_data(:,9) == 2;
    end
    
    if exist('metrics_dir', 'var')
        load([metrics_dir im_name '_am.mat'], 'apex_measures');
    end
    
    vessel_prob = u_load([prob_dir im_name '_pred.mat']);
    vessel_prob = conv2(g', g, vessel_prob, 'same');

    %Decide if these are hits or no
    [~, d_hits, d_cans] =...
        evaluate_apex_candidates(vessels.cluster_centres/2, candidate_xy(selected_distal,:),... 
        vessels.cluster_radius/2, vessel_prob, [], [], 0);    
    [~, nd_hits nd_cans] =...
        evaluate_apex_candidates(vessels.cluster_centres/2, candidate_xy(selected_non_distal,:),...
        vessels.cluster_radius/2, vessel_prob, [], [], 0);
    
    %Anything that wasn't a hit in d_cans or nd_cans can be added now, as
    %we know neither observer marked these
    co_occurrence(1,1,2,i_im) = sum(~d_cans);
    co_occurrence(1,1,3,i_im) = sum(~nd_cans);

    %Now loop through the vessels that were marked by observers
    markers = vessels.markers(1:2);
    for i_ve = 1:length(vessels.cluster_shapes)

        %Get column index based on observers 1
        idx = find(vessels.cluster_members{i_ve} == markers(1));
        if isempty(idx)
            c_i = 1;
        else
            shapes = vessels.cluster_shapes{i_ve};
            if strcmpi(shapes{idx}, 'NonDistal')
                c_i = 3;
            else
                c_i = 2;
                width1 = vessels.cluster_widths{i_ve}(idx)/2;
            end
        end

        %Get row index based on observer 2
        idx = find(vessels.cluster_members{i_ve} == markers(2));
        if isempty(idx)
            r_i = 1;
        else
            shapes = vessels.cluster_shapes{i_ve};
            if strcmpi(shapes{idx}, 'NonDistal')
                r_i = 3;
            else
                r_i = 2;
                width2 = vessels.cluster_widths{i_ve}(idx)/2;
            end
        end
        
        %Get 3rd dimension based on software (but only if 1 or 2 marked it)
        feature = args.width_feature;
        if r_i > 1 || c_i > 1
            if d_hits(i_ve)
                q_i = 2;
                if exist('metrics_dir', 'var')
                    if isfield(apex_measures, 'distal') && isfield(apex_measures.distal, feature)
                        width_c = apex_measures.distal.(feature)(d_hits(i_ve)) + args.width_fudge;
                    elseif isfield(apex_measures, feature)
                        width_c = apex_measures.(feature)(d_hits(i_ve)) + args.width_fudge;
                    else
                        width_c = 0;
                    end
                else
                    width_c = 0;
                end
            elseif nd_hits(i_ve)
                q_i = 3;
            else
                q_i = 1;
            end

            co_occurrence(r_i,c_i,q_i,i_im) = co_occurrence(r_i,c_i,q_i,i_im) + 1;
            
            if exist('metrics_dir', 'var')
                if r_i == 2 && c_i == 2 
                    if q_i == 2
                        a_b_widths{i_im} = [a_b_widths{i_im}; width1 width2];
                        a_c_widths{i_im} = [a_c_widths{i_im}; width1 width_c];
                        b_c_widths{i_im} = [b_c_widths{i_im}; width2 width_c];
                        detected_widths = [detected_widths; (width1 + width2)/2]; %#ok
                    else
                        missed_widths = [missed_widths; (width1 + width2)/2]; %#ok
                    end
                end
            end
            
        end
    end
end
%%
a_v_b = sum(sum(co_occurrence,4),3);
total_a = sum(a_v_b,1);
total_b = sum(a_v_b,2);
display('A vs B');
display(['% A distal missed by B = ' num2str(100*a_v_b(1,2) / total_a(2), 3) ', (' num2str(a_v_b(1,2)) ')']);
display(['% A distal marked as distal by B = ' num2str(100*a_v_b(2,2) / total_a(2), 3) ', (' num2str(a_v_b(2,2)) ')']);
display(['% A distal marked as non-distal by B = ' num2str(100*a_v_b(3,2) / total_a(2), 3) ', (' num2str(a_v_b(3,2)) ')']);
display('---');
display('B vs A');
display(['% B distal missed by A = ' num2str(100*a_v_b(2,1) / total_b(2), 3) ', (' num2str(a_v_b(2,1)) ')']);
display(['% B distal marked as distal by A = ' num2str(100*a_v_b(2,2) / total_b(2), 3) ', (' num2str(a_v_b(2,2)) ')']);
display(['% B distal marked as non-distal by A = ' num2str(100*a_v_b(2,3) / total_b(2), 3) ', (' num2str(a_v_b(2,3)) ')']);
display(['Accuracy = ' num2str(100*sum(diag(a_v_b)) / sum(a_v_b(:)),3)]);
display('-----------------------------------------------------------------');
%%
a_v_c = squeeze(sum(sum(co_occurrence,4),1))';
total_a = sum(a_v_c,1);
total_c = sum(a_v_c,2);
display('A vs C');
display(['% A distal missed by C = ' num2str(100*a_v_c(1,2) / total_a(2), 3) ', (' num2str(a_v_c(1,2)) ')']);
display(['% A distal marked as distal by C = ' num2str(100*a_v_c(2,2) / total_a(2), 3) ', (' num2str(a_v_c(2,2)) ')']);
display(['% A distal marked as non-distal by C = ' num2str(100*a_v_c(3,2) / total_a(2), 3) ', (' num2str(a_v_c(3,2)) ')']);
display('---');
display('C vs A');
display(['% C distal missed by A = ' num2str(100*a_v_c(2,1) / total_c(2), 3) ', (' num2str(a_v_c(2,1)) ')']);
display(['% C distal marked as distal by A = ' num2str(100*a_v_c(2,2) / total_c(2), 3) ', (' num2str(a_v_c(2,2)) ')']);
display(['% C distal marked as non-distal by A = ' num2str(100*a_v_c(2,3) / total_c(2), 3) ', (' num2str(a_v_c(2,3)) ')']);
display(['Accuracy = ' num2str(100*sum(diag(a_v_c)) / sum(a_v_c(:)),3)]);
display('-----------------------------------------------------------------');
%%
b_v_c = squeeze(sum(sum(co_occurrence,4),2))';
total_b = sum(b_v_c,1);
total_c = sum(b_v_c,2);
display('B vs C');
display(['% B distal missed by C = ' num2str(100*b_v_c(1,2) / total_b(2), 3) ', (' num2str(b_v_c(1,2)) ')']);
display(['% B distal marked as distal by C = ' num2str(100*b_v_c(2,2) / total_b(2), 3) ', (' num2str(b_v_c(2,2)) ')']);
display(['% B distal marked as non-distal by C = ' num2str(100*b_v_c(3,2) / total_b(2), 3) ', (' num2str(b_v_c(3,2)) ')']);
display('---');
display('C vs B');
display(['% C distal missed by B = ' num2str(100*b_v_c(2,1) / total_c(2), 3) ', (' num2str(b_v_c(2,1)) ')']);
display(['% C distal marked as distal by B = ' num2str(100*b_v_c(2,2) / total_c(2), 3) ', (' num2str(b_v_c(2,2)) ')']);
display(['% C distal marked as non-distal by B = ' num2str(100*b_v_c(2,3) / total_c(2), 3) ', (' num2str(b_v_c(2,3)) ')']);
display(['Accuracy = ' num2str(100*sum(diag(b_v_c)) / sum(b_v_c(:)),3)]);
display('-----------------------------------------------------------------');
%%
if exist('metrics_dir', 'var')
    
    feature_str = feature;
    feature_str(feature_str == '_') = ' ';
    feature_str(1) = feature_str(1) - 32;
        
    a_b_widths_all = zeros(0,2);
    a_c_widths_all = zeros(0,2);
    b_c_widths_all = zeros(0,2);
    for i_im = 1:num_images
        a_b_widths_all = [a_b_widths_all; a_b_widths{i_im}]; %#ok
        a_c_widths_all = [a_c_widths_all; a_c_widths{i_im}]; %#ok
        b_c_widths_all = [b_c_widths_all; b_c_widths{i_im}]; %#ok
    end

    figure; 
    subplot(1,2,1); hold all;
    plot(sort(abs(diff(a_b_widths_all,1,2))), linspace(0,1,size(a_b_widths_all,1)), 'linewidth', 2);
    plot(sort(abs(diff(a_c_widths_all,1,2))), linspace(0,1,size(a_c_widths_all,1)), 'linewidth', 2);
    plot(sort(abs(diff(b_c_widths_all,1,2))), linspace(0,1,size(b_c_widths_all,1)), 'linewidth', 2);
    title(['CDF of absolute differences in width prediction (' feature_str ')']);
    legend({'|B - A|', '|C - A|', '|C - B|'}, 'location', 'southeast');
    xlabel('Absolute width difference');
    ylabel('% of data');
    axis([0 20 0 1]);

    subplot(1,2,2); hold all;
    plot(sort(diff(a_b_widths_all,1,2)), linspace(0,1,size(a_b_widths_all,1)), 'linewidth', 2);
    plot(sort(diff(a_c_widths_all,1,2)), linspace(0,1,size(a_c_widths_all,1)), 'linewidth', 2);
    plot(sort(diff(b_c_widths_all,1,2)), linspace(0,1,size(b_c_widths_all,1)), 'linewidth', 2);
    title(['CDF of signed differences in width prediction (' feature_str ')']);
    legend({'B - A', 'C - A', 'C - B'}, 'location', 'southeast');
    xlabel('Width difference');
    ylabel('% of data');
    axis([-20 20 0 1]);
    plot([0 0], [0 1], 'k--');
    plot([-20 20], [0.5 0.5], 'k--');

    figure; hold all;
    plot(mean(a_b_widths_all,2), diff(a_b_widths_all,1,2), 'x');
    title(['Bland-altman plot for capillary widths, A v B (' feature_str ')']);
    xlabel('Mean width of A and B'); ylabel('B - A');
    plot([0 80], [0 0], 'k--');
    
    [~,m,b] = regression(mean(a_b_widths_all,2), diff(a_b_widths_all,1,2), 'one');
    x = [0 80];
    y = m*x + b;
    plot(x, y, 'y--');
    a1 = gca;
    
    figure; hold all;
    plot(mean(a_b_widths_all,2), diff(a_c_widths_all,1,2), 'x');
    title(['Bland-altman plot for capillary widths, A v C (' feature_str ')']);
    xlabel('Mean width of A and B'); ylabel('C - A');
    plot([0 80], [0 0], 'k--');
    [~,m,b] = regression(mean(a_c_widths_all,2), diff(a_c_widths_all,1,2), 'one');
    x = [0 80];
    y = m*x + b;
    plot(x, y, 'y--');
    a2 = gca;
    
    figure; hold all;
    plot(mean(a_b_widths_all,2), diff(b_c_widths_all,1,2), 'x');
    title(['Bland-altman plot for capillary widths, B v C (' feature_str ')']);
    xlabel('Mean width of A and B'); ylabel('C - B');
    plot([0 80], [0 0], 'k--');
    [~,m,b] = regression(mean(b_c_widths_all,2), diff(b_c_widths_all,1,2), 'one');
    x = [0 80];
    y = m*x + b;
    plot(x, y, 'y--');
    a3 = gca;
    
    linkaxes([a1 a2 a3]);
    axis([0 80 -40 40]);
end






















