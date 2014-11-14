function [] = display_automated_markup(varargin)
%DISPLAY_AUTOMATED_MARKUP *Insert a one line summary here*
%   [] = display_automated_markup(varargin)
%
% DISPLAY_AUTOMATED_MARKUP uses the U_PACKARGS interface function
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
% Created: 14-Apr-2014
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    0, ... % non-strict mode
    {'image_names'},         ...
    'data_dir',             [nailfoldroot 'data/2_year_study/'],...
    'image_dir',            'images',...
    'vessel_centre_dir',    'vessel_centres\',...
    'metrics_dir',          'apex_metrics',...
    'candidates_dir',       'apex_maps\local_maxima',...
    'selected_dir',         'selected_apexes',...
    'selected_features', [],...
    'do_distal',        1,...
    'do_nondistal',     1, ...
    'do_people_plots',  0,...
    'fig_dir',          [],...
    'um_per_pix',       1.25,...
    'xls_filename',     'auto_stats.xls',...
    'aam_thresh',       -2e4,...
    'plot_rejected',    1,...
    'plot_r', 3,...
    'plot_c', 1);
clear varargin;

image_dir = [args.data_dir '/' args.image_dir '/'];
%vessel_centre_dir = [args.data_dir '/' args.vessel_centre_dir '/'];
metrics_dir = [args.data_dir '/' args.metrics_dir '/'];
selected_dir = [args.data_dir '/' args.selected_dir '/'];
candidates_dir = [args.data_dir '/' args.candidates_dir '/'];

num_images = length(args.image_names);


do_capillary_type = [args.do_distal args.do_nondistal];
capillary_type = {'distal', 'nondistal'};

plots_per_fig = args.plot_r * args.plot_c;
plot_num = 1;

for i_im = 1:num_images
    im_name = args.image_names{i_im};
    
    if plot_num == 1
        figure;
    end
    if ~isempty(im_name)
        subplot(args.plot_r, args.plot_c, plot_num);
    end
    if plot_num == plots_per_fig
        plot_num = 1;
    else
        plot_num = plot_num + 1;
    end
    
    if isempty(im_name); continue; end
    
    nailfold = u_load([image_dir im_name '.mat']);
    apex_measures = u_load([metrics_dir im_name '_am.mat']);
    
    
    imgray(nailfold);
    title([im_name ': '...
        num2str(size(apex_measures.distal.median_width,1)) ' distal capillaries, ' ...
        num2str(size(apex_measures.nondistal.median_width,1)) ' non-distal capillaries']);
        
    for i_type = 1:2
    
        if ~do_capillary_type(i_type); continue; end
    
        num_cans = size(apex_measures.(capillary_type{i_type}).median_width,1);
        colors = lines(num_cans);

        max_width = 0; 
        max_width_id = 0;
        
        line_width = 3 - i_type;

        for i_can = 1:num_cans


            if isfield(apex_measures.(capillary_type{i_type}), 'apex_aam_fit') &&...
                    ~isempty(apex_measures.(capillary_type{i_type}).apex_aam_fit) && ...
                    isfield(apex_measures.(capillary_type{i_type}), 'apex_aam_score') &&...
                    apex_measures.(capillary_type{i_type}).apex_aam_score(i_can) > args.aam_thresh
                    
                apex_xy = apex_measures.(capillary_type{i_type}).apex_aam_fit(:,:,i_can);
                ax = apex_xy(16,1);
                ay = apex_xy(16,2);
                
                plot(apex_xy(:,1), apex_xy(:,2), '--', 'Color', 'r', 'linewidth', line_width); %colors(i_can,:)
                plot(ax, ay, 'x', 'Color', colors(i_can,:), 'MarkerSize', 4);       
                
            elseif isfield(apex_measures.(capillary_type{i_type}), 'apex_xy') && ~isempty(apex_measures.(capillary_type{i_type}).apex_xy);
                ax = apex_measures.(capillary_type{i_type}).apex_xy(i_can,1);
                ay = apex_measures.(capillary_type{i_type}).apex_xy(i_can,2);
                plot(ax, ay, 'x', 'Color', colors(i_can,:), 'MarkerSize', 4);
            else
                break;
            end

            a_theta = angle(apex_measures.(capillary_type{i_type}).base_orientation(i_can)) / 2;
            nx = sin(a_theta);
            ny = cos(a_theta);

            apex_width_x = ax + apex_measures.(capillary_type{i_type}).width_at_apex(i_can)*[-1 1]*nx/2;
            apex_width_y = ay + apex_measures.(capillary_type{i_type}).width_at_apex(i_can)*[-1 1]*ny/2;
            
            if apex_measures.(capillary_type{i_type}).mean_weighted_width(i_can) > max_width
                max_width = apex_measures.(capillary_type{i_type}).mean_weighted_width(i_can);
                max_width_id = i_can;
            end

            plot(apex_width_x, apex_width_y, '-', 'Color', colors(i_can,:), 'linewidth', line_width);
            
            if i_type == 1
                plot(apex_width_x, apex_width_y, 'o', 'MarkerEdgeColor', colors(i_can,:), 'MarkerSize', 8, 'MarkerFaceColor', colors(i_can,:));
            else
                plot(apex_width_x, apex_width_y, 'o', 'MarkerEdgeColor', colors(i_can,:), 'MarkerSize', 8);
            end
        end
        
        %Plot the maximum width vessel
        ax = apex_measures.(capillary_type{i_type}).apex_xy(max_width_id,1);
        ay = apex_measures.(capillary_type{i_type}).apex_xy(max_width_id,2);
        plot(ax, ay, 'r*', 'MarkerSize', 10);
    end
    if args.plot_rejected
        load([candidates_dir im_name '_candidates.mat'],...
            'candidate_xy');
        load([selected_dir im_name '_sel'],...
            'selected_distal', 'selected_non_distal', 'candidate_class_probs');
        rejected = ~selected_distal & ~selected_non_distal;
        rejected_xy = candidate_xy(rejected,:);
        plot(rejected_xy(:,1), rejected_xy(:,2), 'r.', 'markersize', 2);
            
%             rejected_probs = candidate_class_probs(rejected);
%             for i_r = 1:length(rejected_probs)
%                 plot(rejected_xy(i_r,1), rejected_xy(i_r,1), 'r.', 'markersize', 4*rejected_probs(i_r));
%             end
    end
end
