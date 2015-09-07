function [] = display_individual_markup(nailfold, analysis, varargin)
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
    'selected_features', [],...
    'plot_distal',        1,...
    'plot_nondistal',     1, ...
    'um_per_pix',       1.25,...
    'aam_thresh',       -2e4,...
    'plot_rejected',    1,...
    'plot_view', [],...
    'plot_caxis', []);
clear varargin;




do_capillary_type = [args.plot_distal args.plot_nondistal];
capillary_type = {'distal', 'nondistal'};

figure;
imgray(nailfold);
title([...
    num2str(size(analysis.apex_measures.distal.median_width,1)) ' distal capillaries, ' ...
    num2str(size(analysis.apex_measures.nondistal.median_width,1)) ' non-distal capillaries']);

for i_type = 1:2

    if ~do_capillary_type(i_type); continue; end

    num_cans = size(analysis.apex_measures.(capillary_type{i_type}).median_width,1);
    colors = lines(num_cans);

    max_width = 0; 
    max_width_id = 0;

    line_width = 3 - i_type;

    for i_can = 1:num_cans


        if isfield(analysis.apex_measures.(capillary_type{i_type}), 'apex_aam_fit') &&...
                ~isempty(analysis.apex_measures.(capillary_type{i_type}).apex_aam_fit) && ...
                isfield(analysis.apex_measures.(capillary_type{i_type}), 'apex_aam_score') &&...
                analysis.apex_measures.(capillary_type{i_type}).apex_aam_score(i_can) > args.aam_thresh

            apex_xy = analysis.apex_measures.(capillary_type{i_type}).apex_aam_fit(:,:,i_can);
            ax = apex_xy(16,1);
            ay = apex_xy(16,2);

            plot(apex_xy(:,1), apex_xy(:,2), '--', 'Color', 'r', 'linewidth', line_width); %colors(i_can,:)
            plot(ax, ay, 'x', 'Color', colors(i_can,:), 'MarkerSize', 4);       

        elseif isfield(analysis.apex_measures.(capillary_type{i_type}), 'apex_xy') && ~isempty(analysis.apex_measures.(capillary_type{i_type}).apex_xy);
            ax = analysis.apex_measures.(capillary_type{i_type}).apex_xy(i_can,1);
            ay = analysis.apex_measures.(capillary_type{i_type}).apex_xy(i_can,2);
            plot(ax, ay, 'x', 'Color', colors(i_can,:), 'MarkerSize', 4);
        else
            break;
        end

        a_theta = angle(analysis.apex_measures.(capillary_type{i_type}).base_orientation(i_can)) / 2;
        nx = sin(a_theta);
        ny = cos(a_theta);

        apex_width_x = ax + analysis.apex_measures.(capillary_type{i_type}).width_at_apex(i_can)*[-1 1]*nx/2;
        apex_width_y = ay + analysis.apex_measures.(capillary_type{i_type}).width_at_apex(i_can)*[-1 1]*ny/2;

        if analysis.apex_measures.(capillary_type{i_type}).mean_weighted_width(i_can) > max_width
            max_width = analysis.apex_measures.(capillary_type{i_type}).mean_weighted_width(i_can);
            max_width_id = i_can;
        end

        plot(apex_width_x, apex_width_y, '-', 'Color', colors(i_can,:), 'linewidth', line_width);

        if i_type == 1
            plot(apex_width_x, apex_width_y, 'o', 'MarkerEdgeColor', colors(i_can,:), 'MarkerSize', 4, 'MarkerFaceColor', colors(i_can,:));
        else
            plot(apex_width_x, apex_width_y, 'o', 'MarkerEdgeColor', colors(i_can,:), 'MarkerSize', 4);
        end
    end

    %Plot the maximum width vessel
    if ~isempty(analysis.apex_measures.(capillary_type{i_type}).apex_xy)
        ax = analysis.apex_measures.(capillary_type{i_type}).apex_xy(max_width_id,1);
        ay = analysis.apex_measures.(capillary_type{i_type}).apex_xy(max_width_id,2);
        plot(ax, ay, 'r*', 'MarkerSize', 10);
    end
end
if args.plot_rejected
    rejected = ~analysis.selected_distal & ~analysis.selected_non_distal;
    rejected_xy = analysis.candidate_xy(rejected,:);
    plot(rejected_xy(:,1), rejected_xy(:,2), 'r.', 'markersize', 2);
            
%             rejected_probs = analysis.candidate_class_probs(rejected);
%             for i_r = 1:length(rejected_probs)
%                 plot(rejected_xy(i_r,1), rejected_xy(i_r,1), 'r.', 'markersize', 4*rejected_probs(i_r));
%             end
end
if ~isempty(args.plot_view)
    axis(args.plot_view);
end
if ~isempty(args.plot_caxis);
    caxis([args.plot_caxis]);
end

