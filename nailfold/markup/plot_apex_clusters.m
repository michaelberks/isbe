function [fig_h] = plot_apex_clusters(vessels, fig_h, reduce_size)
%PLOT_APEX_CLUSTERS Plot the marked vessels for a nailfold image
%   [] = plot_apex_clusters(nailfold, nailfold_num, markers_dir, markers_list)
%
% Inputs:
%      vessels - structure computed from CLUSTER_VESSEL_APICES or POST_MERGE_APEXES
%
%      fig_h - if not set or not empty, will create new figure.
%       Alternatively the handle of figure in which to produce the plot
%
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also: CLUSTER_VESSEL_APICES POST_MERGE_APEXES
%
% Created: 06-Feb-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

num_vessels = size(vessels.cluster_centres,1);
num_markers = length(vessels.markers);   

if ~exist('fig_h', 'var') || isempty(fig_h)
    
    %Create figure
    fig_h = figure; hold all; axis equal ij;
else
    figure(fig_h);
end
if ~exist('reduce_size', 'var') || isempty(reduce_size)
    reduce_size = 1;
end

if num_markers > 1
    title_txt = ['Image marked by ' num2str(num_markers) ' observers, grades:'];
    for i_ma = 1:num_markers
        title_txt = [title_txt ' ' vessels.grades{i_ma}]; %#ok
    end
    
else
   title_txt = ['Image marked by 1 observers, grade: ' vessels.grades{1} ];    
end
 title(title_txt);   

%Plot initial markers for each vessel type and create a legend
plot(-1,-1, 'k^', 'visible', 'off'); % Normal
plot(-1,-1, 'ko', 'visible', 'off'); % Angiogenic
plot(-1,-1, 'ks', 'visible', 'off'); % Non-specific
plot(-1,-1, 'kx', 'visible', 'off'); % Undefined
plot(-1,-1, 'k+', 'visible', 'off'); % Non-distal

legend({'Normal shape', 'Angiogenic', 'Non-specific', 'Undefined', 'Non-distal'},...
    'location', 'southeast');

colors = lines(num_vessels);
pi_pts = linspace(0, 2*pi, 20);
for i_ve = 1:num_vessels;

    num_markers = size(vessels.cluster_xy{i_ve}(:,1),1);

    
    
    %Make a circle to show the width of the vessel marked
    outer_circ_x = vessels.cluster_centres(i_ve,1) + vessels.cluster_radius(i_ve)*cos(pi_pts)/2;
    outer_circ_y = vessels.cluster_centres(i_ve,2) + vessels.cluster_radius(i_ve)*sin(pi_pts)/2;
    plot(outer_circ_x/reduce_size, outer_circ_y/reduce_size,...
        'color', colors(i_ve,:), 'linewidth', num_markers);
    
    if num_markers > 1
        %Make a circle to enclose each marked apex in the cluster
        inner_circ_r2 = sum(bsxfun(@minus, vessels.cluster_xy{i_ve}, vessels.cluster_centres(i_ve,:)).^2,2);
        inner_circ_r = sqrt(max(inner_circ_r2));
        inner_circ_x = vessels.cluster_centres(i_ve,1) + inner_circ_r*cos(pi_pts);
        inner_circ_y = vessels.cluster_centres(i_ve,2) + inner_circ_r*sin(pi_pts);
        plot(inner_circ_x/reduce_size, inner_circ_y/reduce_size,...
            '-.', 'color', colors(i_ve,:));

        %Plot each marked point
        for i_ma = 1:num_markers
            switch vessels.cluster_shapes{i_ve}{i_ma}
                case 'Normal'
                    marker_type = '^';

                case 'Angiogenic'
                    marker_type = 'o';

                case 'Non-specific'
                    marker_type = 's';

                case 'Undefined'
                    marker_type = 'x';

                case 'NonDistal'
                    marker_type = '+';
            end

            plot(...
                vessels.cluster_xy{i_ve}(i_ma,1)/reduce_size,...
                vessels.cluster_xy{i_ve}(i_ma,2)/reduce_size,...
                marker_type, 'color', colors(i_ve,:));
        end
        
    end

    %Plot the centre of the cluster
    switch vessels.majority_shapes{i_ve}
        case 'Normal'
            marker_type = '^';

        case 'Angiogenic'
            marker_type = 'o';

        case 'Non-specific'
            marker_type = 's';

        case 'Undefined'
            marker_type = 'x';

        case 'NonDistal'
            marker_type = '+';
    end
    plot(vessels.cluster_centres(i_ve,1)/reduce_size,...
        vessels.cluster_centres(i_ve,2)/reduce_size,...
        marker_type, 'color', colors(i_ve,:), 'markersize', 10);
    
end

    