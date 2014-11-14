function [vessels] = post_merge_apexes(vessels, vessel_prob, dist_thresh, patch_sz2, connect_thresh, n_connect_pts, reduction_factor, do_plot)
%POST_MERGE_APEXES Take an existing vessels structure (i.e. on computed
%from CLUSTER_VESSEL_APICES) and determine if any clusters should be merged
%on the basis of their connectivity in a detection probability map of vessels 
%   [vessels] = post_merge_apexes(vessels, vessel_prob, min_dist, patch_sz, connect_thresh, n_connect_pts)
%
% Inputs:
%      vessels - Existing vessel structure computed from CLUSTER_VESSEL_APICES
%
%      vessel_prob - vessel probability detection image
%
%      dist_thresh - maximum distance between pairs of clusters to be
%       tested for merging
%
%      patch_sz2 - Halfwidth of the patch extracted from vessel_prob when
%       computing the connectivity for a pair of clusters
%
%      connect_thresh - Threshold on the connectivity for a pair of
%       clusters to be merged
%
%      n_connect_pts - Number of points between 0 and 1 to test
%       connectivity
%
%
% Outputs:
%      vessels - New vessels structure, with any vessels merged as
%      necessary. note an addition field 'premerged_subgroups' is added,
%      which records the original composition of any clusters prior to them
%      being merged
%
%
% Example:
%
% Notes:
%
% See also: CLUSTER_VESSEL_APICES PLOT_APEX_CLUSTERS
%
% Created: 16-Oct-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('reduction_factor', 'var') || isempty(reduction_factor)
    reduction_factor = 1;
end
if ~exist('do_plot', 'var') || isempty(do_plot)
    do_plot = false;
end

%Add fields to store data on what gets merged and why
num_vessels = length(vessels.cluster_members);
vessels.premerged_subgroups = cell(num_vessels,1);

%Reduce down the patch_sz and dist_thresh if necessary (it may seem
%unintuitive to do it this way, but it easier to maintain standard defaults
%this way for images we've resized)
dist_thresh = dist_thresh / reduction_factor;
patch_sz2 = patch_sz2 / reduction_factor;

curr_vessel = 1;
plot_num = 1;
patch_sz = patch_sz2*2 + 1;
num_markers = length(vessels.markers);
while curr_vessel <= size(vessels.cluster_centres,1) && size(vessels.cluster_centres,1) > 1

    %Only look for vessels not marked by all markers
    if vessels.num_markers(curr_vessel) < num_markers &&...
        ~any(strcmpi(vessels.majority_shapes{curr_vessel}, {'NonDistal', 'Undefined'}));

        %Find the nearest vessel
        idx = setdiff(1:size(vessels.cluster_centres,1), curr_vessel);           
        dists = sum(bsxfun(@minus, vessels.cluster_centres(idx,:), vessels.cluster_centres(curr_vessel,:)).^2,2);               
        [min_dist min_idx] = min(dists);                
        nearest_vessel = idx(min_idx);

        %Check that the nearest vessel a) is within the minimum
        %distance b)is also not marked by everyone and c) is not
        %marked by a marker of c (does (c) -> (b), probably...)              
        if sqrt(min_dist) < dist_thresh && ...
            vessels.num_markers(nearest_vessel) < num_markers && ...
            isempty(intersect(vessels.cluster_members{curr_vessel}, vessels.cluster_members{nearest_vessel})) &&...
            ~any(strcmpi(vessels.majority_shapes{nearest_vessel}, {'NonDistal', 'Undefined'}));

            %Get the two centres and their midpoint
            centre1 = vessels.cluster_centres(curr_vessel,:) / reduction_factor;
            centre2 = vessels.cluster_centres(nearest_vessel,:) / reduction_factor;
            centre_midpoint = (centre1 + centre2) / 2;

            %Sample a patch about the midpoint
            vessel_prob_patch = sample_window(vessel_prob, patch_sz, round(centre_midpoint(2)), round(centre_midpoint(1)));
            %vessel_prob_patch = conv2(g', g, vessel_prob_patch, 'same');
            
            %Make the centres xy coords relative to the patch frame
            centre1 = centre1 - round(centre_midpoint) + patch_sz2;
            centre2 = centre2 - round(centre_midpoint) + patch_sz2;    

            %Work connectivity
            for i_con = fliplr(linspace(0, 1, n_connect_pts));
                
                con_mask = bwselect(vessel_prob_patch > i_con, centre1(1), centre1(2));
                if con_mask(round(centre2(2)), round(centre2(1)))                            
                    break;
                end
            end
            connectedness = i_con;
            
            if connectedness > connect_thresh
                %Make a cluster subgroup cell
                if isempty(vessels.premerged_subgroups{curr_vessel,:})
                    vessels.premerged_subgroups{curr_vessel,:} = vessels.cluster_xy(curr_vessel,:);
                end
                if isempty(vessels.premerged_subgroups{nearest_vessel,:})
                    vessels.premerged_subgroups{nearest_vessel,:} = vessels.cluster_xy(nearest_vessel,:);
                end
                subgroup = {vessels.premerged_subgroups{nearest_vessel,:}, vessels.premerged_subgroups{nearest_vessel,:}, connectedness};
                vessels.premerged_subgroups{curr_vessel,:} = subgroup;
                
                %Merge the vessels
                vessels.cluster_xy{curr_vessel,:} = [vessels.cluster_xy{curr_vessel,:}; vessels.cluster_xy{nearest_vessel,:}];
                vessels.cluster_members{curr_vessel,:} = [vessels.cluster_members{curr_vessel,:}; vessels.cluster_members{nearest_vessel,:}];
                vessels.cluster_shapes{curr_vessel,:} = [vessels.cluster_shapes{curr_vessel,:}; vessels.cluster_shapes{nearest_vessel,:}];
                vessels.cluster_sizes{curr_vessel,:} = [vessels.cluster_sizes{curr_vessel,:}; vessels.cluster_sizes{nearest_vessel,:}];
                vessels.cluster_widths{curr_vessel,:} = [vessels.cluster_widths{curr_vessel,:}; vessels.cluster_widths{nearest_vessel,:}];   

                %Delete nearest_vessel
                vessels.cluster_xy(nearest_vessel,:) = [];
                vessels.cluster_members(nearest_vessel,:) = [];
                vessels.cluster_shapes(nearest_vessel,:) = [];
                vessels.cluster_sizes(nearest_vessel,:) = [];
                vessels.cluster_widths(nearest_vessel,:) = [];
                vessels.cluster_centres(nearest_vessel,:) = [];
                vessels.cluster_radius(nearest_vessel,:) = [];
                vessels.premerged_subgroups(nearest_vessel,:) = [];

                %Recompute curr_vessel's mean
                valid_vessels = ~strcmpi(vessels.cluster_shapes{curr_vessel,:}, 'NonDistal') & ...
                    ~strcmpi(vessels.cluster_shapes{curr_vessel,:}, 'undefined');  
                if any(valid_vessels)
                    vessels.cluster_centres(curr_vessel,:) = mean(vessels.cluster_xy{curr_vessel,:}(valid_vessels,:),1);
                    vessels.cluster_radius(curr_vessel,:) = max(vessels.cluster_widths{curr_vessel,:}(valid_vessels,:));
                end
                
            end
            
            if do_plot
                %Display
                if plot_num == 1
                    figure;
                end
                subplot(3,4,plot_num); 
                imgray(vessel_prob_patch); 

                title([vessels.majority_shapes{curr_vessel} ', ' vessels.majority_shapes{nearest_vessel}]);
                xlabel(num2str(connectedness));

                plot(centre1(:,1), centre1(:,2), 'x', 'markersize', 20);
                plot(centre2(:,1), centre2(:,2), 'x', 'markersize', 20);

                if plot_num == 12
                    plot_num = 1;
                else
                    plot_num = plot_num + 1;
                end
            end
            
        end
    end
    %Increment the vessel count
    curr_vessel = curr_vessel + 1;
end

%Recompute the majority shapes
if size(vessels.cluster_centres,1) ~= length(vessels.num_markers)
    num_vessels = size(vessels.cluster_centres,1);
    vessels.majority_shapes = cell(num_vessels,1);
    vessels.majority_sizes = cell(num_vessels,1);
    vessels.num_markers = zeros(num_vessels,1);
    
    for i_v = 1:num_vessels;   

        [size_idx sizes] = grp2idx(vessels.cluster_sizes{i_v});
        [mode_idx] = mode(size_idx);
        vessels.majority_sizes{i_v} = sizes{mode_idx};

        [shape_idx shapes] = grp2idx(vessels.cluster_shapes{i_v});
        [mode_idx] = mode(shape_idx);
        vessels.majority_shapes{i_v} = shapes{mode_idx};

        vessels.num_markers(i_v) = length(vessels.cluster_sizes{i_v});

    end
end
