function [vessels] = cluster_vessel_apices(markers, marker_files, min_apex_dist, include_nondistal, make_plot)
%CLUSTER_VESSEL_APICES Group observer markups of the same image to
%determine individual vessels
%   [] = cluster_vessel_apices(nailfold, nailfold_num, markers_dir, markers_list)
%
% Inputs:
%      markers - List of marker ID's
%
%      markers_files - Cell array of full filepaths to the markup files for
%      each observer
%
%      min_apex_dist - The minimum distance allowed between clusters (the
%       larger of this and the marked apex width is used for each individual
%       vessel cluster
%
%      include_nondistal - Set to 1 to include vessels marked as non-distal
%       by and observer in the clustering process. If included, non-distal
%       vessels can be added to an existing distal cluster, but will not
%       contribute to the cluster's mean position or mean width. Clusters that 
%       only contain non-distal vessels have their width set to the min_apex_dist  
%
%
% Outputs:
%      vessels: structure containing details of the clustered vessels
%
% Example:
%
% Notes:
%
% See also: POST_MERGE_APEXES PLOT_APEX_CLUSTERS
%
% Created: 06-Feb-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

if ~exist('min_apex_dist', 'var') || isempty(min_apex_dist)
    min_apex_dist = 15;
end

if ~exist('include_nondistal', 'var') || isempty(include_nondistal)
    include_nondistal = false;
end
if ~exist('make_plot', 'var') || isempty(make_plot)
    make_plot = false;
end

vessels.markers = markers;
vessels.grades = [];
vessels.cluster_centres = zeros(0,2);
vessels.cluster_xy = cell(0,1);
vessels.cluster_members = cell(0,1);
vessels.cluster_shapes = cell(0,1);
vessels.cluster_sizes = cell(0,1);
vessels.cluster_widths = cell(0,1);
vessels.cluster_radius = zeros(0,1);
vessels.med_dist = 0;

%Loop through markers_list
for i_ma = 1:length(marker_files)
    
    %Get the marker id number
    m_idx = markers(i_ma);
    
    %Read in vessel markup data - where multiple markup files exist,
    %use the last one
    markup = read_markup_from(marker_files{i_ma});
    vessels.grades{end+1,1} = markup.image_grade;
    
    num_vessels = length(markup.vessels);

    for i_v = 1:num_vessels

        %Check this is a valid vessel
        anchor_xy = markup.vessels(i_v).anchor;

        if isempty(anchor_xy); continue; end                           

        %Check if distal
        is_distal = markup.vessels(i_v).ncm_vessel_properties.is_distal;

        if is_distal
            %Find nearest marked apex to the anchor point
            num_apices = length(markup.vessels(i_v).apices);
            
            %A bit clumsy, but first loop through and tidy up any duplicate
            %and presumably unintended marks
            
            %First throaway any empty apices
            junk = false(num_apices,1);
            for i_a = 1:num_apices
                if isempty(markup.vessels(i_v).apices(i_a).inner_point)
                    junk(i_a) = 1;
                end
            end
            markup.vessels(i_v).apices(junk) = [];
            
            %Now throwaway duplicates
            num_apices = sum(~junk);
            junk = false(num_apices,1);
            for i_a1 = 1:num_apices
                apex1 = [...
                    markup.vessels(i_v).apices(i_a1).inner_point; ...
                    markup.vessels(i_v).apices(i_a1).outer_point ];
                
                apex_w1 = sum(diff(apex1).^2);
                for i_a2 = 1:i_a1-1
                    apex2 = [...
                        markup.vessels(i_v).apices(i_a2).inner_point; ...
                        markup.vessels(i_v).apices(i_a2).outer_point ];
                    apex_w2 = sum(diff(apex2).^2);
                    
                    do_cross = line_cross(apex1, apex2);
                    nearby = sum((mean(apex1)-mean(apex2)).^2) < max(apex_w1, apex_w2);
                
                    if do_cross || nearby
                        
                        if sum(apex1.^2) > sum(apex2.^2)
                            junk(i_a2) = 1;
                        else
                            junk(i_a1) = 1;
                        end
                    end
                end
            end
            markup.vessels(i_v).apices(junk) = [];           
            num_apices = sum(~junk);
            
            %Finally we can crack on with assigning the apexes to clusters
            for i_a = 1:num_apices
                
                apex_xy = (markup.vessels(i_v).apices(i_a).inner_point +...
                    markup.vessels(i_v).apices(i_a).outer_point) / 2;
                apex_width = sqrt(sum( (markup.vessels(i_v).apices(i_a).inner_point -...
                    markup.vessels(i_v).apices(i_a).outer_point).^2 ));
                
                vessel_shape = markup.vessels(i_v).ncm_vessel_properties.shape;
                vessel_size = markup.vessels(i_v).ncm_vessel_properties.size;
                
                if m_idx == 1
                    vessels = add_new_vessel(vessels, apex_xy, 1, vessel_size, vessel_shape, apex_width);
                else
                    [vessels vessel_num] = find_matching_vessel(vessels, apex_xy);
                    if vessel_num
                        vessels = add_apex_to_vessel(vessels, apex_xy, vessel_num, m_idx, vessel_size, vessel_shape, apex_width);
                    else
                        vessels = add_new_vessel(vessels, apex_xy, m_idx, vessel_size, vessel_shape, apex_width);
                    end
                end                
            end
        elseif include_nondistal
            %should we include non-distal vessels?
            apex_xy = anchor_xy;
            apex_width = min_apex_dist;
            
            vessel_shape = 'NonDistal';
            vessel_size = 'NonDistal';

            if m_idx == 1
                vessels = add_new_vessel(vessels, apex_xy, 1, vessel_size, vessel_shape, apex_width);
            else
                [vessels vessel_num] = find_matching_vessel(vessels, apex_xy);
                if vessel_num
                    vessels = add_apex_to_vessel(vessels, apex_xy, vessel_num, m_idx, vessel_size, vessel_shape, apex_width);
                else
                    vessels = add_new_vessel(vessels, apex_xy, m_idx, vessel_size, vessel_shape, apex_width);
                end
            end
        end
    end
    
    %Now do one final loop to merge any groups that started off separate bu
    %are similar enough to group
    i_ap1 = 1;
    num_apices = size(vessels.cluster_centres,1);
    while i_ap1 < num_apices
        apex1 = vessels.cluster_centres(i_ap1,:);
        members1 = vessels.cluster_members{i_ap1};
        
        %Loop through vessels next in the list from a1
        i_ap2 = i_ap1 + 1;
        while i_ap2 <= num_apices
            
            if isempty(intersect(members1,vessels.cluster_members{i_ap2}));
                min_apex_dist_i = max([...
                    vessels.cluster_radius(i_ap1),...
                    vessels.cluster_radius(i_ap2),...
                    min_apex_dist]);

                apex2 = vessels.cluster_centres(i_ap2,:);

                %if a1 and a2 not sufficiently separated...
                if sum((apex1-apex2).^2) < min_apex_dist_i^2;


                    %merge them, and update a1
                    vessels = merge_vessels(vessels, i_ap1, i_ap2);
                    apex1 = vessels.cluster_centres(i_ap1,:);
                    members1 = vessels.cluster_members{i_ap1};
                    
                    %decrement the total number of apexes - note we don't
                    %increment i_ap2 as the next vessels will have 'shifted
                    %down', so i_ap2 is already pointing at the next one
                    num_apices = num_apices - 1;
                else
                    %otherwise, increment i_ap2
                    i_ap2 = i_ap2 + 1;
                end
            else
                %otherwise, increment i_ap2
                i_ap2 = i_ap2 + 1;
            end
        end
        i_ap1 = i_ap1 + 1;
    end
    
    for i_ve = 1:size(vessels.cluster_centres,1)
        if length(unique(vessels.cluster_members{i_ve})) < length(vessels.cluster_members{i_ve})
            0;
        end
    end
end

if isempty(vessels.markers)
    %Noone has marked this image - return!
    return;
end

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

[grade_idx grades] = grp2idx(vessels.grades);
[mode_idx] = mode(grade_idx);
vessels.majority_grade = grades{mode_idx};

if make_plot
    plot_apex_clusters(vessels);     
end
    
%main function end
%%
%%-------------------------------------------------------------------------
%aux functions
function vessels = add_new_vessel(vessels, new_apex, marker, vessel_size, vessel_shape, apex_width)

    %Add the new vessel to the end of the list
    vessels.cluster_centres(end+1,:) = new_apex;
    vessels.cluster_xy{end+1,:} = new_apex;
    vessels.cluster_members{end+1,:} = marker;
    vessels.cluster_sizes{end+1,:} = {vessel_size};
    vessels.cluster_shapes{end+1,:} = {vessel_shape};
    vessels.cluster_widths{end+1,:} = apex_width;
    vessels.cluster_radius(end+1,:) = apex_width;
    
    %Sort the vessels from left to right in the image
    [vessels.cluster_centres sort_idx] = sortrows(vessels.cluster_centres);
    vessels.cluster_xy = vessels.cluster_xy(sort_idx);
    vessels.cluster_members = vessels.cluster_members(sort_idx);
    vessels.cluster_shapes = vessels.cluster_shapes(sort_idx);
    vessels.cluster_sizes = vessels.cluster_sizes(sort_idx);
    vessels.cluster_widths = vessels.cluster_widths(sort_idx);
    vessels.cluster_radius = vessels.cluster_radius(sort_idx);
    
    %Compute the mean distance between the vessels
    if size(vessels.cluster_centres, 1) > 1
        vessels.med_dist = median(sum(diff(vessels.cluster_centres).^2, 2));
    end
%%-------------------------------------------------------------------------

function vessels = add_apex_to_vessel(vessels, new_apex, vessel_num, marker, vessel_size, vessel_shape, apex_width)

    %Check if this marker has already been assigned to this vessel
    if ismember(marker, vessels.cluster_members{vessel_num,:})
        
        %if it has, work out whether the new or existing apex is closest to
        %the cluster centre
        idx = find(vessels.cluster_members{vessel_num,:} == marker);
        cluster_xy = vessels.cluster_xy{vessel_num,:};
        old_apex = cluster_xy(idx,:);
        cluster_xy(idx,:) = [];
        mean_xy = mean(cluster_xy);
        d_old = sum((mean_xy - old_apex).^2);
        d_new = sum((mean_xy - new_apex).^2);
        old_shape = vessels.cluster_shapes{vessel_num,:}{idx};
        old_size = vessels.cluster_sizes{vessel_num,:}{idx};
        old_width = vessels.cluster_radius(vessel_num);
        
        if d_new < d_old
            %Swap in the new apex and recalculate the cluster centre
            vessels.cluster_xy{vessel_num,:}(idx,:) = new_apex;
            vessels.cluster_shapes{vessel_num,:}{idx} = vessel_shape;
            vessels.cluster_sizes{vessel_num,:}{idx} = vessel_size;
            vessels.cluster_widths{vessel_num,:}(idx) = apex_width;
            
            %Don't change cluster mean pos radius for non-distal vessels
            valid_vessels = ~strcmpi(vessels.cluster_shapes{vessel_num,:}, 'NonDistal') & ...
                ~strcmpi(vessels.cluster_shapes{vessel_num,:}, 'Undefined'); 
            if any(valid_vessels)
                vessels.cluster_centres(vessel_num,:) = mean(vessels.cluster_xy{vessel_num,:}(valid_vessels,:),1);
                vessels.cluster_radius(vessel_num,:) = max(vessels.cluster_widths{vessel_num,:}(valid_vessels,:));
            end
            
            %Create a new vessel with the old apex
            vessels = add_new_vessel(vessels, old_apex, marker, old_size, old_shape, old_width);
        else
            %create a new vessel with the new apex
            vessels = add_new_vessel(vessels, new_apex, marker, vessel_size, vessel_shape, apex_width);
        end
    else
        %Add the new apex to the vessel and recalculate the cluster centre       
        vessels.cluster_xy{vessel_num,:}(end+1,:) = new_apex;
        vessels.cluster_members{vessel_num,:}(end+1,:) = marker;
        vessels.cluster_shapes{vessel_num,:}{end+1,:} = vessel_shape;
        vessels.cluster_sizes{vessel_num,:}{end+1,:} = vessel_size;
        vessels.cluster_widths{vessel_num,:}(end+1,1) = apex_width;
        
        %Don't change cluster mean pos radius for non-distal vessels
        valid_vessels = ~strcmpi(vessels.cluster_shapes{vessel_num,:}, 'NonDistal') & ...
            ~strcmpi(vessels.cluster_shapes{vessel_num,:}, 'Undefined');
        
        if any(valid_vessels)
            vessels.cluster_centres(vessel_num,:) = mean(vessels.cluster_xy{vessel_num,:}(valid_vessels,:),1);
            vessels.cluster_radius(vessel_num,:) = max(vessels.cluster_widths{vessel_num,:}(valid_vessels,:));
        end
        
    end
%%-------------------------------------------------------------------------

function [vessels vessel_num] = find_matching_vessel(vessels, new_apex)
   dists = sum(bsxfun(@minus,  vessels.cluster_centres, new_apex).^2, 2);
   [min_dist, vessel_num] = min(dists);
   if min_dist > vessels.cluster_radius(vessel_num,:) %min_apex_dist %(vessels.med_dist / 4)
       vessel_num = 0;
   end
%%-------------------------------------------------------------------------

function [vessels] = merge_vessels(vessels, v1, v2)

    %Add v2 into v1
    vessels.cluster_xy{v1,:} = [vessels.cluster_xy{v1,:}; vessels.cluster_xy{v2,:}];
    vessels.cluster_members{v1,:} = [vessels.cluster_members{v1,:}; vessels.cluster_members{v2,:}];
    vessels.cluster_shapes{v1,:} = [vessels.cluster_shapes{v1,:}; vessels.cluster_shapes{v2,:}];
    vessels.cluster_sizes{v1,:} = [vessels.cluster_sizes{v1,:}; vessels.cluster_sizes{v2,:}];
    vessels.cluster_widths{v1,:} = [vessels.cluster_widths{v1,:}; vessels.cluster_widths{v2,:}];   
    
    %Delete v2
    vessels.cluster_xy(v2,:) = [];
    vessels.cluster_members(v2,:) = [];
    vessels.cluster_shapes(v2,:) = [];
    vessels.cluster_sizes(v2,:) = [];
    vessels.cluster_widths(v2,:) = [];
    vessels.cluster_centres(v2,:) = [];
    vessels.cluster_radius(v2,:) = [];
    
    %Recompute v1's mean
    valid_vessels = ~strcmpi(vessels.cluster_shapes{v1,:}, 'NonDistal') & ...
        ~strcmpi(vessels.cluster_shapes{v1,:}, 'undefined');  
    if any(valid_vessels)
        vessels.cluster_centres(v1,:) = mean(vessels.cluster_xy{v1,:}(valid_vessels,:),1);
        vessels.cluster_radius(v1,:) = max(vessels.cluster_widths{v1,:}(valid_vessels,:));
    end
    
    
%%-------------------------------------------------------------------------
    




