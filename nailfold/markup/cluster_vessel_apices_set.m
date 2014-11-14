function [] = cluster_vessel_apices_set(varargin)
%CLUSTER_VESSEL_APICES_SET *Insert a one line summary here*
%   [] = cluster_vessel_apices_set(varargin)
%
% CLUSTER_VESSEL_APICES_SET uses the U_PACKARGS interface function
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
% Created: 05-Dec-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, 0, ...
    {'selected_images',...
    'image_id_data'},...
    'test_dir',             'test',...
    'data_dir',             [nailfoldroot 'data/rsa_study/'],...
    'cluster_dir',          'apex_clusters',...
    'post_cluster_dir',     'apex_clusters_merged',...
    'apex_gt_dir',          'apex_gt',...
    'marker_list',          [],...
    'discard_repeats',      1,...
    'min_apex_dist',        20,... %The minimum distance allowed between separate clusters (the larger out of this value or the apex width is used for each cluster)
    'include_nondistal',    1,... %Record the position of non-distal vessels, and merge these with vessel marked as distal by a different observer if necessary
    'vessel_prob_dir',      'predictions/detection/rf_classification/257273/',...
    'dist_thresh',          100,... %What's the largest separation allowed between cluster centres
    'patch_sz2',            50,... %What's the halfwidth of the patch in which we look for connectivity
    'connect_thresh',       0.5,... %What's the threshold for deciding 2 clusters are connected?
    'n_connect_pts',        20,...%How many points (spaced between 0 and 1) to test the connectivity of pairs);
    'reduction_factor',     1, ... %Should we resize the ground truth (because we've resize images)
    'do_post_merge',       1 ...
); 
clear varargin;

%Set paths to directories containing the images, image markup, and where
%you want to save the clustered markup
cluster_dir = [args.data_dir args.test_dir '/' args.cluster_dir '/'];
post_cluster_dir = [args.data_dir args.test_dir '/' args.post_cluster_dir '/'];
apex_gt_dir = [args.data_dir args.test_dir '/' args.apex_gt_dir '/'];
vessel_prob_dir = [args.data_dir args.test_dir '/' args.vessel_prob_dir '/'];
create_folder(cluster_dir);
if args.do_post_merge
    create_folder(post_cluster_dir);
end
create_folder(apex_gt_dir);

for i_im = 1:length(args.image_id_data.im_names)
    
    if ~args.selected_images(i_im)
        continue;
    end
    
    display(['Processing image : ' num2str(i_im)]);
    
    im_name = args.image_id_data.im_names{i_im};
    marker_files = args.image_id_data.marker_files{i_im};    
    marker_idx = args.image_id_data.marker_idx{i_im};
    
    %Only use markers included on the marker lists
    if ~isempty(args.marker_list)
        markers = args.image_id_data.markers{i_im};
        keep = ismember(markers, args.marker_list);
        marker_idx = marker_idx(keep);
        marker_files = marker_files(keep);
    end
    
    %Don't use the repeat markings for a single observer - at moment this
    %must be true for cluster_vessel_apices to work properly...
    if args.discard_repeats
        [marker_idx keep] = unique(marker_idx, 'first');
        marker_files = marker_files(keep);
    end
    
    %Compute initial clustering of vessel markups and save it
    [vessels] = cluster_vessel_apices(marker_idx, marker_files, args.min_apex_dist, args.include_nondistal);
    if isempty(vessels.markers)
        %Noone has marked this image - continue
        continue;
    end
    %Otherwise save the intial vessel cluster and move on
    save([cluster_dir im_name '_apex_clusters.mat'], 'vessels');
    
    %For images marked by at least 2 observer, see if we need to post-merge
    %any vessels based on the vessel probability images
    if args.do_post_merge
        markers_per_image = length(vessels.markers);
        if markers_per_image > 1
            vessel_prob = u_load([vessel_prob_dir im_name '_pred.mat']);
            load([cluster_dir im_name '_apex_clusters.mat'], 'vessels');

            [vessels] = post_merge_apexes(vessels, vessel_prob, args.dist_thresh, args.patch_sz2,...
                args.connect_thresh, args.n_connect_pts, args.reduction_factor);   
        end
        save([post_cluster_dir im_name '_apex_clusters.mat'], 'vessels');
    end

    %Convert the vessels structure into a ground truth structure
    grades = vessels.grades;
    majority_grade = vessels.majority_grade; %#ok
    gradeable = ~any(...
        strcmpi(grades, 'Ungradeable_Quality') |...
        strcmpi(grades, 'Ungradeable_Condition') ); %#ok
    
    apex_xy = vessels.cluster_centres / args.reduction_factor; %#ok
    apex_widths = vessels.cluster_radius / args.reduction_factor; %#ok
    num_im_markers = length(vessels.markers); %#ok
    num_apex_markers = vessels.num_markers; %#ok
    is_non_distal = strcmpi(vessels.majority_shapes, 'NonDistal');
    is_undefined = strcmpi(vessels.majority_shapes, 'Undefined');
    is_distal = ~is_non_distal & ~is_undefined; %#ok

    apex_shape = vessels.majority_shapes; %#ok
    apex_size = vessels.majority_sizes; %#ok
       
    %Save a new apex GT structure
    save([apex_gt_dir im_name '_gt.mat'], 'gradeable', 'apex_xy', 'apex_widths', 'num_im_markers', 'num_apex_markers',...
        'is_non_distal', 'is_undefined', 'is_distal', 'grades', 'apex_size', 'apex_shape', 'majority_grade');
end
