function [] = cluster_repeat_markings(varargin)
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
    'test_dir',             'test_half',...
    'data_dir',             [nailfoldroot 'data/rsa_study/'],...
    'cluster_dir',          'repeat_apex_clusters',...
    'marker_list',          [],...
    'discard_repeats',      1,...
    'min_apex_dist',        20,... %The minimum distance allowed between separate clusters (the larger out of this value or the apex width is used for each cluster)
    'include_nondistal',    1,... %Record the position of non-distal vessels, and merge these with vessel marked as distal by a different observer if necessary
    'vessel_prob_dir',      'predictions/detection/rf_classification/296655/',...
    'dist_thresh',          100,... %What's the largest separation allowed between cluster centres
    'patch_sz2',            50,... %What's the halfwidth of the patch in which we look for connectivity
    'connect_thresh',       0.5,... %What's the threshold for deciding 2 clusters are connected?
    'n_connect_pts',        20,...%How many points (spaced between 0 and 1) to test the connectivity of pairs);
    'reduction_factor',     2 ... %Should we resize the ground truth (because we've resize images)
); 
clear varargin;

%Set paths to directories containing the images, image markup, and where
%you want to save the clustered markup
cluster_dir = [args.data_dir args.test_dir '/' args.cluster_dir '/'];
vessel_prob_dir = [args.data_dir args.test_dir '/' args.vessel_prob_dir '/'];

%Make a new folder for each marker
for i_ma = 1:length(args.image_id_data.marker_list)
    create_folder([cluster_dir args.image_id_data.marker_list{i_ma}]);
end

num_images = length(args.image_id_data.im_names);
for i_im = 1:num_images
    
    if ~args.selected_images(i_im)
        continue;
    end
    
    m = args.image_id_data.marker_idx{i_im};
    [uni_m uni_i] = unique(m);
    n = length(m);
    
    %Check if this image has any repeat markings by the same observer
    if length(uni_m) < n
        
        %Get the image name 
        im_name = args.image_id_data.im_names{i_im};
        
        %Now workout which markers have repeated this image
        repeated_by = m(setdiff(1:n, uni_i));
        for i_ma = 1:length(repeated_by)
            
            %Cluster the two marker files for each marker
            keep = args.image_id_data.marker_idx{i_im} == repeated_by(i_ma);          
            if sum(keep) ~= 2
                error('Wrong number of matching marker files');
            end
            marker_files = args.image_id_data.marker_files{i_im}(keep);
            [vessels] = cluster_vessel_apices([1 2], marker_files, args.min_apex_dist, args.include_nondistal);
            
            %Post merge the vessels
            vessel_prob = u_load([vessel_prob_dir im_name '_pred.mat']);
            [vessels] = post_merge_apexes(vessels, vessel_prob, args.dist_thresh, args.patch_sz2,...
                args.connect_thresh, args.n_connect_pts, args.reduction_factor); %#ok
            
            %Save the output
            cluster_dir_i = [cluster_dir args.image_id_data.marker_list{repeated_by(i_ma)} '/'];
            save([cluster_dir_i im_name '_apex_clusters.mat'], 'vessels');
        end
    end
end
