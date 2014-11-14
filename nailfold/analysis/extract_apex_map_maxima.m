function [] = extract_apex_map_maxima(varargin)
%EXTRACT_VESSEL_CENTRES *Insert a one line summary here*
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
args = u_packargs(varargin, 0, ... % the user's input
    {'image_names'},       ...
    'data_dir',             [nailfoldroot 'data/rsa_study/master_set/'],...
    'fov_mask_dir',         'fov_masks',...
    'centre_dir',           'vessel_centres',...
    'apex_map_dir',         'apex_maps',...
    'candidates_dir',       'local_maxima',...
    'exclusion_zone',       20,...
    'transform_scores',     0,...
    'separate_trees',  0,...
    'apex_class_thresh',    0.5,...
    'discard_pts',          0,...
    'base_width',           20,...
    'use_island_max',       1,...
    'island_thresh',        0,...
    'overwrite',            0);

%Form full directory paths and create folder for HoGs
fov_mask_dir = [args.data_dir '/' args.fov_mask_dir '/'];
centre_dir = [args.data_dir '/' args.centre_dir '/'];
apex_map_dir = [args.data_dir '/' args.apex_map_dir '/'];
candidates_dir = [args.data_dir '/' args.candidates_dir '/'];
create_folder(candidates_dir);

num_images = length(args.image_names);

%Loop though each image
for i_im = 1:num_images
    
    im_name = args.image_names{i_im} ;      
    if ~args.overwrite && exist([candidates_dir im_name '_candidates.mat'], 'file')
        display(['Image ' num2str(i_im) ', ' ...
            candidates_dir im_name '_candidates.mat'...
            ' already exists. Skipping (switch overwrite to 1 if necessary)']);
        continue;
    end
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);
    
    %Load in the image mask and the apex map
    load([apex_map_dir im_name '_pred.mat']);        
    f_mask = u_load([fov_mask_dir im_name '_f_mask.mat']);

    %See if we need to transform the map
    if args.transform_scores
        load([centre_dir im_name '_vc.mat']);
        
        %Discard edge points
        if args.discard_pts
            [discard_pts] = discard_edge_preds(vessel_centre, f_mask);
            include_pts = ~discard_pts ;
        else
            include_pts = true(size(apex_class_pred));
        end
        
        %Apply threshold and transform points
        include_pts = include_pts & (apex_class_pred > args.apex_class_thresh);
        [apex_offset_map] = ...
            transform_apex_offset_preds(apex_class_pred, apex_offset_x_pred, apex_offset_y_pred,...
                vessel_centre, nrows, ncols, args.base_width, include_pts, args.separate_trees);
        save([apex_map_dir im_name '_pred.mat'], ...
            'apex_offset_map', 'apex_class_pred', 'apex_offset_x_pred', 'apex_offset_y_pred');
    end

    if args.use_island_max
        [candidate_xy candidate_scores] = ...
           island_image_maxima(apex_offset_map, args.island_thresh, args.exclusion_zone, f_mask, 0, 0); %#ok
    else
        %Compute local maxima and save
        [candidate_xy candidate_scores] = ...
            local_image_maxima(apex_offset_map, args.exclusion_zone, f_mask, 0, 0); %#ok
    end
       
    save([candidates_dir im_name '_candidates.mat'], 'candidate_xy', 'candidate_scores');
end        
    
    
