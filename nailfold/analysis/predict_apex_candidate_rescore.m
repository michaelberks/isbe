function [] = predict_apex_candidate_rescore(varargin)
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
    {'image_names',...
    'apex_class_rf'},       ...
    'data_dir',             [nailfoldroot 'data/rsa_study/master_set/'],...
    'candidates_dir',       'apex_maps/local_maxima',...
    'rescore_dir',          'apex_maps/local_maxima',...
    'hog_dir',              'apex_hogs',...
    'reselect_maxima',      0,...
    'island_thresh',        0,...
    'apex_map_dir',         []);

%Form full directory paths and create folder for HoGs
rescore_dir = [args.data_dir '/' args.rescore_dir '/'];
hog_dir = [args.data_dir '/' args.hog_dir '/'];
candidates_dir = [args.data_dir '/' args.candidates_dir '/'];
create_folder(rescore_dir);

if args.reselect_maxima
    apex_map_dir = [args.data_dir '/' args.apex_map_dir '/'];
end
num_images = length(args.image_names);

%Loop though each image
for i_im = 1:num_images
    
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);
    
    im_name = args.image_names{i_im} ;  
    
    load([hog_dir im_name '_hog.mat'], 'candidates_hogs');
    load([candidates_dir im_name '_candidates'], 'candidate_xy');
    [~,votes] = random_forest_class_predict(args.apex_class_rf, candidates_hogs);
    candidate_rescores = votes(:,2) / length(args.apex_class_rf.trees); clear votes;
    
    if args.reselect_maxima
        load([apex_map_dir im_name '_pred.mat'], 'apex_offset_map');
        [candidate_xy, candidate_rescores, discarded_xy, discarded_rescores] = ...
            select_island_maxima(candidate_xy, candidate_rescores, apex_offset_map>args.island_thresh); %#ok
        
    else
        discarded_xy = []; %#ok
        discarded_rescores = []; %#ok
    end
    
    save([rescore_dir im_name '_candidates.mat'], ...
        'candidate_xy', 'candidate_rescores', 'discarded_xy', 'discarded_rescores');
end          
    
    
