function [] = extract_apex_candidate_class(varargin)
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
    'selected_ims',         [],...
    'data_dir',             [nailfoldroot 'data/rsa_study/set12g/'],...
    'prob_dir',             'rf_classification/296655/',...
    'candidates_dir',       'apex_maps/local_maxima',...
    'apex_gt_dir',          'apex_gt',...
    'label_dir',            'labels',...
    'prob_sigma',           2,...
    'overwrite',            0);

%Form full directory paths
prob_dir = [args.data_dir 'predictions/detection/' args.prob_dir '/'];
candidates_dir = [args.data_dir '/' args.candidates_dir '/'];
apex_gt_dir = [args.data_dir '/' args.apex_gt_dir '/'];
label_dir = [args.data_dir '/' args.label_dir '/'];
create_folder(label_dir);

num_images = length(args.image_names); 

if args.prob_sigma
    g_prob = gaussian_filters_1d(args.prob_sigma);
    g_prob = g_prob / sum(g_prob);
end

%Loop though each image
for i_im = 1:num_images
    
    im_name = args.image_names{i_im} ;  
    
    if ~args.overwrite && exist([label_dir im_name '_label.mat'], 'file')
        display(['Skipping image ' num2str(i_im) ',' label_dir im_name '_label.mat already exists. Re-run with overwrite=1 if necessary']); 
        continue;
    end
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);  
    
    %Load in images and vessel markup
    vessel_prob = u_load([prob_dir im_name '_pred.mat']);
    load([candidates_dir im_name '_candidates'],...
            'candidate_xy');
        
    if isempty(candidate_xy)
        continue;
    end
    load([apex_gt_dir im_name '_gt.mat'], 'apex_xy', 'apex_widths', 'is_distal', 'is_non_distal');
    
    if args.prob_sigma
        vessel_prob = conv2(g_prob', g_prob, vessel_prob, 'same');
    end

    num_candidates = size(candidate_xy,1);
            
    [~, ~, candidate_detections] =...
        evaluate_apex_candidates(apex_xy, candidate_xy, apex_widths, vessel_prob, ...
        [], [], 0);
    clear vessel_prob;
    
    candidates_class = candidate_detections > 0;
    candidate_labels = ones(num_candidates,1);
    apex_idx = candidate_detections(candidates_class);
    candidate_labels(candidates_class) = candidate_labels(candidates_class)...
        + is_distal(apex_idx) + 2*is_non_distal(apex_idx); %#ok
    
    save([label_dir im_name '_label.mat'], 'candidates_class', 'candidate_labels');
end          
    
    
