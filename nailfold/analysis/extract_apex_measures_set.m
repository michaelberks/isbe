function [] = extract_apex_measures_set(varargin)
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
args = u_packargs(varargin,... % the user's input
    0, ... % non-strict mode
    {'image_names'},         ...
    'data_dir',             [nailfoldroot 'data/rsa_study/master_set/'],...
    'image_dir',            'images',...
    'prob_dir',             'rf_classification/296655',...
    'ori_dir',              'rf_regression/296621',...
    'width_dir',            'rf_regression/297037',...
    'candidates_dir',       'apex_maps\local_maxima',...
    'displacements_dir',    [],...
    'selected_dir',         [],...
    'metrics_dir',          'apex_metrics',...
    'aam_dir',              'aam',...
    'do_aam',               1,...
    'model_dir',            [nailfoldroot 'data/rsa_study/models/apex_templates/'],...
    'aam_name',             'aam/orig/2/vessel_apex_orig.smd',...
    'width_predictor',      [],...
    'prob_sigma',           2,...
    'ori_sigma',            0,...
    'width_sigma',          2,...
    'do_distal',            1,...
    'do_nondistal',         1,...
    'overwrite',            0,...
    'xs',                   -31.5,...
    'xe',                   31.5,...
    'ys',                   -31.5,...
    'ye',                   95.5,...
    'plot', 0);

if ~isempty(args.width_predictor)
    width_predictor = u_load(args.width_predictor);
else
    width_predictor = [];
end

image_dir = [args.data_dir '/' args.image_dir '/'];
prob_dir = [args.data_dir 'predictions/detection/' args.prob_dir '/'];
ori_dir = [args.data_dir 'predictions/orientation/' args.ori_dir '/'];
width_dir = [args.data_dir 'predictions/width/' args.width_dir '/'];
candidates_dir = [args.data_dir '/' args.candidates_dir '/'];
displacements_dir = [args.data_dir '/' args.displacements_dir '/'];
selected_dir = [args.data_dir '/' args.selected_dir '/'];
aam_dir = [args.data_dir args.aam_dir '/'];

if ~isempty(args.metrics_dir)
    metrics_dir = [args.data_dir '/' args.metrics_dir '/'];
    create_folder(metrics_dir);
end


num_images = length(args.image_names);

x = args.xs:args.xe;
y = repmat((args.ys:args.ye)', 1, length(x));
x = repmat(x, size(y,1), 1);

xy = [x(:) y(:)];
patch_sz = size(x);

if args.do_aam
    load([args.model_dir 'mean_shape.mat'], 'mean_shape'); 
    num_aam_pts = size(mean_shape,1);
end

do_capillary_type = [args.do_distal args.do_nondistal];
capillary_type = {'distal', 'nondistal'};
    
apex_measures.distal = [];
apex_measures.nondistal = [];

%Loop though each image
for i_im = 1:num_images
    im_name = args.image_names{i_im};
    
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);  
    
    %Load in images and vessel markup
    vessel_im = u_load([image_dir im_name '.mat']);
    vessel_prob = u_load([prob_dir im_name '_pred.mat']);
    vessel_ori = u_load([ori_dir im_name '_pred.mat']);
    vessel_width = u_load([width_dir im_name '_pred.mat']);
    load([candidates_dir im_name '_candidates.mat'],...
        'candidate_xy', 'candidate_rescores');
    load([displacements_dir im_name '_dis.mat'],...
        'candidate_displacements');
    load([selected_dir im_name '_sel'],...
        'selected_distal', 'selected_non_distal', 'candidate_class_probs');
    
    if ~args.overwrite && ~isempty(args.metrics_dir) && exist([metrics_dir im_name '_am.mat'], 'file')
        load([metrics_dir im_name '_am.mat'], 'apex_measures');
    end
    
    for i_type = 1:2
    
        if ~do_capillary_type(i_type); continue; end
        
        if i_type == 1
            selected_idx = selected_distal;
        else
            selected_idx = selected_non_distal;
        end

    
        [measures_struc] = extract_apex_measures(...
            vessel_prob, vessel_ori, vessel_width, vessel_im, candidate_xy(selected_idx,:),... 
            'prob_sigma',           args.prob_sigma,...
            'ori_sigma',            args.ori_sigma,...
            'width_sigma',          args.width_sigma,...
            'num_c_pts', 20,...
            'xy', xy,...
            'patch_sz', patch_sz, ...
            'num_ori_bins', 36,...
            'connect_threshold', 0.5,...
            'width_predictor', width_predictor,...
            'plot', args.plot); %#ok
        
        %Copying the fields in this way means we don't overwrite existing
        %data
        fnames = fieldnames(measures_struc);
        for i_f = 1:length(fnames)
            apex_measures.(capillary_type{i_type}).(fnames{i_f}) = measures_struc.(fnames{i_f});
        end
        
        apex_measures.(capillary_type{i_type}).candidate_scores = candidate_rescores(selected_idx,:);  %#ok
        
        if ~isempty(candidate_class_probs) %#ok
            apex_measures.(capillary_type{i_type}).candidate_class_probs = candidate_class_probs(selected_idx,:);
        else
            apex_measures.(capillary_type{i_type}).candidate_class_probs = [];
        end
        if ~isempty(candidate_displacements) %#ok
            apex_measures.(capillary_type{i_type}).candidate_displacements = candidate_displacements(selected_idx,:);
        else
            apex_measures.(capillary_type{i_type}).candidate_displacements = [];
        end

        if args.do_aam
            initialise_aam_candidates([image_dir im_name '.mat'], vessel_width, vessel_ori, candidate_xy(selected_idx,:), ...
                    'aam_dir', [aam_dir im_name '/'],...
                    'mean_shape', mean_shape,...
                    'base_width', 15,...
                    'width_sigma', 2,...
                    'ori_sigma', 0,...
                    'max_num_candidates', inf,...
                    'thetas', linspace(-pi/8, pi/8, 20),...
                    'scales', [0.8 0.9 1.0 1.1 1.2],...
                    'debug', 0);

            fit_aam_to_candidates(...
                'aam_dir', [aam_dir im_name '/'],...
                'aam_exe', 'ncm_sandpit_mb',...
                'aam_path', [args.model_dir args.aam_name],...
                'delete_candidate_patches', 1);

            num_cans = sum(selected_idx);
            apex_measures.(capillary_type{i_type}).apex_aam_fit = zeros(num_aam_pts,2,num_cans);
            apex_measures.(capillary_type{i_type}).apex_aam_score = zeros(num_cans,1);
            for i_can = 1:num_cans
                load([aam_dir im_name '\apex' zerostr(i_can, 4) '_aam.mat'], 'apex_candidate');
                apex_measures.(capillary_type{i_type}).apex_aam_fit(:,1,i_can) = apex_candidate.fitted_vessel_xy(:,1) + apex_candidate.sc;
                apex_measures.(capillary_type{i_type}).apex_aam_fit(:,2,i_can) = apex_candidate.fitted_vessel_xy(:,2) + apex_candidate.sr;
                apex_measures.(capillary_type{i_type}).apex_aam_score(i_can) = apex_candidate.model_score;
                delete([aam_dir im_name '\apex' zerostr(i_can, 4) '_aam.mat']);
            end
            %rmdir([aam_dir im_name '/'], 's');
        end
    end
    if ~isempty(args.metrics_dir)
        save([metrics_dir im_name '_am.mat'], 'apex_measures');    
    end
end    

    
    
