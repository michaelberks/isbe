function [] = detect_capillaries(nailfold_image, steps, varargin)
%DETECT_CAPILLARIES *Insert a one line summary here*
%   [] = detect_capillaries(varargin)
%
% DETECT_CAPILLARIES uses the U_PACKARGS interface function
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
% Created: 14-Mar-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0',...
    'nailfold_mask',            [],...
    'model_root',               'C:\isbe\nailfold\models\',...
    'detection_forest_path',    'vessel\detection\rf_classification\296655\predictor.mat',...
    'orientation_forest_path',  'vessel\orientation\rf_regression\296621\predictor.mat',...
    'width_forest_path',        'vessel\width\rf_regression\297037\predictor.mat',...
    'vessel_decomp_args',       'vessel\detection\rf_classification\296655\job_args.mat',...
    'hog_class_forest_path',    'apex\classification\set12g_half_296655\rf.mat',...
    'hog_off_x_forest_path',    'apex\offset_x\set12g_half_296655\rf.mat',...
    'hog_off_y_forest_path',    'apex\offset_y\set12g_half_296655\rf.mat',...
    'hog_rescore_forest_path',  'apex\rescoring\corrected_miccai_all\rf.mat',...
    'class_map_path',           'apex\final_MAP\miccai_class_MAP\class_map.mat',...
    'border_width',             32,...
    'max_size',                 1000,...
    'separate_trees',           0,...
    'image_sigma',              0,...
    'detection_sigma',          1,...
    'ori_sigma',                0,...
    'width_sigma',              1,...
    'exclusion_zone',           20,...
    'transform_scores',         0,...
    'apex_class_thresh',        0.5,...
    'discard_pts',              0,...
    'base_width',               20,...
    'num_cells',                8,...
    'cell_sz',                  8,... %Size of HoG cells in blocks
    'block_sz',                 [2 2],...%Size of blocks in cells
    'num_ori_bins',             9,... %Number of bins in orientation histograms
    'num_ori_bins_rescore',     12,... %Number of bins in orientation histograms
    'norm_method',              'none',... %Method for local normalisation 'none'
    'norm_method_rescore',      'l1-sqrt',... %Method for local normalisation 'none'
    'block_spacing',            8,... %Separation between blocks in pixels - controls level of overlap
    'gradient_operator',        [-1 0 1],...
    'spatial_sigma',            0, ...
    'angle_wrap',               1,...
    'min_candidates',           3,...
    'initial_thresh',           0.3,...
    'grid_spacing',             8,...
    'use_snake',                0,...
    'alpha',                    1,...
    'beta',                     1,...
    'snake_delta_y',            30,...
    'snake_res_y',              1,...
    'snake_spacing',            50,...
    'strong_vessel_thresh',     0.8,...
    'angle_discard_thresh',     75*pi/180,...
    'do_final_cull',            1,...
    'do_post_merge',            1,...
    'merge_dist_thresh',        120,...
    'merge_connect_thresh',     0.25,...
    'merge_n_pts',              20,...
    'do_aam',                   0,...
    'aam_model_path',           'aam/orig/2/vessel_apex_orig.smd',...
    'apex_width_forest_path',   [],...[nailfoldroot 'models/apex/width/rf.mat'],...
    'do_distal',                1,...
    'do_nondistal',             1,...
    'apex_patch_xs',            -31.5,...
    'apex_patch_xe',            31.5,...
    'apex_patch_ys',            -31.5,...
    'apex_patch_ye',            95.5,...
    'save_path',                [],...
    'plot',                     0);

clear varargin;

if ~isa(nailfold_image, 'double')
    nailfold_image = double(nailfold_image);
end

%Make a mask for the nailfold_image if one hasn't been supplied
if ~isempty(args.nailfold_mask)
    nailfold_mask = args.nailfold_mask;
    args = rmfield(args, 'nailfold_mask');
else
    nailfold_mask = make_nailfold_mosaic_mask(nailfold_image);
end

%Flag if we're saving outputs
do_save = ~isempty(args.save_path);

%Set-up which steps we'll do (default is all)
if ~exist('steps', 'var') || isempty(steps)
    steps = 1:7;
end
first_step = min(steps);
%
%% ------------------------------------------------------------------------
if ismember(1, steps)
    display('Processing step 1: Pixel level predictions of vessel presence, orientation and width');
    
    %Pixel level predictions of vessel presence, orientation and width
    rfs = cell(3,1);
    rfs{1} = u_load([args.model_root args.detection_forest_path]);
    rfs{2} = u_load([args.model_root args.orientation_forest_path]);
    rfs{3} = u_load([args.model_root args.width_forest_path]);
    rf_args = u_load([args.model_root args.vessel_decomp_args]);

    vessel_predictions = predict_image(... % non-strict mode
        'image_in', nailfold_image,...
        'decomposition_args', rf_args.decomposition_args,...
        'predictor', rfs, ...
        'prediction_type', {'rf_classification', 'rf_regression', 'rf_regression'},...
        'output_type', {'detection', 'orientation', 'width'},...
        'use_probs', 0,...
        'mask', nailfold_mask,...
        'tree_mask', [], ...
        'num_trees', [], ...
        'max_size', 1024,...
        'incremental_results', 0);
    
    if args.plot
        figure; imgray(vessel_predictions(:,:,1));
        figure; imgray(complex2rgb(vessel_predictions(:,:,2)));
    end
    if do_save
        save(args.save_path, 'vessel_predictions');
    end
end

%% ------------------------------------------------------------------------
if ismember(2, steps)
    display('Processing step 2: Extracting of vessel centres');
    
    %Extraction of vessel centres
    if first_step == 2
        load(args.save_path);
    end
    
    if args.image_sigma
        g = gaussian_filters_1d(args.image_sigma);
        g = g / sum(g);
        nailfold_image = conv2(g', g, nailfold_image, 'same');
    end
    if args.detection_sigma
        g = gaussian_filters_1d(args.detection_sigma);
        g = g / sum(g);
        vessel_predictions(:,:,1) = conv2(g', g, vessel_predictions(:,:,1), 'same');
    end
    if args.ori_sigma
        g = gaussian_filters_1d(args.ori_sigma);
        g = g / sum(g);
        vessel_predictions(:,:,2) = conv2(g', g, vessel_predictions(:,:,2), 'same');
    end 
    if args.width_sigma
        g = gaussian_filters_1d(args.width_sigma);
        g = g / sum(g);
        vessel_predictions(:,:,3) = conv2(g', g, vessel_predictions(:,:,3), 'same');
    end

    warning('off', 'ASYM:unexpectedArgument');

    [vessel_centre] =...
        extract_vessel_centres(vessel_predictions(:,:,1), vessel_predictions(:,:,2), vessel_predictions(:,:,3));

    [~, vessel_centre] = discard_edge_preds(vessel_centre, nailfold_mask, args.border_width);
    if do_save
        save(args.save_path,...
            'vessel_centre', '-append');    
    end
end
%% ------------------------------------------------------------------------
if ismember(3, steps)
    display('Processing step 3: Prediction of apex locations');
    
    %Prediction of apex locations
    if first_step == 3
        load(args.save_path);
    end
    
    apex_class_rf = u_load([args.model_root args.hog_class_forest_path]);
    apex_offset_x_rf = u_load([args.model_root args.hog_off_x_forest_path]);
    apex_offset_y_rf = u_load([args.model_root args.hog_off_y_forest_path]);

    % Get HoG args from main args
    hog_args.cell_sz = [args.cell_sz args.cell_sz];
    hog_args.block_sz = args.block_sz;
    hog_args.num_ori_bins = args.num_ori_bins;
    hog_args.norm_method = args.norm_method;
    hog_args.block_spacing = args.block_spacing;
    hog_args.gradient_operator = args.gradient_operator;
    hog_args.spatial_sigma = args.spatial_sigma;
    hog_args.angle_wrap = args.angle_wrap;

    % Get patch size and form template x,y coordinates for the patch
    patch_sz = args.num_cells*args.cell_sz;
    patch_sz = patch_sz + 2; %Account for padding
    patch_sz2 = (patch_sz - 1)/2;

    % Get hog size from the hog_args struct
    % Set up x,y coordinates for template patch
    x = repmat(-patch_sz2:patch_sz2, patch_sz, 1);
    y = x';
    xy = [x(:) y(:)];

    [apex_offset_map apex_class_pred apex_offset_x_pred apex_offset_y_pred] = ...
        predict_apex_offsets(...
            'apex_class_rf', apex_class_rf,...
            'apex_offset_x_rf', apex_offset_x_rf,...
            'apex_offset_y_rf', apex_offset_y_rf,...
            'vessel_feature_im', vessel_predictions(:,:,1), ...
            'vessel_centre', vessel_centre, ...
            'separate_trees', args.separate_trees,...
            'smoothing_sigma', 0,...
            'num_cells', args.num_cells,...
            'hog_args', hog_args,...
            'xy', xy,...
            'apex_class_thresh', args.apex_class_thresh,...
            'max_size', args.max_size,...
            'base_width', args.base_width,...
            'include_pts', []); %#ok

    %Compute local maxima and save
    [candidate_xy candidate_scores] = ...
        local_image_maxima(apex_offset_map, args.exclusion_zone, nailfold_mask, 0, 0); %#ok

    %Save this map
    if do_save
        save(args.save_path, 'apex_offset_map',...
            'apex_class_pred', 'apex_offset_x_pred', 'apex_offset_y_pred',...
            'candidate_xy', 'candidate_scores', '-append');

    end 
end
%% ------------------------------------------------------------------------
if ismember(4, steps)
    display('Processing step 4: Rescoring of apex candidates');
    
    %Rescoring of apex candidates
    if first_step == 4
        load(args.save_path);
    end
    
    %Only change we make to HoG args is to stopping wrapping angles
    hog_args.cell_sz = [args.cell_sz args.cell_sz];
    hog_args.block_sz = args.block_sz;
    hog_args.block_spacing = args.block_spacing;
    hog_args.gradient_operator = args.gradient_operator;
    hog_args.spatial_sigma = args.spatial_sigma;
    hog_args.angle_wrap = args.angle_wrap;
    hog_args.num_ori_bins = args.num_ori_bins_rescore;
    hog_args.norm_method = args.norm_method_rescore;
    hog_sz = args.num_cells*args.num_cells*hog_args.num_ori_bins;
    
    % Get patch size and form template x,y coordinates for the patch
    patch_sz = args.num_cells*args.cell_sz;
    patch_sz = patch_sz + 2; %Account for padding
    patch_sz2 = (patch_sz - 1)/2;

    % Get hog size from the hog_args struct
    % Set up x,y coordinates for template patch
    x = repmat(-patch_sz2:patch_sz2, patch_sz, 1);
    y = x';
    xy = [x(:) y(:)];

    if isempty(candidate_xy)
        candidates_hogs = [];
        candidate_oris = []; %#ok
        candidate_widths = [];
    else

        candidate_oris = interp2(vessel_predictions(:,:,2), candidate_xy(:,1), candidate_xy(:,2));
        candidate_widths = interp2(vessel_predictions(:,:,3), candidate_xy(:,1), candidate_xy(:,2));
        candidate_oris = angle(candidate_oris) / 2; 

        num_candidates = size(candidate_xy,1);
        candidates_hogs = zeros(num_candidates, hog_sz);

        figure; imgray(nailfold_image); hold all; a1 = gca;
        for i_can = 1:num_candidates

            %EXtract patch and compute HoG
            ori_c = candidate_oris(i_can);
            width_c = candidate_widths(i_can);

            %Get scale relative to base width a make rotation matrix
            rot = [cos(ori_c) -sin(ori_c); sin(ori_c) cos(ori_c)];
            scale = width_c / args.base_width;

            %Transform points given scale and angle and translate to
            %candidate position
            xya = xy * rot * scale;
            xa = reshape(xya(:,1) + candidate_xy(i_can,1), patch_sz, patch_sz);
            ya = reshape(xya(:,2) + candidate_xy(i_can,2), patch_sz, patch_sz);
            plot(a1, xa, ya, 'b.', 'markersize', 2);

            %Sample vessel prob patch
            vessel_feature_patch = interp2(nailfold_image, xa, ya, '*linear', 0);
            [hog] = compute_HoG(vessel_feature_patch, hog_args);       
            candidates_hogs(i_can,:) = hog(:)';
            
            if args.plot && i_can <= 20
                plot_hog(vessel_feature_patch, hog, args.num_cells, hog_args.angle_wrap);
            end
            
            %figure; imgray(vessel_feature_patch); colorbar;
        end
    end
    
    rescore_rf = u_load([args.model_root args.hog_rescore_forest_path]);
    [~,votes] = random_forest_class_predict(rescore_rf, candidates_hogs);
    candidate_rescores = votes(:,2) / length(rescore_rf.trees); clear votes;
    
    if do_save
        save(args.save_path,...
            'candidates_hogs', 'candidate_oris', 'candidate_widths',...
            'candidate_rescores', '-append');    
    end
end
    
%% ------------------------------------------------------------------------
if ismember(5, steps)
    display('Processing step 5: Computing spatial likelihood of distal row');
    
    %Computing spatial likelihood of distal row
    if first_step == 5
        load(args.save_path);
    end
    
    [nrows ncols] = size(nailfold_image);

    valid_candidates = candidate_rescores > args.initial_thresh;
    if sum(valid_candidates) < args.min_candidates
        candidate_displacements = inf(size(candidate_rescores));
    else
        %Create grid of points spaced over the image
        xx = 1:args.grid_spacing:ncols;
        yy = (1:args.grid_spacing:nrows)';
        grid_x = repmat(xx, length(yy), 1);
        grid_y = repmat(yy, 1, length(xx));

        %Compute weighted kernel estimates of the spatial distribution of
        %candidates over this grid
        [location_distribution] = build_2d_kernel_distribution(...
           candidate_xy(valid_candidates,:),...
           [grid_x(:) grid_y(:)],...
           candidate_rescores(valid_candidates,:), 1);

        location_distribution.D_f = reshape(location_distribution.D_f, size(grid_x));

        %Compute the max y-coordinate for each x-coordinate of the
        %candidates, and the displacement to this value
        [~, y_max] = max(location_distribution.D_f);
        x_i = xx;
        y_i = yy(y_max);

        if args.use_snake

            x_sampled = linspace(1, x_i(end), x_i(end) / (args.snake_spacing));
            y_sampled = interp1(x_i, y_i, x_sampled, 'linear');

            [snake_pnts] = mb_snake(...
                [x_sampled(:) y_sampled(:)] / args.grid_spacing,...
                args.alpha, args.beta, 0, 1, args.snake_delta_y, args.snake_res_y,...
                location_distribution.D_f * 1e6);

            x_i = snake_pnts(:,1)*args.grid_spacing;
            y_i = snake_pnts(:,2)*args.grid_spacing;     

        end

        candidate_polyfit = interp1(x_i, y_i, candidate_xy(:,1), 'linear');
        candidate_displacements = candidate_xy(:,2) - candidate_polyfit;

        if args.plot
            figure; imgray(location_distribution.D_f);
            plot(x_i/args.grid_spacing, y_i/args.grid_spacing, 'y');
            plot(candidate_xy(:,1)/args.grid_spacing, candidate_xy(:,2)/args.grid_spacing, 'rx');
            plot(...
                [candidate_xy(:,1) candidate_xy(:,1)]'/args.grid_spacing,...
                [candidate_xy(:,2) candidate_polyfit]'/args.grid_spacing, '-');
        end
    end    
    if do_save
        save(args.save_path,...
            'candidate_displacements', '-append');    
    end
end
%% ------------------------------------------------------------------------
if ismember(6, steps)
    display('Processing step 6: Selecting apexes from set of all candidates');
    
    %Selecting apexes from set of all candidates
    if first_step == 6
        load(args.save_path);
    end
    
    if all(candidate_displacements < inf)
        %Use class probs model to determine whether a candidate should be kept as a
        %distal or non-distal vessel
        status = 0; %#ok
        class_map = u_load([args.model_root args.class_map_path]);
        candidate_class = interp2(class_map.x, class_map.y, class_map.post_class,...
            candidate_rescores, candidate_displacements, 'linear', 0);
        candidate_class_probs = interp2(class_map.x, class_map.y, class_map.post_probs,...
            candidate_rescores, candidate_displacements, 'linear', 0);

        selected_distal = candidate_class == 1;
        selected_non_distal = candidate_class == 2;

        if args.do_post_merge

            [merged_with] = post_merge_candidates(candidate_xy, candidate_rescores, candidate_widths, ...
                 vessel_predictions(:,:,1), args.merge_dist_thresh, args.merge_connect_thresh, args.merge_n_pts, 0);

            for i_can = 1:length(merged_with)

                if merged_with(i_can)

                    if selected_distal(i_can)

                        if selected_distal(merged_with(i_can))
                            %Two distals merged together, discard the one with
                            %lower score
                            selected_distal(i_can) = 0;

                        elseif selected_non_distal(merged_with(i_can))
                            %Distal, merged with a non-distal. The
                            %non-distal has higher score, so must have been
                            %reject on its location - do nothing

                        else
                            %Distal, merged with something not selected as
                            %a non-distal, but has a higher vessel score
                            %than it. Swap the distal status
                            selected_distal(i_can) = 0;
                            selected_distal(merged_with(i_can)) = 1;
                        end
                    elseif selected_non_distal(i_can)

                        if selected_distal(merged_with(i_can))
                            %Non-distal, merged with a distal. Discard the
                            %non-distal altogether
                            selected_non_distal(i_can) = 0;

                        elseif selected_non_distal(merged_with(i_can))
                            %2 non-distals merged together, keep the one
                            %with the higher score
                            selected_non_distal(i_can) = 0;

                        else
                            %Non-istal, merged with something not selected as
                            %a non-distal, but has a higher vessel score
                            %than it. Swap the distal status
                            selected_non_distal(i_can) = 0;
                            selected_non_distal(merged_with(i_can)) = 1;
                        end
                    end
                end
            end

        end

        fell_at_the_last = false(size(candidate_rescores));
        %Finally, we discard any vessels that lie on top of one another
        if args.do_final_cull
            if any(selected_distal)
                %Which candidates belong to the distal row?
                [distal_idx] = select_distal_candidates(...
                    candidate_xy(selected_distal,:), candidate_rescores(selected_distal,:),...
                    10, args.angle_discard_thresh, 0, 0);

                %Work out which ones we've discarded (with indices relative to the
                %original list)
                distal_idx_tf = true(sum(selected_distal),1);
                distal_idx_tf(distal_idx) = 0;
                fell_at_the_last(selected_distal) = distal_idx_tf;

                %
                selected_distal(fell_at_the_last) = 0;
                selected_non_distal(fell_at_the_last) = 1;
            end
        end
    else
        %Minimum number of candidates not met to estimate displacements,
        %apply strong threshold
        selected_distal = candidate_rescores > args.strong_vessel_thresh;
        selected_non_distal = false(size(candidate_rescores));
        candidate_class = []; %#ok
        candidate_class_probs = [];
        fell_at_the_last = []; %#ok
        status = -1; %#ok
    end

    if do_save
        save(args.save_path,...
            'selected_distal', 'selected_non_distal', 'candidate_class_probs',...
            'candidate_class', 'status', '-append');    
    end
end

%% ------------------------------------------------------------------------
if ismember(7, steps)
    display('Processing step 7: Measuring vessel properties at selected apexes');
    %Measuring vessel properties at selected apexes
    if first_step == 7
        load(args.save_path);
    end
    
    if ~isempty(args.apex_width_forest_path)
        apex_width_rf = u_load(args.apex_width_forest_path);
    else
        apex_width_rf = [];
    end

    x = args.apex_patch_xs:args.apex_patch_xe;
    y = repmat((args.apex_patch_ys:args.apex_patch_ye)', 1, length(x));
    x = repmat(x, size(y,1), 1);

    xy = [x(:) y(:)];
    patch_sz = size(x);

    do_capillary_type = [args.do_distal args.do_nondistal];
    capillary_type = {'distal', 'nondistal'};

    apex_measures.distal = [];
    apex_measures.nondistal = [];

    for i_type = 1:2

        if ~do_capillary_type(i_type); continue; end

        if i_type == 1
            selected_idx = selected_distal;
        else
            selected_idx = selected_non_distal;
        end

        [measures_struc] = extract_apex_measures_new(...
            vessel_predictions(:,:,1), vessel_predictions(:,:,2), vessel_predictions(:,:,3),...
            nailfold_image, candidate_xy(selected_idx,:),... 
            'prob_sigma',           0,... %Already smoothed
            'ori_sigma',            0,... %Already smoothed
            'width_sigma',          0,... %Already smoothed
            'num_c_pts',            20,...
            'xy',                   xy,...
            'patch_sz',             patch_sz, ...
            'num_ori_bins',         36,...
            'connect_threshold',    0.5,...
            'width_predictor',      apex_width_rf,...
            'plot',                 args.plot && i_type==1);

        %Copying the fields in this way means we don't overwrite existing
        %data
        fnames = fieldnames(measures_struc);
        for i_f = 1:length(fnames)
            apex_measures.(capillary_type{i_type}).(fnames{i_f}) = measures_struc.(fnames{i_f});
        end

        apex_measures.(capillary_type{i_type}).candidate_scores = candidate_rescores(selected_idx,:);

        if ~isempty(candidate_class_probs)
            apex_measures.(capillary_type{i_type}).candidate_class_probs = candidate_class_probs(selected_idx,:);
        else
            apex_measures.(capillary_type{i_type}).candidate_class_probs = [];
        end
        if ~isempty(candidate_displacements)
            apex_measures.(capillary_type{i_type}).candidate_displacements = candidate_displacements(selected_idx,:);
        else
            apex_measures.(capillary_type{i_type}).candidate_displacements = [];
        end

        if args.do_aam
            load([args.model_dir 'mean_shape.mat'], 'mean_shape'); 
            num_aam_pts = size(mean_shape,1);
            aam_dir = [args.data_dir args.aam_dir '/'];

            initialise_aam_candidates([image_dir im_name '.mat'], ...
                vessel_predictions(:,:,3), vessel_predictions(:,:,2), candidate_xy(selected_idx,:), ...
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
                'aam_path', [args.model_dir args.aam_model_path],...
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

    if do_save
        save(args.save_path,...
            'apex_measures', '-append');    
    end
end  
%--------------------------------------------------------------------------
%End of main function
%--------------------------------------------------------------------------

