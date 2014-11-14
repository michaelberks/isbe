function [probability_image] = classify_image_rgb(varargin)
%CLASSIFY_IMAGE *Insert a one line summary here*
%   [probability_image] = classify_image(image_in,forest)
%
% Inputs:
%      image_in - Image to classify
%
%      forest - Random forest classifier (model)
%
%
% Outputs:
%      probability_image - Pixel-wise probability of belonging to a bar
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 27-Jan-2010
% Author: Michael Berks
% Email : michael.berks@postgrad.man.ac.uk
% Phone : +44 (0)161 275 1241
% Copyright: (C) University of Manchester

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    {'image_in',... % the mandatory arguments
    'sampling_args',...
    'forest'}, ...
    'decomp_type', 'dt',...
    'forest_type', 'isbe',...
    'use_probs', 0,...
    'mask', [],...
    'num_trees', [], ...
    'max_size', 128);
clear varargin;
sampling_args = args.sampling_args;

%Get size of image;
ROW = size(args.image_in,1);
COL = size(args.image_in,2);

%compute number of parts we need to break image into
r_parts = ceil(ROW / args.max_size);
c_parts = ceil(COL / args.max_size);

probability_image = zeros(ROW, COL);

%Copy arguments used in sampling funtion into new structure
switch args.decomp_type
    case 'dt'
        if length(sampling_args.num_levels) == 1
            sampling_args.levels = 1:sampling_args.num_levels;
        else
            sampling_args.levels = sampling_args.num_levels;
        end
        %Compute dual-tree transform of image
        dt_r = compute_dual_tree(args.image_in(:,:,1), sampling_args.levels(end), sampling_args.use_nag);
        dt_g = compute_dual_tree(args.image_in(:,:,2), sampling_args.levels(end), sampling_args.use_nag);
        dt_b = compute_dual_tree(args.image_in(:,:,3), sampling_args.levels(end), sampling_args.use_nag);
        sampling_args = rmfield(sampling_args, {'num_levels', 'use_nag'});
        args = rmfield(args, 'image_in');
        
    case 'mono'
        %Compute monogenic transform of image
        pad_w = floor(sampling_args.win_size/2);
        args.image_in = padarray(args.image_in, [pad_w pad_w], 'replicate');
        [local_amp local_phase local_ori] = monogenic(args.image_in,...
            sampling_args.num_levels,...
            sampling_args.min_wavelength, 2,...
            sampling_args.onf, 1);
        sampling_args = rmfield(sampling_args, {'num_levels', 'min_wavelength', 'onf'});
        args = rmfield(args, 'image_in');
        
    case 'g2d'
        pad_w = floor(sampling_args.win_size/2);
        args.image_in = padarray(args.image_in, [pad_w pad_w], 'replicate');
        g2d_responses = compute_gaussian_2nd_derivatives(...
            args.image_in,  sampling_args.sigma_range);
        args = rmfield(args, 'image_in');
        
    case 'clover'
        pad_w = floor(sampling_args.win_size/2);
        args.image_in = padarray(args.image_in, [pad_w pad_w], 'replicate');
        g2d_responses = compute_clover_responses(...
            args.image_in,  sampling_args.sigma_range);
        args = rmfield(args, 'image_in');
    
    case 'g1d'
        pad_w = floor(sampling_args.win_size/2);
        args.image_in = padarray(args.image_in, [pad_w pad_w], 'replicate');
        g1d_responses = compute_gaussian_1st_derivatives(...
            args.image_in,  sampling_args.sigma_range);
        args = rmfield(args, 'image_in');
        
    case 'g2di'
        pad_w = floor(sampling_args.win_size/2);
        args.image_in = padarray(args.image_in, [pad_w pad_w], 'replicate');
        g2d_responses = compute_gaussian_2nd_derivatives_d(...
            args.image_in,...
            sampling_args.sigma_range(1),...
            sampling_args.sigma_range(2));
        sampling_args = rmfield(sampling_args, {'sigma_range'});
        args = rmfield(args, 'image_in');
        
    case 'linop'
        pad_w = floor(sampling_args.win_size/2);
        args.image_in = padarray(args.image_in, [pad_w pad_w], 'replicate');
    case 'pixel'
        pad_w = floor(sampling_args.win_size/2);
        args.image_in = padarray(args.image_in, [pad_w pad_w], 'replicate');
    case 'radial'
        %Do nothing just makes sure we don't trigger otherwise...
    otherwise
        warning(['decomposition type: ', args.decomp_type, ' not recognised, using DT-CWT']); %#ok
        %Compute dual-tree transform of image
        dt_r = compute_dual_tree(args.image_in(:,:,1), sampling_args.levels(end), sampling_args.use_nag);
        dt_g = compute_dual_tree(args.image_in(:,:,1), sampling_args.levels(end), sampling_args.use_nag);
        dt_b = compute_dual_tree(args.image_in(:,:,1), sampling_args.levels(end), sampling_args.use_nag);
        sampling_args = rmfield(sampling_args, {'num_levels', 'use_nag'});
        args = rmfield(args, 'image_in');
end

%Go through each segment
for rp = 1:r_parts
    for cp = 1:c_parts
        
        sr = 1 + (rp-1)*args.max_size;
        er = min(rp*args.max_size, ROW);
        sc = 1 + (cp-1)*args.max_size;
        ec = min(cp*args.max_size, COL);
        
        %Get rows/cols subscripts for this part
        [part_cols part_rows] = meshgrid(sc:ec,sr:er);
        part_idx = sub2ind([ROW COL], part_rows, part_cols);
        
        %Check whether we've been given a mask to select specific pixels
        if ~isempty(args.mask)
            
            %Throw away pixels not belonging to the mask
            part_rows(~args.mask(part_idx)) = [];
            part_cols(~args.mask(part_idx)) = [];
            part_idx(~args.mask(part_idx)) = [];
            
            %check whether there's any pixels left to process
            if isempty(part_rows)
                continue;
            end
        end
        
        %Now sample the decomposition data from the partition rows and columns
        %according to the sampling arguments
        switch args.decomp_type
            case 'dt' 
                test_data_r = sample_dt_data(dt_r, part_rows(:), part_cols(:), sampling_args);
                test_data_g = sample_dt_data(dt_g, part_rows(:), part_cols(:), sampling_args);
                %test_data_b = sample_dt_data(dt_b, part_rows(:), part_cols(:), sampling_args);

            case 'mono'
                test_data_r = sample_monogenic_data(...
                    local_amp, local_phase, local_ori, part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);
                test_data_g = sample_monogenic_data(...
                    local_amp, local_phase, local_ori, part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);
                test_data_b = sample_monogenic_data(...
                    local_amp, local_phase, local_ori, part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);
                
            case {'g2d', 'clover'}
                test_data_r = sample_g2d_data(...
                    g2d_responses(:,:,:,1),... 
                    g2d_responses(:,:,:,2),...
                    g2d_responses(:,:,:,3), part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);
                test_data_g = sample_g2d_data(...
                    g2d_responses(:,:,:,1),... 
                    g2d_responses(:,:,:,2),...
                    g2d_responses(:,:,:,3), part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);
                test_data_b = sample_g2d_data(...
                    g2d_responses(:,:,:,1),... 
                    g2d_responses(:,:,:,2),...
                    g2d_responses(:,:,:,3), part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);
                
            case 'g1d'
                test_data_r = sample_g2d_data(...
                    g1d_responses(:,:,:,1),... 
                    g1d_responses(:,:,:,2), part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);
                test_data_g = sample_g2d_data(...
                    g1d_responses(:,:,:,1),... 
                    g1d_responses(:,:,:,2), part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);
                test_data_b = sample_g2d_data(...
                    g1d_responses(:,:,:,1),... 
                    g1d_responses(:,:,:,2), part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);
            
            case 'g2di'
                test_data_r = sample_g2d_data_d(g2d_responses, part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);
                test_data_g = sample_g2d_data_d(g2d_responses, part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);
                test_data_b = sample_g2d_data_d(g2d_responses, part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);
                
            case 'linop'
                test_data_r = sample_linop_data(args.image_in, part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);
                test_data_g = sample_linop_data(args.image_in, part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);
                test_data_b = sample_linop_data(args.image_in, part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);

            case 'pixel'
                %Sample pixels from image
                test_data_r = sample_pixel_data(args.image_in, part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);
                test_data_g = sample_pixel_data(args.image_in, part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);
                test_data_b = sample_pixel_data(args.image_in, part_rows(:)+pad_w, part_cols(:)+pad_w, sampling_args);
            
            case 'radial'
                %Sample pixels from image
                test_data_r = sample_radial_data(part_idx(:), sampling_args);
                test_data_g = sample_radial_data(part_idx(:), sampling_args);
                test_data_b = sample_radial_data(part_idx(:), sampling_args);

            otherwise
                test_data_r = sample_dt_data(dt_r, part_rows(:), part_cols(:), sampling_args);
                test_data_g = sample_dt_data(dt_r, part_rows(:), part_cols(:), sampling_args);
                test_data_b = sample_dt_data(dt_r, part_rows(:), part_cols(:), sampling_args);
        end       
        
        %Now classify/predict the data using chosen RF method
        switch args.forest_type
            case 'breiman'
                [labels votes] = classRF_predict([test_data_r test_data_g test_data_b], args.forest);
                %prob_image_part = reshape(votes(:,1) / args.forest.ntree, size(part_cols));
                prob_image_part = votes(:,1) / args.forest.ntree;
                
            case 'isbe'
                if ~isempty(args.num_trees)
                    if args.num_trees<=length(args.forest.trees)
                        args.forest.trees(args.num_trees+1:end)=[];
                    else
                        error('The number of trees bigger than the total number of trees in the random forests! ');
                    end
                end
                if isnumeric(args.forest.classname) || islogical(args.forest.classname)
                    true_idx = find(args.forest.classname);
                else
                    true_idx = find(str2double(args.forest.classname));
                end
                
                if args.use_probs
                    %[labels votes] = mb_random_forest_prob_predict(args.forest, [test_data_r test_data_g test_data_b]);
                    [labels votes] = mb_random_forest_prob_predict(args.forest, [test_data_r test_data_g]);
                else
                    %[labels votes] = mb_random_forest_class_predict(args.forest, [test_data_r test_data_g test_data_b]);
                    [labels votes] = mb_random_forest_class_predict(args.forest, [test_data_r test_data_g]);
                end
                %prob_image_part = reshape(votes(:,true_idx) / length(args.forest.trees), size(part_cols)); %#ok
                prob_image_part = votes(:,true_idx) / length(args.forest.trees); %#ok
                
            case 'isbe_boot'
                if ~isempty(args.num_trees)
                    if args.num_trees<=length(args.forest.trees)
                        args.forest.trees(args.num_trees+1:end)=[];
                    else
                        error('The number of trees bigger than the total number of trees in the random forests! ');
                    end
                end
                if isnumeric(args.forest.classname) || islogical(args.forest.classname)
                    true_idx = find(args.forest.classname);
                else
                    true_idx = find(str2double(args.forest.classname));
                end
                [labels votes] = mb_random_forest_class_predict_boot(args.forest, [test_data_r test_data_g test_data_b]);
                %prob_image_part = reshape(votes(:,true_idx) / length(args.forest.trees), size(part_cols)); %#ok
                prob_image_part = votes(:,true_idx) / length(args.forest.trees); %#ok
                
                
            case 'regression'
                if ~isempty(args.num_trees)
                    if args.num_trees<=length(args.forest.trees)
                        args.forest.trees(args.num_trees+1:end)=[];
                    else
                        error('The number of trees bigger than the total number of trees in the random forests! ');
                    end
                end
                %[orientations] = mb_random_forest_reg_predict(args.forest, [test_data_r test_data_g test_data_b]);
                %prob_image_part = reshape(orientations, size(part_cols));
                prob_image_part = mb_random_forest_reg_predict(args.forest, [test_data_r test_data_g test_data_b]);
                
            case 'orientation'
                if ~isempty(args.num_trees)
                    if args.num_trees<=length(args.forest.trees)
                        args.forest.trees(args.num_trees+1:end)=[];
                    else
                        error('The number of trees bigger than the total number of trees in the random forests! ');
                    end
                end
                %[orientations] = mb_random_forest_reg_predict(args.forest, [test_data_r test_data_g test_data_b]);
                %prob_image_part = reshape(orientations, size(part_cols));
                prob_image_part = mb_random_forest_ori_predict(args.forest, [test_data_r test_data_g test_data_b]);
                
            case 'linear_regression'
                %Compute output - args.forest should contain sets of
                %regression coefficients for a number of precomputed linear
                %regressions. linear_regression_ouput then applies each
                %individual regressor and takes the mean of these
                %predictions as the output. Note I've not yet created
                %linear_regression_predict so feel free to use a different
                %function name!
                prob_image_part = linear_regressor_predict(args.forest, [test_data_r test_data_g test_data_b]);

			case 'boosted_regression'
                %Compute output using boosted regressor
                prob_image_part = boosted_regressor_predict(args.forest, [test_data_r test_data_g test_data_b]);
				
            otherwise
                error('Type of random forest not recognised');
        end
        
        %probability_image(sr:er, sc:ec) = prob_image_part;
        probability_image(part_idx) = prob_image_part;
        
    end
end
%%
% 'num_levels', 5,...
%     'feature_shape', 'rect',...
%     'feature_type', 'all',...
%     'do_max', 0,...
%     'rotate', 0,...
%     'win_size', 3,...
%     'num_angles', 8,...
%     'min_wavelength', 4,...
%     'onf', 0.65,...
%     'pca', [],...