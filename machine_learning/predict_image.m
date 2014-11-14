function [image_out] = predict_image(varargin)
%PREDCIT_IMAGE predict output for an image using some prebuilt classifier/regressor
%   [probability_image] = predict_image(image_in,forest)
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
    'decomposition_args',...
    'predictor'}, ...
    'prediction_type', 'rf_regression',...
    'output_type', 'orientation',...
    'use_probs', 0,...
    'mask', [],...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 128,...
    'incremental_results', 0);
clear varargin;

%Get arguments needed to decompose image into the required feature vectors
decomposition_args = args.decomposition_args;

%Combine arguments to control the prediction
prediction_args.prediction_type = args.prediction_type;
prediction_args.use_probs = args.use_probs;
prediction_args.output_type = args.output_type;
prediction_args.num_trees = args.num_trees;

%Get size of image;
[ROW COL im_chs] = size(args.image_in);
    
%compute number of parts we need to break image into
r_parts = ceil(ROW / args.max_size);
c_parts = ceil(COL / args.max_size);

% Preallocate space for output image
out_dims = 1;
if strcmp(args.prediction_type, 'rf_classification') && length(args.predictor.classname) > 2;
    out_dims = length(args.predictor.classname);
end     
image_out = zeros(ROW*COL,out_dims);

%See if we have a tree mask and check it's the right size - we'll assume
%this has only been supplied for RF methods
if ~isempty(args.tree_mask)
    if ndims(args.tree_mask) == 3
        if any(size(args.tree_mask) ~= [ROW COL length(args.predictor.trees)])
            error('Tree mask is wrong size');
        else
            args.tree_mask = resize(args.tree_mask, ROW*COL, []);
        end
    else
        if any(size(args.tree_mask) ~= [ROW*COL length(args.predictor.trees)])
            error('Tree mask is wrong size');
        end
    end
end

% Pad image (we do this for all decomposition types now, so we don't need
% to worry about updating this function when we add new types)

if strcmpi(decomposition_args.decomp_type, 'dt')
    pad_w = 0;
else
    pad_w = floor(decomposition_args.win_size/2);
    args.image_in = padarray(args.image_in, [pad_w pad_w], 'replicate');
end

% Workout feature length from decomposition parameters
samples_per_channel = get_samples_per_channel(decomposition_args);
[n_channels, channels] = get_channels(decomposition_args);
sample_vector_size = samples_per_channel * n_channels;

% Compute responses to image filters
responses = cell(n_channels, 1);
for ch = 1:n_channels
    %select the desired channel
    if strcmpi(decomposition_args.rgb_channel, 'rgb')
        if (im_chs == 3)
            channel = rgb2gray(args.image_in);
        else
            % Image is monochrome
            channel = args.image_in;
        end
    else
        channel = args.image_in(:,:,channels(ch));
    end
    responses{ch} = compute_filter_responses(channel, decomposition_args);
end

%Go through each segment
rows = 1:min(args.max_size, ROW);
for rp = 1:r_parts

    cols = 1:min(args.max_size, COL);
    for cp = 1:c_parts
        % Get rows/cols subscripts for this part
        [part_cols, part_rows] = meshgrid(cols, rows);
        part_idx = sub2ind([ROW COL], part_rows, part_cols);
        
        % Move across one block for next iteration
        cols = cols + args.max_size;
        cols(cols > COL) = [];
        
        % Throw away pixels not belonging to a mask (if specified)
        if ~isempty(args.mask)
            part_rows(~args.mask(part_idx)) = [];
            part_cols(~args.mask(part_idx)) = [];
            part_idx(~args.mask(part_idx)) = [];
        end
        
        % Check whether there are any pixels left to process
        if isempty(part_rows), continue; end

        % Account for any padding added earlier
        part_rows = part_rows + pad_w;
        part_cols = part_cols + pad_w;
        
		% Now sample from the filter responses in each channel
        test_data = zeros(length(part_rows(:)), sample_vector_size);
        data_cols = 1:samples_per_channel;
        for ch = 1:n_channels
            test_data(:, data_cols) = sample_image_features(responses{ch}, ...
                                             part_rows(:), part_cols(:), ...
                                             decomposition_args);
            data_cols = data_cols + samples_per_channel;
        end
        
        %Extract the tree mask for this part (if a tree mask exists)
        if ~isempty(args.tree_mask)
            tree_mask = args.tree_mask(part_idx,:);
        else
            tree_mask = [];
        end
        
        % Apply the predictor to the test data
        image_out(part_idx,:) = apply_predictor(args.predictor, test_data, prediction_args, tree_mask);
        
    end
    
    % Move down one block
    rows = rows + args.max_size;
    rows(rows > ROW) = [];
end
image_out = reshape(image_out, ROW, COL, out_dims);

%% Helper functions
function [prob_image_part] = apply_predictor(predictor, test_data, prediction_args, tree_mask)

switch prediction_args.prediction_type
    case {'rf_classification','rf_regression'}
        if ~isempty(prediction_args.num_trees)
            if prediction_args.num_trees <= length(predictor.trees)
                predictor.trees(args.num_trees+1:end) = [];
            else
                error('The number of trees bigger than the total number of trees in the random forests! ');
            end
        end
end

%Now classify/predict the data using chosen RF method
switch prediction_args.prediction_type
    case 'rf_classification'
        [labels votes] = random_forest_class_predict(predictor,...
            test_data, tree_mask, prediction_args.use_probs); %#ok
        
        if length(predictor.classname) > 2
            prob_image_part = votes;
        else
            if isnumeric(predictor.classname) || ...
               islogical(predictor.classname)
                true_idx = find(predictor.classname);
            else
                true_idx = find(str2double(predictor.classname));
            end
            %Note that the rows of votes won't always sum to n_trees if
            %we're using a tree mask to compute OOB predictions
            prob_image_part = votes(:,true_idx) ./ sum(votes,2);
        end

    case 'rf_regression'
        switch prediction_args.output_type
            case {'orientation', 'centre_orientation'}
                is_orientation = true;
            otherwise
                is_orientation = false;
        end
        prob_image_part = random_forest_reg_predict(predictor, test_data, is_orientation, tree_mask);

        
    case 'linear_regression'
        prob_image_part = linear_regressor_predict(predictor, test_data);

    case 'boosted_regression'
        prob_image_part = boosted_regressor_predict(predictor, test_data);

    otherwise
        error(['Unknown regressor type: ',args.forest_type]);
end
