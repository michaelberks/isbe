function [training_data training_labels] = generate_training_data(varargin)
%SAMPLE_VESSEL_DT_DATA *Insert a one line summary here*
%   [] = sample_mammo_training_data(varargin)
%
% SAMPLE_VESSEL_DT_DATA uses the U_PACKARGS interface function
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
% Created: 17-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

f_debug = (nargin == 0 && nargout == 0);
if f_debug, test_script(); return; end

[training_data training_labels] = func(varargin{:});


%% The function
function [training_data training_labels] = func(varargin)
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    ... % Mandatory arguments
   {'sampling_args', ...
    'decomposition_args'},...
    'quiet', true);
clear varargin;

% Copy arguments used in sampling function into new structure
decomposition_args = args.decomposition_args;
sampling_args = args.sampling_args;

% Workout feature length from decomposition parameters
samples_per_channel = get_samples_per_channel(decomposition_args);
[n_channels, channels] = get_channels(decomposition_args);
sample_vector_size = samples_per_channel * n_channels;

% Preferred code (not yet implemented):
% samples_per_pixel = get_samples_per_pixel(decomposition_args);
% n_pixels = sampling_args.win_size^2;
% [n_channels, channels] = get_channels(sampling_args);
% sample_vector_size = samples_per_pixel * n_pixels * n_channels;

% Create memory storage for data
n_samples = sampling_args.num_samples; 
training_data = zeros(n_samples, sample_vector_size);
training_labels = zeros(n_samples, 1);

% Get list of all the images, fov masks and vessel masks - we're going to
% assume these all match up for now
switch sampling_args.image_type
    case {'real'}
        image_lists = create_image_lists(sampling_args);
        n_images = length(image_lists);
        sampled_pts = cell(n_images,1);
        
    case {'line', 'grain'} 
        % Synthetic image types
        
        % Estimated number of images we'll need.
        % (The way thay the number of points sampled per image is later chosen 
        % ensures that we have all the points that we need.)
        image_lists = [];
        n_images = sampling_args.num_samples / sampling_args.pts_per_image;
        
    otherwise
        error(['Image type ', sampling_args.image_type, ' not recognized']);
end

% Allow user to specify a maximum number of images from which to sample
max_n_images_is_valid = isfield(sampling_args, 'max_n_images') && ...
                        ~isempty(sampling_args.max_n_images) && ...
                        isnumeric(sampling_args.max_n_images) && ...
                        sampling_args.max_n_images > 0;
if max_n_images_is_valid
    n_images = min(n_images, sampling_args.max_n_images);
end

% Go back to beginning of the image lists (if they are used).
get_next_image('reset');

% Sample data from every image
n_sampled_so_far = 0;


if ~args.quiet; tb = timebar('limit', n_samples, 'title', 'Sampling data'); end
for img_index = 1:n_images
    % Work out the number of samples to take
    N = n_samples - n_sampled_so_far;
    p = 1 / (n_images - img_index + 1);
    samples_per_image = sample_from_binomial(N, p, 1);
    
    if (samples_per_image == 0)
        continue;
    end
    
    % NB: We could put the resampling code here
    
    % Retrieve the input image, its ground truth labelling and 
    % foreground/background maps.
    [img, img_label, fg_map, bg_map, synth_params] = ...
        get_next_image(sampling_args, image_lists);

    % If parameters are nonzero (i.e. the parameters of a synthetic image)
    % then store them for posterity.
    if ~isempty(synth_params)
        if ~exist('parameters','var')
            % First one - make n_images copies that will be overwritten.
            parameters(1:n_images) = synth_params;
        else
            parameters(img_index) = synth_params;
        end
    end

    % Choose sample locations in foreground and background - this function
    % now checks we have enough samples and returns the best it can
    [fg_indices, bg_indices] = ...
        sample_from_maps(fg_map, bg_map, samples_per_image, sampling_args);
    
    samples_per_image = numel(fg_indices) + numel(bg_indices);
    display(['Samples in image ' num2str(img_index) ': ', ...
                                 num2str(samples_per_image)]);                             
    
    if exist('sampled_pts', 'var')
        %save which points we've sampled
        sampled_pts{img_index} = [fg_indices; bg_indices];
    end

    %Check again if we found anything in the maps, if so, skedadle
    if (samples_per_image == 0)
        continue;
    end
    
    % Copy values from ground truth label images
    switch sampling_args.output_type
        case {'detection', 'centre_detection', ...
              'junction_detection', 'junction_centre_detection',...
              'orientation', 'centre_orientation', 'width'}
            % Sample labels from label image
            labels = [img_label(fg_indices); img_label(bg_indices)];
            
        case {'mixed_orientation', 'mixed_centre_orientation'}
            % TO DO: Needs aligning with get_*_image_label()
            
            % We could deal with this in get_real_image_label() so that the
            % ground truth is a single plane image where the orientation
            % at each pixel has already been sampled from the options.
            % Then we could ditch the switch statement altogether.
            
            %now randomly select one of the mixed orientations for each
            %point
            	
            %Load in orientation ground truth
%             ori_gt = load([args.ori_dir ori_list(this_im).name]);
            labels = zeros(num_samples_image,1);
            for ii = 1:samples_per_image
                mixed_oris = ...
                    img_label.mixed_oris{img_label.mixed_idx(image_idx(ii))};
                
                labels(ii) = mixed_oris(ceil(rand*length(mixed_oris)));
            end
            
        otherwise
            error(['Output type ', sampling_args.output_type, ' not recognized']);
    end   
    
    % Compute range of rows in training data
    first_sample = n_sampled_so_far + 1;
    last_sample  = n_sampled_so_far + samples_per_image;

    % Save the sample labels for this image
    training_labels(first_sample:last_sample) = labels;
    
    % Convert indices to row/col subscripts
    [rows cols] = ind2sub([size(img,1) size(img,2)], ...
                          [fg_indices; bg_indices]);

    data_cols = 1:samples_per_channel;
    for ch = 1:n_channels
        %select the desired channel
        if strcmpi(decomposition_args.rgb_channel, 'rgb')
            if (size(img,3)==3)
                channel = rgb2gray(img);
            else
                % Image is monochrome
                channel = img;
            end
        else
            channel = img(:,:,channels(ch));
        end
        
        % Get training feature vectors and store with rest of training data
        responses = compute_filter_responses(channel, decomposition_args);
        training_data(first_sample:last_sample, data_cols) = ...
            sample_image_features(responses, rows, cols, decomposition_args);

        data_cols = data_cols + samples_per_channel;
    end

    % Update the current sample count
    n_sampled_so_far = n_sampled_so_far + samples_per_image;
    
	% update timebar
	if ~args.quiet; timebar(tb, 'advance', samples_per_image); end
end
if ~args.quiet; timebar(tb,'close'); clear tb; end

% All training features/labels are now in. Save the parameter lists if they
% exist - a bit messy to have this here, but to keep the tree training code
% general the generate training data function can only return X & y. So we
% use the model wrapper function to pass a path that the model will be
% saved to, so the parameters/sampling points data can also be saved to
% this path
if ~isempty(sampling_args.sampled_data_dir)
    create_folder(sampling_args.sampled_data_dir);
    par_list = dir([sampling_args.sampled_data_dir '/sampled_pts*.mat']);
    tree_num = length(par_list) + 1;
    if exist('sampled_pts','var') && ~isempty(sampled_pts)
        save([sampling_args.sampled_data_dir '/sampled_pts' zerostr(tree_num, 3) '.mat'], 'sampled_pts');
    end
    if exist('parameters','var') && ~isempty(parameters)
        save([sampling_args.sampled_data_dir 'parameters' zerostr(tree_num, 3) '.mat'], ...
             'parameters');
    end
    if ~exist('image_lists','file') && ~isempty(image_lists)
        save([sampling_args.sampled_data_dir '/image_lists.mat'], 'image_lists');
    end
end


%% Test script
function test_script()
clc;
args = default_args();
args.image_root = [asymmetryroot,'data'];
args.output_type = 'detection';
args.decomp_type = 'dt';
args.num_samples = 100;
args.sampling_method = 'generate_training_data';

% Test with synthetic data
args.image_type = 'line';
args.pts_per_image = 20;
args.bg_dir = prettypath([args.image_root, '/', args.bg_dir]);

newargs = [];
newargs.sampling_args = get_sampling_args_from(args);
newargs.sampling_args = get_output_args_from(args, newargs.sampling_args);
newargs.decomposition_args = get_decomposition_args_from(args);

warning('off', 'ASYM:unexpectedArgument');
[X, y] = generate_training_data(newargs);
warning('on', 'ASYM:unexpectedArgument');

if false
    % Test with real data
    args.image_type = 'real';
    args.image_root = [asymmetryroot,'data/retinograms/drive/training'];
    sampling_args = get_sampling_args_from(args);
    image_lists = create_image_lists(sampling_args);
    image_lists(1)
    warning('off', 'ASYM:unexpectedArgument');
    [X, y] = generate_training_data(args);
    warning('on', 'ASYM:unexpectedArgument');
end