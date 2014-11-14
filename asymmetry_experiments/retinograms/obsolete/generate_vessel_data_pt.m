function [training_data training_labels] = generate_vessel_data(varargin)
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
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    {'num_samples', ...
    'image_dir',...
    'foveal_mask_dir',...
    'vessel_mask_dir'},... % the mandatory arguments
    'prediction_type', 'detection',...
    'ori_dir', [],...
    'width_dir', [],...
    'selected_images', [], ...
    'rgb_channel', 'g',...
    'win_size', 3,...
    'num_levels', 4,...
    'feature_type', 'conj',...
    'feature_shape', 'rect',...
    'rotate', 0,...
    'do_max', 0,...
    'decomp_type', 'dt',...
    'bg_ratio', 1,...
    'num_angles', 8,...
    'min_wavelength', 4,...
    'onf', 0.65,...
    'sigma_range', [1 2 4 8],...
    'pca', [],...
    'use_nag', 0);
clear varargin;

% total number of feature vectors to sample
num_samples = args.num_samples; 

% Copy arguments used in sampling funtion into new structure
sampling_args = get_sampling_args_from(args);

% See if principal components are specified, and if so, load them. 
% This defines the sample feature length
if ~isempty(args.pca);
    % pca should have fields mean and modes
    size_sample_vector = size(args.pca.modes,2);
    sampling_args.pca = args.pca; 
	args = rmfield(args, 'pca');
else
	% Workout feature length from sampling arguments
    samples_per_channel = get_samples_per_channel(sampling_args);
    [num_channels, channels] = get_channels(args);
	size_sample_vector = num_channels * samples_per_channel;    
end

% Create memory storage for data
training_data = zeros(args.num_samples, size_sample_vector);
training_labels = zeros(args.num_samples, 1);

% Get list of all the images, foveal masks and vessel masks - we're going to
% assume these all match up for now
imlists = create_vessel_imlists(args);
num_images = length(imlists);

% Loop through each image sampling data
curr_sample = 1;
for kk = 1:num_images
    % Work out the number of samples to take
    N = num_samples + 1 - curr_sample;
    p = 1 / (num_images + 1 - kk);
    num_samples_image = sample_from_binomial(N, p, 1);
    
    if (num_samples_image == 0)
        continue;
    end
        
    % Get (possibly weighted) foreground and background pixel sets
    [fg_mask, bg_mask] = get_vessel_masks(imlists(kk), args);
   
    % Check we have enough foreground samples to choose from
    total_v_pts = sum(fg_mask(:));
    num_samples_image = min(num_samples_image, total_v_pts); %just in case...
    display(['Samples in image ' num2str(kk) ': ', num2str(num_samples_image)]);

    % Choose sample positions in foreground and background
    [fg_idx, bg_idx] = sample_vessel_labels(fg_mask, bg_mask, ...
                                            num_samples_image, args);

    %
    % Get relevant output data for given prediction type
    %----------------------------------------------------------------------
    switch args.prediction_type
        case {'detection', 'centre_detection', 'logistic_classification', ...
              'foveal_detection', ...
              'junction_detection', 'junction_centre_detection'}
            % Combine the vessel/background samples
            labels = [true(numel(fg_idx), 1); false(numel(bg_idx), 1)];
            
        case {'orientation', 'linear_regression', 'logistic_regression', ...
              'boosted_regression', 'centre_orientation'}
            % Load in orientation_map and save sample labels
            ori_map = u_load([args.ori_dir imlists(kk).ori_gt]);
            labels = ori_map(fg_idx);
            
        case {'mixed_orientation', 'mixed_centre_orientation'}
            %now randomly select one of the mixed orientations for each
            %point
            %DEBUG: figure; imgray(vessel_mask);
            	
            %Load in orientation ground truth
            ori_map = u_load([args.ori_dir imlists(kk).ori_gt]);
            ori_sample = zeros(num_samples_image,1);
            for ii = 1:num_samples_image
                mixed_oris = ori_gt.mixed_oris{ori_gt.mixed_idx(image_idx(ii))};
                ori_sample(ii) = mixed_oris(ceil(rand*length(mixed_oris)));
                
                %DEBUG: [y x] = ind2sub(size(vessel_mask), image_idx(ii));
                %DEBUG: quiver(ones(length(mixed_oris),1)*x,ones(length(mixed_oris),1)*y,real(mixed_oris),-imag(mixed_oris), 'color','r');
                %DEBUG: quiver(ones(length(mixed_oris),1)*x,ones(length(mixed_oris),1)*y,-real(mixed_oris),imag(mixed_oris), 'color','r');
            end
            
            %Save the sample orientations for this image
            labels = ori_sample;
            
        case {'width'}
            % Load in width map and sample at points
            width_map = u_load([args.width_dir imlists.width]);
            labels = width_map(fg_idx);
            
        otherwise
            if strcmp(get_username, 'ptresadern')
                error(['Unknown output type: ',args.prediction_type]);
            end
    end   
    
    % Compute last row of training data matrices (first row is curr_sample)
    last_sample = curr_sample + num_samples_image-1;

    % Save the sample labels for this image
    training_labels(curr_sample:last_sample) = labels;

    
    %
    % 4. Sample the decomposition coefficients for each selected channel
    %----------------------------------------------------------------------
    
    %Load in retinogram
    ret = u_load([args.image_dir imlists(kk).raw]);
    
    %Convert to row/col
    [rows cols] = ind2sub([size(ret,1) size(ret,2)], [fg_idx; bg_idx]);

    for ch = 1:num_channels
        %select the desired channel
        if strcmpi(args.rgb_channel, 'rgb')
            channel = rgb2gray(ret);
        else
            channel = ret(:,:,channels(ch));
        end
        
        % Get training feature vectors
        features = sample_image_features(channel, rows, cols, sampling_args);

        % Store with rest of training data
        data_cols = (1:samples_per_channel) + (ch-1)*samples_per_channel;
        training_data(curr_sample:last_sample, data_cols) = features;
    end
    
    %Update the current sample count
    curr_sample = curr_sample + num_samples_image;
end


