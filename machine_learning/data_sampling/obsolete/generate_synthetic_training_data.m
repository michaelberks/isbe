function [training_data training_labels parameters] = generate_training_data(varargin)
%
% GENERATE_LINE_TRAINING_DATA - 
% USAGE:
%
% Inputs:
%      num_samples - the total number of samples
%
%
% Outputs:
%
% Example:
%
% Notes:
% See also:
%
% Created: 08-February-2010
% Author: Michael Berks
% Email : michael.berks@manchester.ac.uk
% Phone : +44 (0)161 275 7669
% Copyright: (C) University of Manchester
%

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    {'num_samples', ...
    'bg_dir'},... % the mandatory arguments
	'image_type','line',...
    'decomp_type', 'dt',...
    'num_bgs', [],...
    'bg_stem', [],...
    'bg_zeros', [],...
	'bg_fmt', 'mat',...
    'detection_type', 'detection',...
    'pts_per_image', 500,...
    'bg_ratio', 1,...
    'width_range', [4 16],...
    'orientation_range', [0 360],...
    'contrast_range', [4 8],...
    'decay_rate', 4,...
    'line_type', 'sin',...
    'normalise', 1,...
    'num_levels', 5,...
    'feature_shape', 'rect',...
    'feature_type', 'all',...
    'do_max', 0,...
    'rotate', 0,...
    'win_size', 3,...
    'num_angles', 8,...
    'min_wavelength', 4,...
    'onf', 0.65,...
    'sigma_range', [1 2 4 8],...
    'pca', [],...
    'use_nag', 1,...
    'save_path', [], ...
    'plot', 0,...
    'quiet', 1);
clear varargin;

%Copy arguments used in sampling funtion into new structure
switch args.decomp_type
    case 'dt'
        if length(args.num_levels) == 1
            sampling_args.levels = 1:args.num_levels;
        else
            sampling_args.levels = args.num_levels;
        end
        sampling_args.feature_shape = args.feature_shape;
        sampling_args.feature_type = args.feature_type;
        sampling_args.do_max = args.do_max;
        sampling_args.rotate = args.rotate;
        sampling_args.win_size = args.win_size;
    case 'mono'
        sampling_args.win_size = args.win_size;
    case {'g2d','clover','haar'}
        sampling_args.win_size = args.win_size;
    case 'linop'
        sampling_args.win_size = args.win_size;
        sampling_args.num_levels = args.num_levels;
        sampling_args.num_angles = args.num_angles;
        sampling_args.do_max = args.do_max;
        sampling_args.rotate = args.rotate;
    case 'pixel'
        sampling_args.win_size = args.win_size;
    otherwise
        warning(['decomposition type: ', args.decomp_type, ' not recognised, using DT-CWT']); %#ok
        sampling_args.feature_shape = args.feature_shape;
        sampling_args.feature_type = args.feature_type;
        sampling_args.do_max = args.do_max;
        sampling_args.rotate = args.rotate;
        sampling_args.win_size = args.win_size;
end

%Constants computed from arguments
pad_w = floor(args.win_size/2); %half size of window size

%work out size of feature vectors

%See if principal components are specified, and if so, load them. This
%defines the sample feature length
if ~isempty(args.pca);
    %pca should have fields mean and modes
    size_sample_vector = size(args.pca.modes,2);
    sampling_args.pca = args.pca; args = rmfield(args, 'pca');
    
else %Workout feature length from sampling arguments
    
    switch args.decomp_type
        case 'dt'
            size_sample_vector = args.num_levels*compute_dt_feature_size(sampling_args);
            
        case 'mono'
            size_sample_vector = 3*args.num_levels*args.win_size^2;
            
        case {'g2d','clover'}
            size_sample_vector = 3*length(args.sigma_range)*args.win_size^2;
            
        case 'haar'
            size_sample_vector = 2*length(args.sigma_range)*args.win_size^2;
            
        case 'linop'
            size_sample_vector = ...
                (args.num_angles - args.do_max*(args.num_angles-1))*args.num_levels*args.win_size^2;
            
        case 'pixel'
            size_sample_vector = args.win_size^2;
            
        otherwise
            size_sample_vector = args.num_levels*compute_dt_feature_size(sampling_args);
    end    
end   

%Workout whether we need to search for the listing of background patches
if isempty(args.bg_stem) || isempty(args.num_bgs)
    %get directory listing of backgrounds
    bg_list = dir([args.bg_dir,'*.',args.bg_fmt]);
    args.num_bgs = length(bg_list);
elseif isempty(args.bg_zeros) 
    args.bg_zeros = floor(log10(args.num_bgs)) + 1;
end

%Create memory storage for data
training_data = zeros(args.num_samples, size_sample_vector);
training_labels = zeros(args.num_samples, 1);

curr_sample = 1;
curr_image = 1;

%Create images filling up sample data until desired number of samples reached
if ~args.quiet; tb = timebar('limit',args.num_samples,'title','Sampling data'); end
while true 
    %reaching num_samples is checked below, loop terminates with a return statement
    
%% Generate a synthetic image to sample for training data
    
		switch args.image_type
			case 'line',
				% Load in a preprocessed mammographic background

				% randomly select background - will load a matrix 'bg'
				bg_idx = ceil(args.num_bgs*rand);
				if exist('bg_list', 'var')
					bg_filename = [args.bg_dir bg_list(bg_idx).name];
				else   
					bg_filename = [args.bg_dir args.bg_stem zerostr(bg_idx, args.bg_zeros),'.',args.bg_fmt];
				end
				
				% if background not found then skip and try another
				if ~exist(bg_filename,'file')
					warning([bg_filename,' not found']);
					continue;
				end

				% load background according to specified format
				switch args.bg_fmt
					case 'mat', bg = u_load(bg_filename);
					case 'png', bg = double(imread(bg_filename));
				end									

				% generate the image
				[image_in,label,label_centre,label_orientation,params] = ...
					generate_line_image(bg,args);
				
				% Store bar parameters
				params.bg_idx = bg_idx; %#ok
				parameters(curr_image) = params; %#ok
				
			case 'grain',
				% generate a 'grain' image
				[image_in,label,label_centre,label_orientation,params] = ...
					generate_grain_image(args);
				
				% Store grain parameters
				parameters(curr_image) = params;
				
			otherwise,
				% do something
		end

    %get size of background
    [row col] = size(image_in);
	
%% Select the foreground and background pixels to sample from in this image
    
    %Make an edge mask to avoid edge of image
    if pad_w <1
        pad_w=1;
    end
    edge_mask = true(row, col);
    edge_mask([1:pad_w end-pad_w+1:end], :) = false;
    edge_mask(:, [1:pad_w end-pad_w+1:end]) = false;
    
    switch args.detection_type
        case {'detection','logistic_classification'}
            %Need bg pixels as well
            bar_idx = find(label_centre & edge_mask); %Only use line centre
            bg_idx = find(~label_centre & edge_mask);
            n_pts = length(bar_idx) + length(bg_idx);
        case {'orientation','linear_regression','logistic_regression','boosted_regression'}
            %Use label to get indices for this image
            bar_idx = find(label & edge_mask); %Use whole line width
            n_pts = length(bar_idx);
        case 'width'
            %Only use line centre
            bar_idx = find(label_centre & edge_mask); 
            n_pts = length(bar_idx);
    end
    clear edge_mask;
    
    %Compute the number of samples to take from this image    
    num_samples_image = sample_from_binomial(args.num_samples,...
							args.pts_per_image/args.num_samples, 1);

	%Check we've not selected too many points - either because there aren't
    %that many in the image or we've already reached our target
    num_samples_image = min([num_samples_image n_pts args.num_samples-curr_sample+1]);

	if ~num_samples_image
        continue;
    end
    
    switch args.detection_type
        case {'detection','logistic_classification'}
            
            %Workout number of background and foreground pixels required
            num_fg = floor(num_samples_image / (1 + args.bg_ratio));
            num_bg = ceil(num_samples_image * args.bg_ratio / (1 + args.bg_ratio));

            %Take a random sample of the foreground pixels
            rand_idx_fg = randperm(length(bar_idx));
            bar_idx = bar_idx(rand_idx_fg(1:num_fg));

            %Now take a random sample of the background pixels
            rand_idx_bg = randperm(length(bg_idx));
            bg_idx = bg_idx(rand_idx_bg(1:num_bg)); 
            clear rand_idx_bg rand_idx_fg;

            %Combine the foreground/background samples
            image_idx = [bar_idx; bg_idx];

            %Save the sample labels for this image
            training_labels(curr_sample:num_samples_image+curr_sample-1) = label_centre(image_idx);

            clear bar_idx bg_idx rand_idx_fg rand_idx_bg
            
        case {'orientation','linear_regression','logistic_regression','boosted_regression'}
            %Take required random sample of indices
            rand_idx = randperm(n_pts);
            image_idx = bar_idx(rand_idx(1:num_samples_image));

            %Save the sample labels for this image
            training_labels(curr_sample:num_samples_image+curr_sample-1) = ...
                complex(...
                    cosd(2*label_orientation(image_idx)),...
                    sind(2*label_orientation(image_idx)));

            clear bar_idx rand_idx
            
        case 'width'
            %Take required random sample of indices
            rand_idx = randperm(n_pts);
            image_idx = bar_idx(rand_idx(1:num_samples_image));

            %Save the sample labels for this image
            training_labels(curr_sample:num_samples_image+curr_sample-1) = width;

            clear bar_idx rand_idx
	end
    
    clear label label_centre label_orientation;
    
    %Convert to row/col
    [rows cols] = ind2sub([row col], image_idx);
    
    %If asked display the sampling locations
    %display(['Num samples in image ' num2str(curr_image) ': ' num2str(num_samples_image)]);
    if args.plot
        figure; imagesc(label); axis image; colormap(gray(256)); hold on;
        plot(cols, rows, 'rx');
    end
    
    %Normalise the image if necessary
    if args.normalise
        im_mean = mean(image_in(:));
        im_std = std(image_in(:));
        image_in = (image_in - im_mean) / im_std;
    end
     
%% Sample the decomposition coefficients for this image
    
    switch args.decomp_type
        case 'dt' 
            %Compute dual-tree transform of image
            dt = compute_dual_tree(image_in, args.num_levels, args.use_nag);
            
            %Sample DT coefficients from specified rows and cols according to
            %sampling arguments
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                sample_dt_data(dt, rows, cols, sampling_args);
			clear dt;
            
        case 'mono'
            [local_amp local_phase local_ori] = ...
							monogenic(image_in, args.num_levels, args.min_wavelength, 2, args.onf, 1);
            clear image_in;
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                sample_monogenic_data(local_amp, local_phase, local_ori, rows, cols, sampling_args);
            clear local_amp local_phase;
            
        case 'g2d'
            g2d_responses = compute_gaussian_2nd_derivatives(image_in,  args.sigma_range);
            clear image_in;
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                sample_g2d_data(...
                g2d_responses(:,:,:,1),... 
                g2d_responses(:,:,:,2),...
                g2d_responses(:,:,:,3), rows, cols, sampling_args);
            clear local_amp local_phase;
			
		case 'clover'
            clover_responses = compute_clover_responses(image_in, args.sigma_range);
            clear image_in;
            training_data(curr_sample:curr_sample+num_samples_image-1,:) = ...
                sample_g2d_data(...
                clover_responses(:,:,:,1),... 
                clover_responses(:,:,:,2),...
                clover_responses(:,:,:,3), rows, cols, sampling_args);
            clear local_amp local_phase;
            
		case 'haar'
            haar_responses = compute_haar_responses(image_in, args.sigma_range);
            clear image_in;
            training_data(curr_sample:curr_sample+num_samples_image-1,:) = ...
                sample_haar_data(...
                haar_responses(:,:,:,1),... 
                haar_responses(:,:,:,2), rows, cols, sampling_args);
            clear local_amp local_phase;

		case 'linop'
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                sample_linop_data(image_in, rows, cols, sampling_args);
            clear image_in;
            
        case 'pixel'
            %Sample pixels from image
            training_data(curr_sample:num_samples_image+curr_sample-1, :) = ...
                    sample_pixel_data(image_in, rows, cols, sampling_args);
            
        otherwise
            %Compute dual-tree transform of image
            dt = compute_dual_tree(image_in, args.num_levels, args.use_nag);
            
            %Sample DT coefficients from specified rows and cols according to
            %sampling arguments
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                sample_dt_data(dt, rows, cols, sampling_args);
            clear dt;
    end       

%%  Update counters and check termination criteria
    
    %Increment the sample count
    curr_sample = curr_sample + num_samples_image;
    
    %Workout whether we've reached the diseased number of samples
    if curr_sample > args.num_samples
        %if so, save the sampling parameters the return
        if ~isempty(args.save_path)
            if ~(strcmp(args.save_path(end), '/') || strcmp(args.save_path(end), '\'));  
                args.save_path = [args.save_path '/'];
            end
            if ~exist(args.save_path, 'dir')
                mkdir(args.save_path);
				if ~ispc
					fileattrib(args.save_path,'+w','g');
				end
            end
            par_list = dir([args.save_path '*parameters*']);
            tree_num = length(par_list) + 1;
            save([args.save_path 'parameters' zerostr(tree_num, 3), '.mat'], 'parameters');
		end
        break;
    end
    
    %Otherwise increment the image count and continue
    curr_image = curr_image + 1; 
		
	% update timebar
	if ~args.quiet; timebar(tb,'advance',num_samples_image); end
end
if ~args.quiet; timebar(tb,'close'); clear tb; end
