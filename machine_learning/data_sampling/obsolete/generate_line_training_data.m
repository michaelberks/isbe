function [training_data training_labels parameters] = generate_line_training_data(varargin)
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
    'decomp_type', 'dt',...
    'num_bgs', [],...
    'bg_stem', [],...
    'bg_zeros', [],...
    'bg_mask_dir',[],...
    'detection_type', 'detection',...
    'pts_per_image', 500,...
    'bg_ratio', 1,...
    'width_range', [4 16],...
    'orientation_range', [0 360],...
    'contrast_range', [4 8],...
    'decay_rate', 4,...
    'line_type', 'sin',...
    'normalise', 0,...
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
    'plot', 0);
clear varargin;

%Copy arguments used in sampling funtion into new structure
switch args.decomp_type
    case 'dt'
        if 0%length(args.num_levels) == 1
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
    case {'g2d', 'g2di', 'clover', 'haar', 'g1d'}
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
            size_sample_vector = length(sampling_args.levels)*compute_dt_feature_size(sampling_args);
            
        case 'mono'
            size_sample_vector = 3*args.num_levels*args.win_size^2;
            
        case {'g2d','clover'}
            size_sample_vector = 3*length(args.sigma_range)*args.win_size^2;
            
        case {'haar', 'g1d'}
            size_sample_vector = 2*length(args.sigma_range)*args.win_size^2;
            
        case 'g2di'
            size_sample_vector = 3*args.sigma_range(2)*args.win_size^2;
            
        case 'linop'
            size_sample_vector = ...
                (args.num_angles - args.do_max*(args.num_angles-1))*args.num_levels*args.win_size^2;
            
        case 'pixel'
            size_sample_vector = args.win_size^2;
            
        otherwise
            size_sample_vector = length(sampling_args.levels)*compute_dt_feature_size(sampling_args);
    end    
end   

%Workout whether we need to search for the listing of background patches
if isempty(args.bg_stem) || isempty(args.num_bgs)
    %get directory listing of backgrounds
    bg_list = dir([args.bg_dir, '*.mat']);
    args.num_bgs = length(bg_list);
elseif isempty(args.bg_zeros) 
    args.bg_zeros = floor(log10(args.num_bgs)) + 1;
end

%Create memory storage for data
training_data = zeros(args.num_samples, size_sample_vector);
training_labels = zeros(args.num_samples, 1);

curr_sample = 1;
curr_image = 1;

%Create images filling up sample data until desired number of sampled
%reached
while true 
    %reaching num_samples is checked below, loop terminates with a return statement
    
    %
    % 1. Load in a preprocessed mammographic background
    %----------------------------------------------------------------------
    
    %randomly select background - will load a matrix 'bg'
    bg_idx = ceil(args.num_bgs*rand); 
    if exist('bg_list', 'var')
        bg_name = [args.bg_dir bg_list(bg_idx).name];
        bg_mask_name = [args.bg_mask_dir bg_list(bg_idx).name(1:end-4) '_mask.mat'];
    else   
        bg_name = [args.bg_dir args.bg_stem zerostr(bg_idx, args.bg_zeros), '.mat'];
        bg_mask_name = [args.bg_mask_dir args.bg_stem zerostr(bg_idx, args.bg_zeros), '_mask.mat'];
    end
    bg = double(u_load(bg_name));
    
    %get size of background
    [row col] = size(bg);
    
    %
    % 2. Generate a synthetic linear structure to superimpose on the
    % background
    %----------------------------------------------------------------------
    
    %randomly select orientation
    orientation = args.orientation_range(1) +...
        (args.orientation_range(2)-args.orientation_range(1))*rand;
    
    %randomly select width
    %width = args.width_range(1) + (args.width_range(2)-args.width_range(1))*rand;
    mu = (args.width_range(2) - args.width_range(1)) / (2*log(args.decay_rate));
    width = args.width_range(1) + exp_rand_sample(mu);
    
    %randomly select contrast on exponential distribution
    mu = (args.contrast_range(2) - args.contrast_range(1)) / (2*log(args.decay_rate));
    contrast = args.contrast_range(1) + exp_rand_sample(mu);
    
    %randomly select squashiness of the line
    squash = rand;
    
    %randomly select the radial curvature of the line
    radius = 2^(8*rand + 8);
    
    %make  bar
    switch args.line_type
        case 'ellipse'
            [bar_image, label, label_centre, label_orientation] = create_ellipse_bar(width/2, contrast, orientation, row, col, row/2, col/2);
            
        case 'sin'
            [bar_image, label, label_centre, label_orientation] = create_sin_bar(width/2, contrast, orientation, row, col, squash, row/2, col/2);
        
        case 'curve'
            [bar_image, label, label_centre, label_orientation] = create_sin_curve(width/2, contrast, radius, orientation, squash, row, col, row/2, col/2);
            
        otherwise
            warning(['Bar type: ', args.line_type, ' not recognised']); %#ok
            [bar_image, label, label_centre, label_orientation] = create_ellipse_bar(width/2, contrast, orientation, row, col, row/2, col/2);
    end
    
    clear dummy;
    
    %Store bar parameters
    parameters(curr_image).bg_idx = bg_idx; %#ok
    parameters(curr_image).width = width; %#ok
    parameters(curr_image).contrast = contrast; %#ok
    parameters(curr_image).orientation = orientation; %#ok
    parameters(curr_image).squash = squash; %#ok
    parameters(curr_image).radius = radius; %#ok   
    
    %Add bar to background
    image_in = bar_image + bg; clear bar_image bg;
    
    %
    % 3. Select the foreground and background pixels to sample from in this
    % image
    %----------------------------------------------------------------------
    
    %Make an edge mask to avoid edge of image
    if pad_w <1
        pad_w=1;
    end
    
    %Check to see if a mask exists
    if exist(bg_mask_name, 'file')
        edge_mask = u_load(bg_mask_name);
    else
        edge_mask = true(row, col);
    end   
    
    edge_mask([1:pad_w end-pad_w+1:end], :) = false;
    edge_mask(:, [1:pad_w end-pad_w+1:end]) = false;
    
    switch args.detection_type
        case 'detection'
            %Need bg pixels as well
            bar_idx = find(label_centre & edge_mask); %Only use line centre
            bg_idx = find(~label_centre & edge_mask);
            n_pts = length(bar_idx) + length(bg_idx);
        case 'orientation'
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
    num_samples_image = ...
        sample_from_binomial(args.num_samples, args.pts_per_image/args.num_samples, 1);

    %Check we've not selected too many points - either because there aren't
    %that many in the image or we've already reached our target
    num_samples_image = min([num_samples_image n_pts args.num_samples-curr_sample+1]);
    
    if ~num_samples_image
        continue;
    end
    
    switch args.detection_type
        case 'detection'
            
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
            training_labels(curr_sample:num_samples_image+curr_sample-1) = ...
                label_centre(image_idx);

            clear bar_idx bg_idx rand_idx_fg rand_idx_bg
            
        case 'orientation'
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
            training_labels(curr_sample:num_samples_image+curr_sample-1) = ...
                width;

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
     
    %
    % 4. Sample the decomposition coefficients for this image
    %----------------------------------------------------------------------
    
    switch args.decomp_type
        case 'dt' 
            %Compute dual-tree transform of image
            dt = compute_dual_tree(image_in, sampling_args.levels(end), args.use_nag);
            
            %Sample DT coefficients from specified rows and cols according to
            %sampling arguments
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                sample_dt_data(dt, rows, cols, sampling_args);
            clear dt;
            
        case 'mono'
            [local_amp local_phase local_ori] = monogenic(image_in, args.num_levels, args.min_wavelength, 2, args.onf, 1);
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
           clear g2d_responses;
           
        case 'g2di'       
            g2d_responses = compute_gaussian_2nd_derivatives_d(image_in,  args.sigma_range(1), args.sigma_range(2));
            clear image_in;
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                sample_g2d_data_d(g2d_responses, rows, cols, sampling_args);
            clear g2d_responses;
            
        case 'clover'
            clover_responses = compute_clover_responses(image_in, args.sigma_range);
            clear channel;
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
               sample_g2d_data(...
               clover_responses(:,:,:,1),... 
               clover_responses(:,:,:,2),...
               clover_responses(:,:,:,3), rows, cols, sampling_args);
           clear clover_responses;
           
        case 'haar'
            haar_responses = compute_haar_responses(image_in, args.sigma_range);
            clear channel;
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
               sample_haar_data(...
               haar_responses(:,:,:,1),... 
               haar_responses(:,:,:,2), rows, cols, sampling_args);
           clear haar_responses;
       case 'g1d'
            g1d_responses = compute_gaussian_1st_derivatives(image_in,  args.sigma_range);
            clear channel;
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
               sample_g1d_data(...
               g1d_responses(:,:,:,1),... 
               g1d_responses(:,:,:,2), rows, cols, sampling_args);
           clear g1d_responses;
            
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
            dt = compute_dual_tree(image_in, sampling_args.levels(end), args.use_nag);
            
            %Sample DT coefficients from specified rows and cols according to
            %sampling arguments
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                sample_dt_data(dt, rows, cols, sampling_args);
            clear dt;
    end       

    %
    % 5. Update counters and check termination criteria
    %----------------------------------------------------------------------
    
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
            end
            par_list = dir([args.save_path '*parameters*']);
            tree_num = length(par_list) + 1;
            save([args.save_path 'parameters' zerostr(tree_num, 3), '.mat'], 'parameters');
        end
        return;
    end
    
    %Otherwise increment the image count and continue
    curr_image = curr_image + 1;   
    
end
