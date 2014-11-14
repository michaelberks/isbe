function [training_data training_labels parameters] = generate_dt_rotation_coeffs(varargin)
%
% SAMPLE_TRAINING_DATA_ALLSUBBANDS create training data which use all
% subbands coefficients of DTCWT
% USAGE:
% [training_data training_labels parameters] =...
%   sample_training_data_allsubbands('num_samples', 1000, 'bg_dir',...
%   'E:\DTCWTmini\data\normal_smooth512\', 'save_path',
%   'E:\DTCWTmini\data\synthetic_512\', 'plot', 1);
%
% Inputs:
%      num_samples - the total number of samples
%      bg_dir -  the patch of background patch
%      halfwidth - halfwidth of Gaussian profile at half its maximum height
%
%      contrast - maximum height of (scaled) Gaussian profile
%
%      orientation - orientation of bar in image in degrees
%
%      row - number of rows in image
%
%      col - number of columns in image
%
%
% Outputs:
%      image_out - image containing spicule
%
%      label - label of spicule (1) vs background (0)
%
%      label_centre - the label of centre line (1) vs background (0)
%
% Example:
%
% Notes:
% See also: sample_bar_training_data
%
% Created: 08-February-2010
% Author: Zezhi Chen and Michael Berks
% Email : zezhi.chen@manchester.ac.uk
% Phone : +44 (0)161 275 7669
% Copyright: (C) University of Manchester
%


% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    {'num_samples', ...
    'bg_dir'},... % the mandatory arguments
    'num_bgs', [],...
    'bg_stem', [],...
    'bg_zeros', [],...
    'bg_ratio', 1,...
    'fg_ratio', 0.5,...
    'width_range', [4 16],...
    'orientation_range', [0 180],...
    'contrast_range', [4 16],...
    'decay_rate', 4,...
    'line_type', 'sin',...
    'num_levels', 5,...
    'feature_shape', 'rect',...
    'feature_type', 'all',...
    'do_max', 0,...
    'rotate', 0,...
    'win_size', 3,...
    'pca', [],...
    'use_nag', 1,...
    'save_path', [], ...
    'plot', 0);
clear varargin;

%Copy arguments used in sampling funtion into new structure
sampling_args.feature_shape = args.feature_shape;
sampling_args.feature_type = args.feature_type;
sampling_args.do_max = args.do_max;
sampling_args.rotate = args.rotate;
sampling_args.win_size = args.win_size;


%Constants computed from arguments
pad_w = floor(args.win_size/2); %half size of window size

%work out size of feature vectors

%See if principal components are specified, and if so, load them. This
%defines the sample feature length
if ~isempty(args.pca);
    %pca should have fields mean and modes
    pca = u_load(args.pca);
    size_sample_vector = size(pca.modes,2);
    sampling_args.pca = pca;
    
else %Workout feature length from sampling arguments
    
    if args.do_max
        num_bands = 1;
    else
        num_bands = 6;
    end
    switch args.feature_type
        case 'all'
            num_features = 2;

        case {'real', 'mag', 'phase', 'conj'}    
            num_features = 1;

        otherwise
            warning(['Feature type: ', args.feature_type, ' not recognised, using phase and magnitude (feature_type = ''all'')']); %#ok
            num_features = 2;
    end
    switch args.feature_shape
        case 'rect'
            size_sample_vector = args.win_size*args.win_size*args.num_levels*num_bands*num_features;

        case 'clock'
            switch args.feature_type
                case {'real', 'phase'}    
                    size_sample_vector = 84*args.num_levels;
                    
                otherwise %inlcudes all, mag
                    size_sample_vector = (84*num_features - 6)*args.num_levels;
            end
        
        otherwise
            warning(['Feature shape: ', args.feature_shape, ' not recognised, using square windows (feature_type = ''rect'')']); %#ok
            size_sample_vector = args.win_size*args.win_size*args.num_levels*num_bands*num_features;
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
        bg = u_load([args.bg_dir bg_list(bg_idx).name]);
    else   
        bg = u_load([args.bg_dir args.bg_stem zerostr(bg_idx, args.bg_zeros), '.mat']);
    end
    
    %get size of background
    [row col] = size(bg);
    
    %
    % 2. Generate a synthetic linear structure to superimpose on the
    % background
    %----------------------------------------------------------------------
    
    %randomly select orientation between 0 and 180 degrees
    orientation = args.orientation_range(1) +...
        (args.orientation_range(2)-args.orientation_range(1))*rand;
    
    %randomly select width between 4 and 32
    width = args.width_range(1) + (args.width_range(2)-args.width_range(1))*rand;
    
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
            [bar_image, dummy, label] = create_ellipse_bar(width/2, contrast, orientation, row, col, row/2, col/2);
            
        case 'sin'
            [bar_image, dummy, label] = create_sin_bar(width/2, contrast, orientation, row, col, squash, row/2, col/2);
        
        case 'circle'
            [bar_image, dummy, label] = create_sin_curve(width/2, contrast, radius, orientation, squash, row, col, row/2, col/2);
            
        otherwise
            warning(['Bar type: ', args.line_type, ' not recognised']); %#ok
            [bar_image, dummy, label] = create_ellipse_bar(width/2, contrast, orientation, row, col, row/2, col/2);
    end
    
    clear dummy;
    
    %Store bar parameters
    parameters(curr_image).bg_idx = bg_idx; %#ok
    parameters(curr_image).width = width; %#ok
    parameters(curr_image).contrast = contrast; %#ok
    parameters(curr_image).orientation = orientation; %#ok
    parameters(curr_image).squash = squash; %#ok
    
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
    edge_mask = true(row, col);
    edge_mask([1:pad_w end-pad_w+1], :) = false;
    edge_mask(:, [1:pad_w end-pad_w+1]) = false;
    
    %Use label to get indices for this image
    bar_idx = find(label & edge_mask); %use all bar pixels
    bg_idx = find(~label & edge_mask);
    clear edge_mask;
    
    %Now select a random selection of the background pixels
    rand_idx_fg = randperm(length(bar_idx));
    bar_idx = bar_idx(rand_idx_fg(1:round(length(rand_idx_fg)*args.fg_ratio)));
    
    rand_idx = randperm(length(bg_idx));
    bg_idx = bg_idx(rand_idx(1:length(bar_idx)*args.bg_ratio)); clear rand_idx rand_idx_fg;
    
    %Get indices for this image and convert to row,col
    image_idx = [bar_idx; bg_idx];
    num_samples_image = length(image_idx);
    [rows cols] = ind2sub([row col], image_idx);
    
    %If asked display the sampling locations
    %display(['Num samples in image ' num2str(curr_image) ': ' num2str(num_samples_image)]);
    if args.plot
        figure; imagesc(label); axis image; colormap(gray(256)); hold on;
        plot(cols(1:length(bar_idx)), rows(1:length(bar_idx)), 'rx');
        plot(cols(length(bar_idx)+1:end), rows(length(bar_idx)+1:end), 'gx');
    end
    
    %Save the sample labels for this image
    training_labels(curr_sample:num_samples_image+curr_sample-1)= label(image_idx);
    clear bar_idx bg_idx image_idx label;
    
    
    %
    % 4. Sample the DT-CWT coefficients for this image
    %----------------------------------------------------------------------
    
    % Create DT-CWT of image
    dual_tree = dtwavexfm2b(image_in, args.num_levels); clear image_in;
    
    %Try computing NAG interpolation - this may not work eg. if NAG libraries
    %aren't installed
    if args.use_nag
        try 
            [dt.knot_mag dt.knot_im dt.knot_re dt.dt_dims] = dt_interp_nag_in(dual_tree);
        catch
            err = lasterror;
            display('Probably executing NAG interpolation. Will use matlab interpolation');
            display(err.message);
            dt = dual_tree;
        end
    else
        dt = dual_tree;
    end
    clear dual_tree;

    %Sample DT coefficients from specified rows and cols according to
    %sampling arguments
    training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
        sample_dt_data(dt, rows, cols, sampling_args);

    %Workout whether we've reached the diseased number of samples
    if curr_sample >= args.num_samples
        %if so, save the sampling parameters the return
        if ~isempty(args.save_path)
            if ~(strcmp(args.save_path(end), '/') || strcmp(args.save_path(end), '\'));  
                args.save_path = [args.save_path '/']; %#ok
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
    
    %Otherwise increment the sample count and image count and continue
    curr_sample = curr_sample + num_samples_image;
    curr_image = curr_image + 1;   
    
end
