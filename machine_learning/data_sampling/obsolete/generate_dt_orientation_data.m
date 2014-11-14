function [training_data training_oris parameters] = generate_dt_orientation_data(varargin)
%GENERATE_DT_ORIENTATION_DATA *Insert a one line summary here*
%   [] = sample_bar_training_data(varargin)
%
% GENERATE_DT_ORIENTATION_DATA uses the U_PACKARGS interface function
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
% Created: 28-Jan-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    {'num_samples', ...
    'bg_dir'},... % the mandatory arguments
    'num_bgs', [],...
    'bg_stem', [],...
    'bg_zeros', [],...
    'pts_per_image', 200,...
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
        case {'all', 'conj'}
            num_features = 2;

        case {'real', 'mag', 'phase'}    
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
training_oris = zeros(args.num_samples, 1);

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
            [bar_image label, label_centre, label_orientation] = create_ellipse_bar(width/2, contrast, orientation, row, col, row/2, col/2);
            
        case 'sin'
            [bar_image, label, label_centre, label_orientation] = create_sin_bar(width/2, contrast, orientation, row, col, squash, row/2, col/2);
        
        case 'circle'
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
    bar_idx = find(label & edge_mask); clear edge_mask;
    n_pts = length(bar_idx);
    
    %Compute the number of samples to take from this image    
    num_samples_image = ...
        sample_from_binomial(args.num_samples, args.pts_per_image/args.num_samples, 1);

    %Check we've not selected too many points - either because the aren't
    %that many in the image or we've already reached our target
    num_samples_image = min([num_samples_image n_pts args.num_samples-curr_sample +1]);
    
    if ~num_samples_image
        continue;
    end
    
    %Take required random sample of indices and convert to row/col
    rand_idx = randperm(n_pts);
    image_idx = bar_idx(rand_idx(1:num_samples_image));
    [rows cols] = ind2sub([row col], image_idx);

    if args.plot
        figure; imagesc(label); axis image; colormap(gray(256)); hold on;
        plot(cols, rows, 'rx');
    end

    %Save the orientation of line in the labels container
    training_oris(curr_sample:num_samples_image+curr_sample-1) = label_orientation(image_idx);
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