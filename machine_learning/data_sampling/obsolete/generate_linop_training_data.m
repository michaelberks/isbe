function [training_data training_labels parameters] = generate_linop_training_data(varargin)
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
    'detection_type', 'detection',...
    'win_size', 3,...
    'num_levels', 6,...
    'num_angles', 5,...
    'do_max', 0,...
    'pts_per_image', 400,...
    'bg_ratio', 1,...
    'width_range', [4 16],...
    'orientation_range', [0 360],...
    'contrast_range', [4 16],...
    'decay_rate', 4, ...
    'line_type', 'sin',...
    'save_path', [], ...
    'plot', 0);
clear varargin;

win_size = args.win_size; % the size of scanning window
num_levels = args.num_levels; %number of levels of DT-CWT from which to extract coefficients
num_angles = args.num_angles; %number of levels of DT-CWT from which to extract coefficients

%Constants computed from arguments
pad_w = floor(win_size/2); %half size of window size
win_idx = -pad_w:pad_w;
do_orientation = strcmpi(args.detection_type, 'orientation');

%work out size of feature vectors
if args.do_max
    size_sample_vector = win_size*win_size*num_levels;
else
    size_sample_vector = win_size*win_size*num_levels*num_angles;
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
            [bar_image, label, label_centre, label_orientation] = create_ellipse_bar(width/2, contrast, orientation, row, col, row/2, col/2);
            
        case 'sin'
            [bar_image, label, label_centre, label_orientation] = create_sin_bar(width/2, contrast, orientation, row, col, squash, row/2, col/2);
        
        case 'curve'
            [bar_image, label, label_centre, label_orientation] = create_sin_curve(width/2, contrast, radius, orientation, squash, row, col, row/2, col/2);
            
        otherwise
            warning(['Bar type: ', args.line_type, ' not recognised']); %#ok
            [bar_image, label, label_centre, label_orientation] = create_sin_bar(width/2, contrast, orientation, row, col, squash, row/2, col/2);
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
    edge_mask([1:pad_w end-pad_w+1:end], :) = false;
    edge_mask(:, [1:pad_w end-pad_w+1:end]) = false;
    
    if do_orientation
        %Use label to get indices for this image
        bar_idx = find(label & edge_mask); %Use whole line width
        n_pts = length(bar_idx);
        
    else %Need bg pixels as well
        bar_idx = find(label_centre & edge_mask); %Only use line centre
        bg_idx = find(~label_centre & edge_mask);
        n_pts = length(bar_idx) + length(bg_idx); 
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
    
    if do_orientation
        %Take required random sample of indices
        rand_idx = randperm(n_pts);
        image_idx = bar_idx(rand_idx(1:num_samples_image));
        
        %Save the sample labels for this image
        training_labels(curr_sample:num_samples_image+curr_sample-1) = label_orientation(image_idx);
        
        clear bar_idx rand_idx
    else
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
    end
    clear label label_centre label_orientation;
    
    %Convert to row/col
    [rows cols] = ind2sub([row col], image_idx); clear image_idx;
    
    %If asked display the sampling locations
    %display(['Num samples in image ' num2str(curr_image) ': ' num2str(num_samples_image)]);
    if args.plot
        figure; imagesc(label); axis image; colormap(gray(256)); hold on;
        plot(cols, rows, 'rx');
    end
     
    %
    % 4. Sample the linop coefficients for this image
    %----------------------------------------------------------------------
    
    %Make copies of sample rows and cols at positions of local window patch
    rr = repmat(rows*ones(1,win_size) + ones(num_samples_image,1)*win_idx, 1, win_size);
    cc = kron(cols*ones(1,win_size) + ones(num_samples_image,1)*win_idx, ones(1,win_size));
    
    %Now get linop samples from the image
    linop_samples = line_operator_octave_subset(image_in, num_angles, num_levels, rr, cc); clear image_in rr cc;
    
    if args.do_max
        %get the maximum response across orientations
        linop_samples = squeeze(max(linop_samples, [], 3)); clear dt;
    end
    
    %save sample into main training data
    training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
        reshape(linop_samples, num_samples_image, []);    
    clear linop_samples;
    
    %Increment the sample count
    curr_sample = curr_sample + num_samples_image;
    
    %Workout whether we've reached the diseased number of samples
    if curr_sample > args.num_samples
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
    
    %Otherwise increment the image count and continue
    curr_image = curr_image + 1;   
end
