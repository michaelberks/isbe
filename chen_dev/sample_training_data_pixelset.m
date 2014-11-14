function [training_data training_labels parameters] = sample_training_data_pixelset(varargin)
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
    'win_size', 3,...
    'num_levels', 6,...
    'do_max', 0,...
    'num_bgs', [],...
    'bg_stem', [],...
    'bg_zeros', [],...
    'bg_ratio', 1,...
    'fg_ratio', 0.5,...
    'width_range', [4 16],...
    'orientation_range', [0 180],...
    'contrast_range', [4 16],...
    'decay_rate', 4,...
    'bar_type', 'ellipse',...
    'feature_type', 'all',...
    'save_path', [], ...
    'plot', 0);

win_size = args.win_size; % the size of scanning window
num_levels = args.num_levels; %number of levels of DT-CWT from which to extract coefficients
num_samples = args.num_samples; %total number of feature vectors to sample
bg_ratio = args.bg_ratio; %ration of backgorund pixels to bar pixels
fg_ratio = args.fg_ratio; %ration of sample pixels to the bar pixels
parameters_dir = args.save_path;

%Constants computed from arguments
pad_w = floor(win_size/2); %half size of window size
win_idx = -pad_w:pad_w;

%work out size of feature vectors
if args.do_max
    num_bands = 1;
else
    num_bands = 6;
end
switch args.feature_type
    case {'all', 'rotate'}
        size_sample_vector = 2*win_size*win_size*num_levels*num_bands;

    case {'real', 'mag', 'phase'}
        size_sample_vector = win_size*win_size*num_levels*num_bands;
        
    case {'clock', 'clock_rotate'}
        size_sample_vector = 162*num_levels;
        discard_idx = repmat(logical([ones(1,6) zeros(1,78)]), 1, num_levels);

    otherwise
        warning(['Feature type: ', args.feature_type, ' not recognised, using phase and magnitude (feature_type = ''all'')']); %#ok
        size_sample_vector = 2*win_size*win_size*num_levels;
end

%
if isempty(args.bg_stem) || isempty(args.num_bgs)
    %get directory listing of backgrounds
    bg_list = dir([args.bg_dir, '*.mat']);
    args.num_bgs = length(bg_list);
elseif isempty(args.bg_zeros) 
    args.bg_zeros = floor(log10(args.num_bgs)) + 1;
end

%Create memory storage for data
training_data = zeros(num_samples, size_sample_vector);
training_labels = zeros(num_samples, 1);

curr_sample = 1;
curr_image = 1;

%Create images filling up sample data until training data full
while true %reaching num_samples is checked below, loop terminates with a return statement
    
    %randomly select background - will load a matrix 'bg'
    bg_idx = ceil(args.num_bgs*rand);
    if exist('bg_list', 'var')
        bg = u_load([args.bg_dir bg_list(bg_idx).name]);
    else   
        bg = u_load([args.bg_dir args.bg_stem zerostr(bg_idx, args.bg_zeros), '.mat']);
    end
    
    %get size of background
    [row col] = size(bg);
    
    %randomly select orientation between 0 and 180 degrees
    orientation = args.orientation_range(1) +...
        (args.orientation_range(2)-args.orientation_range(1))*rand;
    
    %randomly select width between 4 and 32
    width = args.width_range(1) + (args.width_range(2)-args.width_range(1))*rand;
    
    %     %randomly select contrast between 8 and 16
    %     contrast = args.contrast_range(1) + (args.contrast_range(2)-args.contrast_range(1))*rand;
    %randomly select contrast on exponential distribution
    mu = (args.contrast_range(2) - args.contrast_range(1)) / (2*log(args.decay_rate));
    contrast = args.contrast_range(1) + exp_rand_sample(mu);
    
    %randomly select squashiness of the line
    squash = rand;
    
    %randomly select the radial curvature of the line
    radius = 2^(8*rand + 8);
    
    %make  bar
    switch args.bar_type
        case 'ellipse'
            [bar_image, dummy, label] = create_ellipse_bar(width/2, contrast, orientation, row, col, row/2, col/2);
            
        case 'sin'
            [bar_image, dummy, label] = create_sin_bar(width/2, contrast, orientation, row, col, squash, row/2, col/2);
        
        case 'circle'
            [bar_image, dummy, label] = create_sin_curve(width/2, contrast, radius, orientation, squash, row, col, row/2, col/2);
            
        otherwise
            warning(['Bar type: ', args.bar_type, ' not recognised']); %#ok
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
    bar_idx = bar_idx(rand_idx_fg(1:round(length(rand_idx_fg)*fg_ratio)));
    
    rand_idx = randperm(length(bg_idx));
    bg_idx = bg_idx(rand_idx(1:length(bar_idx)*bg_ratio)); clear rand_idx rand_idx_fg;
    
    %Get indices for this image and convert to row,col
    image_idx = [bar_idx; bg_idx];
    num_samples_image = length(image_idx);
    [rows cols] = ind2sub([row col], image_idx);
    %display(['Num samples in image ' num2str(curr_image) ': ' num2str(num_samples_image)]);
    
    % Create DT-CWT of image
    dt = dtwavexfm2b(image_in, num_levels); clear image_in;
    
    %Now sample the training data - if using clock method we split here, all other methods
    %use the else clause
    if ~isempty(strfind(args.feature_type, 'clock'))
        
        %Get samples in clock representation
        [dt_samples] = sample_dt_polar13(dt, rows, cols, strcmpi(args.feature_type, 'clock_rotate')); 
        clear dt rows cols;
        
        %Convert into magnitude and phase
        mags = abs(dt_samples(:,:));
        phases = angle(dt_samples(:,:)); clear dt_samples
        
        %Discard surplus magnitude columns
        mags(:,discard_idx) = [];
        
        %Save into the output data
        training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
            [mags phases];
    else
    
        %Make copies of sample rows and cols at positions of local window patch
        rr = repmat(rows*ones(1,win_size) + ones(num_samples_image,1)*win_idx, 1, win_size);
        cc = kron(cols*ones(1,win_size) + ones(num_samples_image,1)*win_idx, ones(1,win_size));

        %Get interpolated dual-tree coefficients
        dt_samples = dt_to_pixel_subset(dt, rr, cc); clear dt rr cc;

        if args.do_max
            %get the maximum response across orientations
            dt_samples = squeeze(max(dt_samples, [], 3)); clear dt;
        end

        %Reshape data so each row corresponds to a single pixel representation
        temp_samples=reshape(dt_samples, num_samples_image, []);
        clear dt_samples

        %Save the data to the main training data array
        switch args.feature_type
            case 'all'
                training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                    [abs(temp_samples) angle(temp_samples)];

            case 'real'
                training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                    real(temp_samples);

            case 'mag'
                training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                    abs(temp_samples);

            case 'phase'
                training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                    angle(temp_samples);

            case 'rotate'

                %Convert samples into phase and magnitude
                mags = abs(temp_samples);
                phases = angle(temp_samples);

                %Find the sub-band that corresponds to the maximum magnitude in
                %each row
                [dummy max_ori] = max(mags, [], 2);
                max_ori = rem(max_ori-1,6)+1;

                %Circular shift the data according to the maximum subband
                for ori = 1:6
                    shift_idx = max_ori == ori;
                    for lev = 1:num_levels
                        cols = 6*(lev-1)+(1:6);
                        mags(shift_idx, cols) = circshift(mags(shift_idx, cols), [0 1-ori]);
                        phases(shift_idx, cols) = circshift(phases(shift_idx, cols), [0 1-ori]);
                    end
                end

                %Save into the output data
                training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                    [mags phases];

            otherwise
                training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                    [abs(temp_samples) angle(temp_samples)];
        end
    end
    
    training_labels(curr_sample:num_samples_image+curr_sample-1)= label(image_idx);
    clear temp_samples;
    
    if args.plot
        figure; imagesc(label); axis image; colormap(gray(256)); hold on;
        plot(cols(1:length(bg_idx)), rows(1:length(bg_idx)), 'rx');
        plot(cols(length(bg_idx)+1:length(image_idx)), rows(length(bg_idx)+1:length(image_idx)), 'gx');
    end
    
    if curr_sample >= num_samples
        %save parameters
        if ~isempty(parameters_dir)
            if ~(strcmp(parameters_dir(end), '/') || strcmp(parameters_dir(end), '\'));  
                parameters_dir = [parameters_dir '/']; %#ok
            end
            if ~exist(parameters_dir, 'dir')
                mkdir(parameters_dir);
            end
            par_list = dir([parameters_dir '*parameters*']);
            tree_num = length(par_list) + 1;
            save([parameters_dir 'parameters' zerostr(tree_num, 3), '.mat'], 'parameters');
        end

        return;
    else
        %increment sample count
        curr_sample = curr_sample + num_samples_image;
    end
    
    clear bar_idx bg_idx image_idx label;
    %increment image count
    curr_image = curr_image + 1;
end

return
