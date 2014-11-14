function [training_data training_oris parameters] = sample_orientation_training_data(varargin)
%SAMPLE_BAR_TRAINING_DATA *Insert a one line summary here*
%   [] = sample_bar_training_data(varargin)
%
% SAMPLE_BAR_TRAINING_DATA uses the U_PACKARGS interface function
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
			 {'num_samples'},... % the mandatory arguments
             'win_size', 3,...
             'num_levels', 5,...
             'do_max', 0,...
             'width_range', [4 16],...
             'orientation_range', [0 180],...
             'contrast_range', [4 16],...
             'line_type', 'ellipse',...
             'feature_type', 'all',...
             'save_path', [],...
             'idx_ratio', 1,...
             'plot', 0);

win_size = args.win_size; % the size of scanning window
num_levels = args.num_levels; %number of levels of DT-CWT from which to extract coefficients
num_samples = args.num_samples; %total number of feature vectors to sample

%Constants computed from arguments
pad_w = floor(win_size/2); %half size of window size

%Compute size of training vector
if args.do_max
    size_sample_vector = 2*win_size*win_size*num_levels;
else
    size_sample_vector = 2*win_size*win_size*num_levels*6;
end

%Set local window offsets
win_idx = -pad_w:pad_w;

% %get directory listing of backgrounds
% bg_list = dir([args.bg_dir, '*.mat']);

%Create memory storage for data
training_data = zeros(num_samples, size_sample_vector);
training_oris = zeros(num_samples, 1);

curr_sample = 1;
curr_image = 1;

%Create images filling up sample data until training data full
while true %reaching num_samples is checked below, loop terminates with a return statement        

    %randomly select background - will load a matrix 'bg'
    bg = zeros(128);

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
    decay_rate = 4;
    mu = (args.contrast_range(2) - args.contrast_range(1)) / (2*log(decay_rate));
    contrast = args.contrast_range(1) + exp_rand_sample(mu);  
    
    %randomly select squashiness of the line
    squash = rand;
    
    %make  bar
    switch args.line_type
        case 'ellipse'
            [bar_image, dummy, label] = create_ellipse_bar(width/2, contrast, orientation, row, col, row/2, col/2);
            
        case 'sin'
            [bar_image, dummy, label] = create_sin_bar(width/2, contrast, orientation, row, col, squash, row/2, col/2);
            
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
    image_in = bar_image + bg;
    
    %Use label to get indices for this image
    bar_idx = find(label); %use all bar pixels
    rand_idx = randperm(length(bar_idx));
    
    %Get indices for this image and convert to row,col
    num_samples_image = round(args.idx_ratio * rand * length(bar_idx));
    
    if num_samples_image
        image_idx = bar_idx(rand_idx(1:num_samples_image));
    else
        continue;
    end
    
    [rows cols] = ind2sub([row col], image_idx);

    % Create DT-CWT of image
    dt = dtwavexfm2b(image_in, num_levels); clear image_in;

    %Make copies of sample rows and cols at positions of local window patch
    rr = repmat(rows*ones(1,win_size) + ones(num_samples_image,1)*win_idx, 1, win_size);
    cc = kron(cols*ones(1,win_size) + ones(num_samples_image,1)*win_idx, ones(1,win_size));

    %Get interpolated dual-tree coefficients
    dt_samples = dt_to_pixel_subset(dt, rr, cc); clear dt rr cc;
    dt_samples = reshape(dt_samples, num_samples_image, []);
    
    
    training_data(curr_sample:num_samples_image+curr_sample-1,:) = [abs(dt_samples) angle(dt_samples)];
    %training_oris(curr_sample:num_samples_image+curr_sample-1) = exp(i*pi*orientation/180);
    training_oris(curr_sample:num_samples_image+curr_sample-1) = orientation;
    
    clear dt_samples temp_samples;

    if args.plot
        figure; imagesc(label); axis image; colormap(gray(256)); hold on;
        plot(cols(1:length(bg_idx)), rows(1:length(bg_idx)), 'rx');
        plot(cols(length(bg_idx)+1:length(image_idx)), rows(length(bg_idx)+1:length(image_idx)), 'gx');
    end

    if curr_sample >= num_samples
        %save parameters
        if ~isempty(args.save_path)
            save(args.save_path, 'parameters');
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
        