function [training_data training_labels parameters] = sample_training_data_max_pixelset(varargin)
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
    'spicule', 0, ...
    'win_size', 3,...
    'num_levels', 6,...
    'bg_ratio', 1,...
    'fg_ratio', 0.5,...
    'width_range', [2 32],...
    'orientation_range', [0 180],...
    'contrast_range', [4 16],...
    'save_path', [],...
    'plot', 0);

spicule = args.spicule; % if the training data include spicules
win_size = args.win_size; % the size of scanning window
num_levels = args.num_levels; %number of levels of DT-CWT from which to extract coefficients
num_samples = args.num_samples; %total number of feature vectors to sample
bg_ratio = args.bg_ratio; %ration of backgorund pixels to bar pixels
fg_ratio = args.fg_ratio; %ration of sample pixels to the bar pixels

%Constants computed from arguments
pad_w = floor(win_size/2); %half size of window size
size_sample_vector = 2*win_size*win_size*num_levels;
win_idx = -pad_w:pad_w;

%get directory listing of backgrounds
bg_list = dir([args.bg_dir, '*.mat']);

%Create memory storage for data
training_data = zeros(num_samples, size_sample_vector);
training_labels = zeros(num_samples, 1);

curr_sample = 1;
curr_image = 1;

%Create images filling up sample data until training data full
while true %reaching num_samples is checked below, loop terminates with a return statement

    %randomly select background - will load a matrix 'bg'
    bg_idx = ceil(length(bg_list)*rand);
    bg = u_load([args.bg_dir, bg_list(bg_idx).name]);

    %get size of background
    [row col] = size(bg);

    %randomly select orientation between 0 and 180 degrees
    orientation = args.orientation_range(1) +...
        (args.orientation_range(2)-args.orientation_range(1))*rand;

    %randomly select width between 4 and 32
    width = args.width_range(1) + (args.width_range(2)-args.width_range(1))*rand;

    %randomly select contrast between 8 and 16
    contrast = args.contrast_range(1) + (args.contrast_range(2)-args.contrast_range(1))*rand;

    %     randomly select bar type
    if spicule
        if width > 6
            spicule = rand > .5;
        else
            spicule = false;
        end
    end

    if spicule
        %make gaussian bar
        [bar_image, dummy, label] = create_ellipse_spicule(width/2, contrast, orientation, row, col);
    else
        %make rectangular bar
        [bar_image, dummy, label] = create_ellipse_bar(width/2, contrast, orientation, row, col);
    end
    clear dummy;

    %Store bar parameters
    parameters(curr_image).bg_idx = bg_idx; %#ok
    parameters(curr_image).width = width; %#ok
    parameters(curr_image).contrast = contrast; %#ok
    parameters(curr_image).orientation = orientation; %#ok
    parameters(curr_image).spicule = spicule; %#ok

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
    bg_idx = bg_idx(rand_idx(1:length(bar_idx)*bg_ratio));
    clear rand_idx rand_idx_fg;

    %Get indices for this image and convert to row,col
    image_idx = [bar_idx; bg_idx];
    num_samples_image = length(image_idx);
    [rows cols] = ind2sub([row col], image_idx);

    % Create DT-CWT of image
    dt = dtwavexfm2b(image_in, num_levels); clear image_in;

    %Make copies of sample rows and cols at positions of local window patch
    rr = repmat(rows*ones(1,win_size) + ones(num_samples_image,1)*win_idx, 1, win_size);
    cc = kron(cols*ones(1,win_size) + ones(num_samples_image,1)*win_idx, ones(1,win_size));

    %Get interpolated dual-tree coefficients
    dt_samples = dt_to_pixel_subset(dt, rr, cc); clear dt;
    clear rr cc;

    %get the maximum
    dt_samples = squeeze(max(dt_samples,[],3));

    temp_samples=reshape(dt_samples, num_samples_image, []);

    training_data(curr_sample:num_samples_image+curr_sample-1,:) = [abs(temp_samples) angle(temp_samples)];

    training_labels(curr_sample:num_samples_image+curr_sample-1)= label(image_idx);
    clear dt_samples;

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

        %Convert into magnitude and phase, This is a fast way, but out off
        %momery
        %         training_data(:,1+end/2:end) = angle(training_data(:,1:end/2));
        %         training_data(:,1:end/2) = abs(training_data(:,1:end/2));
        return;
    else
        %increment sample count
        curr_sample = curr_sample + num_samples_image;
    end

    clear  bar_idx bg_idx image_idx label cols rows;
    %increment image count
    curr_image = curr_image + 1;
end

return
