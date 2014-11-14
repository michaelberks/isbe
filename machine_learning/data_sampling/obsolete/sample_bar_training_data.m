function [training_data training_labels parameters] = sample_bar_training_data(varargin)
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
             'num_levels', 4,...
             'bg_dir', 'C:\isbe\dev\classification\data\normal_smooth128\',...
             'bg_ratio', 4,...
             'spicule_ratio', 0.5,...
             'width_range', [2 16],...
             'orientation_range', [0 180],...
             'contrast_range', [8 16],...
             'save_path', [],...
             'plot', 0);

win_size = args.win_size; % the size of scanning window
num_levels = args.num_levels; %number of levels of DT-CWT from which to extract coefficients
num_samples = args.num_samples; %total number of feature vectors to sample
bg_ratio = args.bg_ratio; %ration of backgorund pixels to bar pixels

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

    %randomly select bar type
    spicule = rand < args.spicule_ratio;

    if spicule    
        %make gaussian bar
        [bar_image, dummy, label] = create_ellipse_spicule_mb(width/2, contrast, orientation, row, col);
    else       
        %make rectangular bar
        [bar_image, dummy, label] = create_ellipse_bar_mb(width/2, contrast, orientation, row, col);
    end
    clear dummy;
    
    %Store bar parameters
    parameters(curr_image).bg_idx = bg_idx; %#ok
    parameters(curr_image).width = width; %#ok
    parameters(curr_image).contrast = contrast; %#ok
    parameters(curr_image).orientation = orientation; %#ok
    parameters(curr_image).spicule = spicule; %#ok
    
    %Add bar to background
    image_in = bar_image + bg;
    
    %Use label to get indices for this image
    bar_idx = find(label); %use all bar pixels
    bg_idx = find(~label);
    
    %Now select a random selection of the background pixels
    rand_idx = randperm(length(bg_idx));
    bg_idx = bg_idx(rand_idx(1:length(bar_idx)*bg_ratio));   
    
    % Create DT-CWT of image
    %dt = dtwavexfm2(image_in,num_levels,'near_sym_b','qshift_b');
    dt = dtwavexfm2b(image_in, num_levels);

    %interpolate DT-CWT to full pixel grid
    dt_inter = dt_to_full_image(dt); clear dt;

    %get the maximum response across orientations
    max_dt = squeeze(max(dt_inter, [], 3)); clear dt;

    %pad the edges of subbands and label
    pad_dt = padarray(max_dt, [pad_w pad_w], 'replicate'); clear max_dt;
    pad_label = padarray(label, [pad_w pad_w], 'replicate');

    %Get indices for this image and convert to row,col
    image_idx = [bar_idx; bg_idx];
    num_samples_image = length(image_idx);
    [rows cols] = ind2sub([row col], image_idx);

    if args.plot
        figure; imagesc(pad_label); axis image; colormap(gray(256)); hold on;
    end
    
    %loop through image pixels selecting points
    for ii = 1:num_samples_image

        rr = rows(ii) + pad_w;
        cc = cols(ii) + pad_w;

        %Extract label and sample vector
        sample_label = pad_label(rr, cc);
        sample_vector = pad_dt(rr+win_idx, cc+win_idx, :);

        %Convert sample vector to mag/phase form and save in main data
        %training_data(curr_sample,:) = [abs(sample_vector(:))' angle(sample_vector(:))'];
        training_data(curr_sample,1:end/2) = sample_vector(:).';
        training_labels(curr_sample,:)= sample_label;
        
        if args.plot
            if sample_label
                plot(cc, rr, 'rx');
            else
                plot(cc, rr, 'gx');
            end
        end
        
        if curr_sample == num_samples
            %save parameters
            if ~isempty(args.save_path)
                save([args.save_path, 'bar_parameters.mat'], 'parameters');
            end
            
            %Convert into magnitude and phase
            training_data(:,1+end/2:end) = angle(training_data(:,1:end/2));
            training_data(:,1:end/2) = abs(training_data(:,1:end/2));
            return;
        else
            %increment sample count
            curr_sample = curr_sample + 1;
        end
    end
    clear pad_dt;
    %increment image count
    curr_image = curr_image + 1;
end
        