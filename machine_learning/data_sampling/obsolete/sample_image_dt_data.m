function [training_data] = sample_image_dt_data(varargin)
%SAMPLE_MAMMO_TRAINING_DATA *Insert a one line summary here*
%   [] = sample_mammo_training_data(varargin)
%
% SAMPLE_MAMMO_TRAINING_DATA uses the U_PACKARGS interface function
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

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    {'num_samples',... % the mandatory arguments
    'image_dir'}, ...
    'win_size', 3,...
    'num_levels', 4,...
    'feature_type', 'all',...
    'feature_shape', 'rect',...
    'rotate', 0,...
    'do_max', 0,...
    'image_list', [], ...
    'image_type', '.mat',...
    'use_nag', 1,...
    'save_path', []);
clear varargin;

win_size = args.win_size; % the size of scanning window
num_samples = args.num_samples; %total number of feature vectors to sample
image_dir = args.image_dir; %directory in which images are stored
image_list = args.image_list;

sampling_args.feature_shape = args.feature_shape;
sampling_args.feature_type = args.feature_type;
sampling_args.do_max = args.do_max;
sampling_args.rotate = args.rotate;
sampling_args.win_size = args.win_size;

size_sample_vector = args.num_levels*compute_dt_feature_size(sampling_args);

if isempty(image_list)
    image_list = dir([image_dir, '*', args.image_type]);
end

%Constants computed from arguments
pad_w = floor(win_size/2); %half size of window size

%Get list of images in image dir
%Create data
training_data = zeros(num_samples, size_sample_vector);

curr_sample = 1;
num_images = length(image_list);
im_order = randperm(num_images);

%loop through each image sampling data
for kk = 1:num_images
    
    %num_samples_image = binornd((num_samples + 1 - curr_sample), 1/(num_images+1-kk), 1, 1);
    num_samples_image = ...
        sample_from_binomial((num_samples + 1 - curr_sample), 1/(num_images+1-kk), 1);
    
    %display(['Samples in image ' num2str(kk) ': ', num2str(num_samples_image)]);
    
    if ~num_samples_image
        continue;
    end
    
    %load in image
    image_in = u_load([image_dir image_list(im_order(kk)).name]);
    [row, col] = size(image_in);
    num_samples_image = min(num_samples_image, numel(image_in)); %just in case...

    
    %Mask off edge of image
    if pad_w <1
        pad_w=1;
    end
    mask = true(row, col);
    mask([1:pad_w end-pad_w+1], :) = false;
    mask(:, [1:pad_w end-pad_w+1]) = false;
    
    %Get indices of all pixels in the mask for this image
    image_idx = find(mask);
    clear mask;
    
    %Now select a random subset of the pixels in the mask and convert to row,col
    rp = randperm(length(image_idx));
    image_idx = image_idx(rp(1:num_samples_image)); clear rp;
    [rows cols] = ind2sub([row col], image_idx);
    
    %Compute dual-tree transform of image
    dt = compute_dual_tree(image_in, args.num_levels, args.use_nag);

    %Sample DT coefficients from specified rows and cols according to
    %sampling arguments
    training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
        sample_dt_data(dt, rows, cols, sampling_args);
    clear dt;

    curr_sample = curr_sample + num_samples_image;
end