function [training_data] = sample_spicule_training_data_g2d(varargin)
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
			 {'image_list',...
             'image_sample_pts',...
             'image_dir'},... % the mandatory arguments
             'win_size', 3,...
             'sigma_range', [1 2 4 8],...
             'save_path', [],...
             'plot', 0);

sampling_args.win_size = args.win_size;
sigma_range = args.sigma_range;
    
num_images = length(args.image_sample_pts);
if length(args.image_list) ~= num_images
    error('Length of image list differs from length of image sample pts');
end

%Constants computed from arguments
pad_w = floor(args.win_size/2); %half size of window size

%Compute size of training vector
size_sample_vector = 3*args.win_size*args.win_size*length(sigma_range);

%Workout number of points to sample
num_samples = 0;
for ii = 1:num_images
    num_samples = num_samples + size(args.image_sample_pts{ii},1);
end

%Create memory storage for data
training_data = zeros(num_samples, size_sample_vector);

curr_sample = 1;

%Go through each image sampling points at specified indices        
for ii = 1:num_images
    %randomly select background - will load a matrix 'bg'
    image_in = u_load([args.image_dir, args.image_list(ii).name]);
    image_in = padarray(image_in, [pad_w pad_w], 'replicate');
    
    %get sample pts for this image
    num_samples_image = size(args.image_sample_pts{ii},1);
    if num_samples_image
        rows = args.image_sample_pts{ii}(:,2); %y-coordinates of sample pts specify rows
        cols = args.image_sample_pts{ii}(:,1); %x-coordinates of sample pts specify columns
    else
        continue;
    end

    display(['Sampling points from image: ' num2str(ii)]);
    
    %Compute gaussian 2nd derivative responses 
    g2d_responses = compute_gaussian_2nd_derivatives(image_in,  sigma_range);
    training_data(curr_sample:num_samples_image+curr_sample-1,:) = sample_g2d_data(...
        g2d_responses(:,:,:,1),... 
        g2d_responses(:,:,:,2),...
        g2d_responses(:,:,:,3), rows+pad_w, cols+pad_w, sampling_args);
 
    %increment sample count
    curr_sample = curr_sample + num_samples_image;

    if args.plot
        figure; imagesc(image_in); axis image; colormap(gray(256)); hold on;
        plot(cols+pad_w, rows+pad_w, 'rx');
    end
    clear image_in;
end
        