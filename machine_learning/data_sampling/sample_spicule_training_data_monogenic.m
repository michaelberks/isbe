function [training_data] = sample_spicule_training_data_monogenic(varargin)
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
             'num_levels', 6,...
             'min_wavelength', 4,...
             'onf', 0.65,...
             'save_path', [],...
             'plot', 0);

win_size = args.win_size; % the size of scanning window
num_levels = args.num_levels; %number of levels of DT-CWT from which to extract coefficients
min_wavelength = args.min_wavelength;
onf = args.onf;
    
num_images = length(args.image_sample_pts);
if length(args.image_list) ~= num_images
    error('Length of image list differs from length of image sample pts');
end

%Constants computed from arguments
pad_w = floor(win_size/2); %half size of window size

%Compute size of training vector
size_sample_vector = 2*win_size*win_size*num_levels;

%Set local window offsets
win_idx = -pad_w:pad_w;

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

    %Make copies of sample rows and cols at positions of local window patch
    rr = repmat(rows*ones(1,win_size) + ones(num_samples_image,1)*win_idx, 1, win_size);
    cc = kron(cols*ones(1,win_size) + ones(num_samples_image,1)*win_idx, ones(1,win_size));
    patch_idx = sub2ind(size(image_in), rr+pad_w, cc+pad_w); clear rr cc;
    
    
    %Now get monogenic samples from the image
    [local_amp local_phase] = monogenic(image_in, num_levels, min_wavelength, 2, onf);
    
    %save sample into main training data
    local_amp = reshape(local_amp(:,:,2:end), [], num_levels);
    local_phase = reshape(local_phase(:,:,2:end), [], num_levels);
    
    for level = 1:num_levels
        amp_cols = (1:win_size^2) + (level-1)*win_size^2;
        training_data(curr_sample:num_samples_image+curr_sample-1, amp_cols) = ...
            reshape(local_amp(patch_idx,level), num_samples_image, win_size^2);
        
        phase_cols = (1:win_size^2) + (num_levels + level-1)*win_size^2;
        training_data(curr_sample:num_samples_image+curr_sample-1, phase_cols) = ...
            reshape(local_phase(patch_idx,level),num_samples_image, win_size^2);
    end
    
    %increment sample count
    curr_sample = curr_sample + num_samples_image;

    if args.plot
        figure; imagesc(image_in); axis image; colormap(gray(256)); hold on;
        plot(cols, rows, 'rx');
    end
    clear image_in;
end
        