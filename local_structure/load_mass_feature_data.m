function [training_data] = load_mass_feature_data(varargin)
%LOAD_MASS_FEATURE_DATA *Insert a one line summary here*
%   [] = sample_mammo_training_data(varargin)
%
% LOAD_MASS_FEATURE_DATA uses the U_PACKARGS interface function
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
    'mammo_idx', ...
    'samples_path'}, ...
    'dist_range', [16 32 64 128, 256],...
    'sigma_range', [1 2 4 8],...
    'angular_res', 1,...
    'quiet', 1,...
    'do_template', 1,...
    'do_scale', 1,...
    'save_path', []);
clear varargin;

if args.quiet
    warning('off', 'load_uint8:missing_variables');
end

%Get constants from args
num_images = length(args.mammo_idx);
num_dists = length(args.dist_range);
num_sigma = length(args.sigma_range);
size_sample_vector = num_dists*num_sigma*args.angular_res + args.do_template + args.do_scale;

%Load in the sample per image data
samples_per_image = u_load([args.samples_path '_samples_per_image.mat']);

%Workout index offsets for each image
idx_offset = [0; cumsum(samples_per_image(1:end-1))];

%Get a random ordering of the images
im_order = randperm(num_images);

%Pre-allocate indices
sample_idx = zeros(args.num_samples, 1);
curr_sample = 1;

%First we work out a set of sampling indices from all the images
for ii = 1:num_images
    
    mm = args.mammo_idx(im_order(ii));
    
    %Check we've got samples to take from this image
    if ~samples_per_image(mm)
        continue;
    end

    %Randomly sample how many points to coose from this image
    num_samples_image = ...
        sample_from_binomial((args.num_samples + 1 - curr_sample), 1/(num_images+1-ii), 1);
    
    %check this isn't more than we can manage
    num_samples_image = min(num_samples_image, samples_per_image(mm));
    
    %Create random permutation to select random samples
    rp = randperm(samples_per_image(mm));
    
    %Save random samples idx in main data structure
    sample_idx(curr_sample:num_samples_image+curr_sample-1) =...
        idx_offset(mm) + rp(1:num_samples_image);
    
    %Increment current sample counter
    curr_sample = curr_sample + num_samples_image;
end

%Discard any samples we have not filled
sample_idx(curr_sample:end) = [];

%Preallocate storage for sampled data
%training_data = zeros(curr_sample-1, size_sample_vector);
training_data = zeros(curr_sample-1, size_sample_vector);

%Now we have sample indices, we can loop through each dimension of the
%pre-sampled data to populate the training data
for ii = 1:num_dists
    for jj = 1:num_sigma
        
        %Load in pre-sampled data
        dimension_data = load_uint8([args.samples_path ...
            '_d' zerostr(args.dist_range(ii), 3) ...
            '_s' zerostr(args.sigma_range(jj),2) '.mat']);

        %Workout columns to be populate by this data
        col_offset = (num_sigma*args.angular_res)*(ii-1) + args.angular_res*(jj-1);
        cc = (1:args.angular_res) + col_offset;
        
        %Sample data at required indices
        training_data(:,cc) = dimension_data(sample_idx, :);
        clear dimension_data;
    end
end

%Finally we need to load the template and scale map data
if args.do_template
    template_data = load_uint8([args.samples_path '_template_data.mat']);
    training_data(1:curr_sample-1,num_dists*num_sigma*args.angular_res + 1) = template_data(sample_idx, :);
    clear template_data;
end
if args.do_scale
    scale_data = load_uint8([args.samples_path '_scale_data.mat']);
    training_data(:,end) = scale_data(sample_idx, :);
    clear scale_data;
end
            