function [training_data] = sample_mass_feature_data(varargin)
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
    'mask_dir',...
    'radial_dir',...
    'template_dir'}, ...
    'dist_range', [16 32 64 128, 256],...
    'sigma_range', [1 2 4 8],...
    'angular_res', 1,...
    'mammo_list', [], ...
    'quiet', 1,...
    'save_path', []);
clear varargin;

if args.quiet
    warning('off', 'load_uint8:missing_variables');
end

%Get constants from args
num_samples = args.num_samples; %total number of feature vectors to sample
mask_dir = args.mask_dir; %directory in which masks are stored
mammo_list = args.mammo_list;
num_dists = length(args.dist_range);
num_sigma = length(args.sigma_range);
size_sample_vector = num_dists*num_sigma*args.angular_res + 2;

%Copy args into sampling args structure
sampling_args.radial_dir = args.radial_dir; %directory in which radial maps are stored
sampling_args.template_dir = args.template_dir; %directory in which mass maps are stored
sampling_args.sigma_range = args.sigma_range;
sampling_args.angular_res = args.angular_res;

%If not supplied, get list of maks
if isempty(mammo_list)
    mammo_list = dir([mask_dir '*.mat']);
end

%Get mammo names
[mammo_names] = get_mammo_info(mammo_list);
num_images = length(mammo_list);
clear mammo_list;

%Get list of masks
[mask_names missing_idx] =...
    match_mammo_names(args.mask_dir, mammo_names);

%Get list of template maps
[template_names template_idx] =...
    match_mammo_names(args.template_dir, mammo_names);
missing_idx = union(missing_idx, template_idx);

%Get list of template map scales
[scale_names scale_idx] =...
    match_mammo_names([args.template_dir '/scales'], mammo_names);
missing_idx = union(missing_idx, scale_idx);

%Get lists of radial maps
radial_names = cell(num_images, num_dists);
for ii = 1:num_dists
    [radial_names(:,ii) missing_rad] = ...
        match_mammo_names(args.radial_dir, mammo_names, zerostr(args.dist_range(ii),3));
    missing_idx = union(missing_idx, missing_rad);
end

%Now work out valid images
if isempty(missing_idx)
    valid_idx = 1:num_images;
else
    warning('mass_feature:incomplete','There was incomplete data for this set of mammograms');
    valid_idx = setdiff(1:num_images, missing_idx)';
    num_images = length(valid_idx);
end

%Preallocate storage for sampled data
training_data = zeros(num_samples, size_sample_vector);
curr_sample = 1;

%Get a random ordering of the images
im_order = randperm(num_images);

%loop through each image sampling data
for ii = 1:num_images
    
    mm = valid_idx(im_order(ii));
    
    %num_samples_image = binornd((num_samples + 1 - curr_sample), 1/(num_images+1-mm), 1, 1);
    num_samples_image = ...
        sample_from_binomial((num_samples + 1 - curr_sample), 1/(num_images+1-ii), 1);
    
    display(['Samples in image ' num2str(ii) ': ', num2str(num_samples_image)]);
    
    if ~num_samples_image
        continue;
    end
    
    %load in mask
    mask = u_load([mask_dir mask_names{mm}]);
    
    num_samples_image = min(num_samples_image, sum(mask(:))); %just in case...
    
    %Get indices of all pixels in the mask for this image
    image_idx = find(mask);
    
    %Now select a random subset of the pixels in the mask and convert to row,col
    rp = randperm(length(image_idx));
    image_idx = image_idx(rp(1:num_samples_image)); clear rp;
    
    %Now populate the training data with samples from the radial and
    %template maps
    sampling_args.radial_names = radial_names(mm,:);
    sampling_args.template_name = template_names{mm};
    sampling_args.scale_name = scale_names{mm};
    training_data(curr_sample:num_samples_image+curr_sample-1, :) = ...
        sample_radial_data(image_idx, sampling_args);

    %Increment current sample counter
    curr_sample = curr_sample + num_samples_image;
end