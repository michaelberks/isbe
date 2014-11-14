function save_radial_data(varargin)
%SAVE_RADIAL_DATA *Insert a one line summary here*
%   [] = sample_mammo_training_data(varargin)
%
% SAVE_RADIAL_DATA uses the U_PACKARGS interface function
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
    {'save_path',... % the mandatory arguments
    'mask_dir',...
    'radial_dir',...
    'template_dir'}, ...
    'samples_per_image', [],...
    'dist_range', [16 32 64 128, 256],...
    'sigma_range', [1 2 4 8],...
    'angular_res', 1,...
    'do_template', 1,...
    'do_scale', 1,...
    'mammo_list', [], ...
    'quiet', 1);
clear varargin;

if args.quiet
    warning('off', 'load_uint8:missing_variables');
end

%Get values from args
samples_per_image = args.samples_per_image; %total number of feature vectors to sample
mask_dir = args.mask_dir; %directory in which masks are stored
mammo_list = args.mammo_list;
num_dists = length(args.dist_range);
num_sigma = length(args.sigma_range);

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

%loop through each mask to work out total number of samples
if isempty(samples_per_image)
    samples_per_image = zeros(num_images,1);
    for ii = 1:num_images
        mm = valid_idx(ii);
        
        %load in mask
        mask = u_load([mask_dir mask_names{mm}]);
        samples_per_image(ii) = sum(mask(:));
    end
end
save([args.save_path '_samples_per_image.mat'], 'samples_per_image');
num_samples = sum(samples_per_image);
            
%Sample and save training data for each dimension of samples vector
for ii = 1:num_dists
    for jj = 1:num_sigma
        
        %Preallocate storage for sampled data
        training_data = zeros(num_samples, args.angular_res);
        curr_sample = 1;

        %Now loop through images
        for kk = 1:num_images

            mm = valid_idx(kk);
    
            %Check we've got samples to take from this image
            if ~samples_per_image(kk)
                continue;
            end
    
            %load in mask
            mask = u_load([mask_dir mask_names{mm}]);

            %Load in radial map
            rad_map = load_uint8([args.radial_dir radial_names{mm,ii}]);
            
            num_bands = size(rad_map, 3);
            if num_bands ~= args.angular_res
                error(['Error in ' rad_name(1).name...
                    ': Number of bands in radial map (' num2str(num_bands)...
                    ') does not match angular resolution (' num2str(args.angular_res) ')']);
            end

            for ll = 1:num_bands
                %smooth the radial map by current sigma
                angle_band = imfilter(rad_map(:,:,ll),...
                    fspecial('gaussian', 5*args.sigma_range(jj), args.sigma_range(jj)), 'symmetric');

                %Extract samples into main data structure
               training_data(curr_sample:samples_per_image(kk)+curr_sample-1, ll) = ...
                   angle_band(mask);
            end
            
            %Increment current sample counter
            curr_sample = curr_sample + samples_per_image(kk);
        end
        
        %Save the data
        data_name = [args.save_path '_d' zerostr(args.dist_range(ii), 3) '_s' zerostr(args.sigma_range(jj),2) '.mat'];
        %save_uint8(data_name, training_data);
        save(data_name, 'training_data');
        clear training_data
    end
end

%Now do the same for the template and scale data

template_data = zeros(num_samples, 1);
scale_data = zeros(num_samples, 1);
curr_sample = 1;

%Loop through images
for kk = 1:num_images

    mm = valid_idx(kk);

    %Check we've got samples to take from this image
    if ~samples_per_image(kk)
        continue;
    end

    %load in mask
    mask = u_load([mask_dir mask_names{mm}]);

    %Load in template map
    template_map = load_uint8([args.template_dir template_names{mm}]);

    %Extract samples into main data structure
    template_data(curr_sample:samples_per_image(kk)+curr_sample-1, 1) = ...
        template_map(mask);
    
    %Load in template map
    scale_map = load_uint8([args.template_dir 'scales/' scale_names{mm}]);

    %Extract samples into main data structure
    scale_data(curr_sample:samples_per_image(kk)+curr_sample-1, 1) = ...
        scale_map(mask);

    %Increment current sample counter
    curr_sample = curr_sample + samples_per_image(kk);
end

%Save the data
% save_uint8([args.save_path '_template_data.mat'], template_data);
% save_uint8([args.save_path '_scale_data.mat'], scale_data);
save([args.save_path '_template_data.mat'], 'template_data');
save([args.save_path '_scale_data.mat'], 'scale_data');
clear template_data scale_data
