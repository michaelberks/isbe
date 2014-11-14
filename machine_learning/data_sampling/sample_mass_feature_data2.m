function [training_data] = sample_mass_feature_data2(varargin)
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
    'save_path', []);
clear varargin;

num_samples = args.num_samples; %total number of feature vectors to sample
mask_dir = args.mask_dir; %directory in which masks are stored
radial_dir = args.radial_dir; %directory in which radial maps are stored
template_dir = args.template_dir; %directory in which mass maps are stored
mammo_list = args.mammo_list;
num_dist = length(args.dist_range);
num_sigma = length(args.sigma_range);
size_sample_vector = num_dist*num_sigma*args.angular_res + 2;

if isempty(mammo_list)
    mammo_list = dir([mask_dir '*.mat']);
end
[mammo_names] = get_mammo_info(mammo_list);

%Get list of images in image dir
%Create data
training_data = zeros(num_samples, size_sample_vector);

curr_sample = 1;
num_images = length(mammo_list);
im_order = randperm(num_images);

%loop through each image sampling data
for mm = 1:num_images
    
    %num_samples_image = binornd((num_samples + 1 - curr_sample), 1/(num_images+1-mm), 1, 1);
    num_samples_image = ...
        sample_from_binomial((num_samples + 1 - curr_sample), 1/(num_images+1-mm), 1);
    display(['Samples in image ' num2str(mm) ': ', num2str(num_samples_image)]);
    
    if ~num_samples_image
        continue;
    end
    
    %load in mask
    mask_name = dir([mask_dir, '*' mammo_names{im_order(mm)} '*.mat']);
    mask = u_load([mask_dir mask_name(1).name]);
    
    num_samples_image = min(num_samples_image, sum(mask(:))); %just in case...
    
    %Get indices of all pixels in the mask for this image
    image_idx = find(mask);
    
    
    %Now select a random subset of the pixels in the mask and convert to row,col
    rp = randperm(length(image_idx));
    image_idx = image_idx(rp(1:num_samples_image)); clear rp;
    
    %Now populate the training data with samples from the radial map
    for ii = 1:num_dist
        
        %load the radial map for this distance
        rad_name = dir([radial_dir, '*' mammo_names{im_order(mm)} '*' zerostr(args.dist_range(ii),3) '*.mat']);
        
        if isempty(rad_name)
            %try loading it again - SAN may just be playing up
            retry = 0;
            while isempty(rad_name) && retry < 10
                rad_name = dir([radial_dir, '*' mammo_names{im_order(mm)} '*' zerostr(args.dist_range(ii),3) '*.mat']);
                retry = retry + 1;
            end
            if isempty(rad_name)
                display(['No radial map found: ' radial_dir, '*' mammo_names{im_order(mm)} '*' zerostr(args.dist_range(ii),3) '*.mat']);
                error('Radial map not found');
            end
        end
        try
            rad_map = load_uint8([radial_dir rad_name(1).name]);
        catch
            rad_map = u_load([radial_dir rad_name(1).name]);
        end
        
        for jj = 1:num_sigma

            [rows cols num_bands] = size(rad_map);
            if any([rows cols] ~= size(mask))
                error(['Error in ' rad_name(1).name ': Size of mask and radial map do not match']);
            elseif num_bands ~= args.angular_res
                error(['Error in ' rad_name(1).name...
                    ': Number of bands in radial map (' num2str(num_bands)...
                    ') does not match angular resolution (' num2str(args.angular_res) ')']);
            end
            
            for kk = 1:num_bands
                %smooth the radial map by current sigma
                angle_band = imfilter(rad_map(:,:,kk),...
                    fspecial('gaussian', 5*args.sigma_range(jj), args.sigma_range(jj)), 'symmetric');
            
                %Compute column and offset of indices for this band
                cc = (num_sigma*num_bands)*(ii-1) + num_bands*(jj-1) + kk; 
            
                %Extract samples into main data structure
                training_data(curr_sample:num_samples_image+curr_sample-1,cc) = ...
                    angle_band(image_idx);
            end
            
%             %Finally save the sum of the angle bands
%             training_data(curr_sample:num_samples_image+curr_sample-1, cc+1) = ...
%                 sum(training_data(...
%                     curr_sample:num_samples_image+curr_sample-1,...
%                     (num_sigma*num_bands)*(ii-1) + num_bands*(jj-1) + 1:cc), 2);
            
        end
        clear rad_map;
    end
    clear mask;
    
    %Now populate the training data with mass template score and scale
    mass_name = dir([template_dir '*' mammo_names{im_order(mm)} '*.mat']);
    mass_map = load_uint8([template_dir mass_name(1).name]);
    training_data(curr_sample:num_samples_image+curr_sample-1,end-1) = ...
        mass_map(image_idx);
    clear mass_map;
    
    scale_name = dir([template_dir, 'scales/*' mammo_names{im_order(mm)} '*.mat']);
    scale_map = load_uint8([template_dir 'scales/' scale_name(1).name]);
    training_data(curr_sample:num_samples_image+curr_sample-1,end) = ...
        scale_map(image_idx);
    clear scale_map;

    %Increment current sample counter
    curr_sample = curr_sample + num_samples_image;
end