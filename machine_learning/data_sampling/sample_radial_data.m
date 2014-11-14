function [radial_samples] = sample_radial_data(image_idx, varargin)
%SAMPLE_RADIAL_DATA *Insert a one line summary here*
%   [dt_samples] = sample_radial_data(dt, rows, cols, rotate)
%
% Inputs:
%
%
% Outputs:
%      radial_samples - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 23-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

% Unpack the arguments:
args = u_packargs(varargin, '0',...
    {'radial_dir',...
     'template_dir',...
     'radial_names'},...
    'template_name', [],...
    'scale_name', [],...
    'sigma_range', [1 2 4 8],...
    'angular_res', 1,...
    'do_template', 1,...
    'do_scale', 1);

num_dist = length(args.radial_names);
num_sigma = length(args.sigma_range);
size_sample_vector = num_dist*num_sigma*args.angular_res + args.do_template + args.do_scale;

radial_samples = zeros(length(image_idx), size_sample_vector);

%Now populate the training data with samples from the radial map
for ii = 1:num_dist

    %load the radial map for this distance
    rad_map = load_uint8([args.radial_dir args.radial_names{ii}]);

    %No get readings for each sigma smoothed copy of the map
    for jj = 1:num_sigma

        num_bands = size(rad_map, 3);
        if num_bands ~= args.angular_res
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
            radial_samples(:,cc) = angle_band(image_idx);
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

if args.do_template
    %Now populate the training data with mass template score and scale
    %load the template map
    template_map = u_load([args.template_dir args.template_name]);
    radial_samples(:,cc+1) = template_map(image_idx);
    clear template_map;
end
if args.do_scale
    %load the scale map
    scale_map = u_load([args.template_dir 'scales/' args.scale_name]);
    radial_samples(:,end) = scale_map(image_idx);
    clear scale_map;
end