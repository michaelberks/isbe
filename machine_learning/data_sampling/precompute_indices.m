function [fg_indices, bg_indices] = precompute_indices(args)
%PRECOMPUTE_INDICES *Insert a one line summary here*
%   [fg_indices, bg_indices] = precompute_indices(args)
%
% Inputs:
%      args - *Insert description of input variable here*
%
%
% Outputs:
%      fg_indices - *Insert description of input variable here*
%
%      bg_indices - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 02-Dec-2015
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
switch args.precompute_method
    case 'sample_to_uniformity'
        [fg_indices, bg_indices] = sample_image_to_uniformity(args.precompute_args);

    otherwise
        error(['Index pre-computation method ' args.precompute_method ' not recognised']);
end


function [fg_indices, bg_indices] = sample_image_to_uniformity(args)

%sample uniformly
r_idx = sample_to_uniformity(args.values, args.num_samples, args.start_value, args.end_value);

num_images = length(args.all_fg_indices);
fg_indices = cell(num_images,1);
bg_indices = cell(num_images,1);
start_sample = 0;

for i_im = 1:num_images
    end_sample = start_sample + length(args.all_fg_indices{i_im});

    r_idx_i = r_idx(r_idx > start_sample & r_idx <= end_sample) - start_sample;
    start_sample = end_sample;

    fg_indices{i_im} = args.all_fg_indices{i_im}(r_idx_i);
    bg_indices{i_im} = [];
end