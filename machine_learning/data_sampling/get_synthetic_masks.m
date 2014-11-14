function [fg_map, bg_map] = ...
    get_synthetic_masks(line_map, centre_map, ...
                        sampling_args, output_type)

% Create a FoV map that masks out k>=1 edge pixels
k = min(floor(sampling_args.win_size / 2), 1);
fov_map = true(size(line_map));
fov_map([1:k end-k+1:end], :) = false;
fov_map(:, [1:k end-k+1:end]) = false;

[fg_map, bg_map] = process_maps(line_map, centre_map, fov_map, ...
                                output_type, sampling_args.shrink_fov);
