function [fg_map, bg_map] = ...
    get_real_masks(images, sampling_args, output_type, img_label)
% Compute foreground and background masks given an image list

% Load FoV mask for some prediction types (used to only do this for
% detection output as we assume all the fg pixels are within the FoV so
% outputs that sample only from the fg don't need the FoV - but we'll load
% for all cases now so the FoV mask can be used to mask any region)
if ~isempty(images.fov_mask)   
    % Load in mask of foveal region
    fov_map = u_load(images.fov_mask);
else
   fov_map = [];
end

% Load in line mask
if isfield(images, 'probability')
    %If probability field isn't empty we need to get a resampling map
    
    if sampling_args.make_resampling_maps
        %Will compute resampling probabilities on the fly from
        %predictions which are saved in the probability dir
        prediction_img = u_load(images.probability);
        fg_mask = u_load(images.fg_mask);
        error_ratio = sampling_args.make_resampling_maps;
        [line_map] = make_resampling_map(...
            prediction_img, fg_mask, fov_map, img_label, error_ratio, output_type);
    else
        %Resampling map is already saved in the probability dir
        [line_map] = u_load(images.probability);
    end
else
    %Otherwise the line map is just the logical fg mask
    line_map = u_load(images.fg_mask);
end

% centre_map to be computed from the fg_map
centre_map = [];

% Process and combine maps to create foreground 
[fg_map, bg_map] = process_maps(line_map, centre_map, fov_map, ...
                                output_type, sampling_args.shrink_fov);
