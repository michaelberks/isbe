function [img, img_label, fg_map, bg_map, parameters, fg_indices, bg_indices] = ...
    get_next_image(sampling_args, varargin)
% Given a sampling method, produce the next image in the sequence, its
% ground truth labelling, and the foreground/background maps.

% For images loaded from file, the index of the filename is stored as a
% persistent variable (Matlab's equivalent of a static).

% For images generated on the fly, a synthetic image is created.

f_debug = (nargin == 0 && nargout == 0);
if f_debug, test_script(); return; end

[img, img_label, fg_map, bg_map, parameters] = ...
    func(sampling_args, varargin{:});


%% The function
function [img, img_label, fg_map, bg_map, parameters] = ...
    func(sampling_args, image_lists)

persistent img_index;

% If the argument is the string 'reset' then go back to the beginning of
% the list and return.
if ischar(sampling_args) && strcmp(sampling_args,'reset')
    img_index = [];
    img = []; img_label = [];
    fg_map = []; bg_map = [];
    parameters = [];
    return;
end

% Persistent variables are initialized to [] so set it to zero manually
if isempty(img_index)
    img_index = 0;
end

output_type = sampling_args.output_type;
switch sampling_args.image_type
    case {'real'}
        if ~exist('image_lists','var')
            error('Must supply an image_list when using real images');
        end
        
        img_index = img_index + 1;
        
        % If we've gone past the end of the list, return empty values
        if (img_index > length(image_lists))
            img = []; img_label = [];
            fg_map = []; bg_map = [];
            parameters = [];
            return;
        end
        
        % Retrieve the next image from the list
        img = u_load(image_lists(img_index).image);
        
        % Get the ground truth labelling
        img_label = get_real_image_label(image_lists(img_index), output_type);
        
        % Get foreground/background pixel maps
        [fg_map, bg_map] = get_real_masks(image_lists(img_index), ... 
                                          sampling_args, output_type, img_label);
        
        % Empty parameter variable
        parameters = [];
        
    case {'line','grain'}
        % Synthesize a new image
        [img, line_map, centre_map, orientation_map, width_map, parameters] = ...
            create_synthetic_image(sampling_args);
        
        % Choose the appropriate label image
        img_label = get_synthetic_image_label(line_map, centre_map, ...
                                              orientation_map, width_map, ...
                                              output_type);
        
        % Get foreground/background pixel maps
        [fg_map, bg_map] = ...
            get_synthetic_masks(line_map, centre_map, ...
                                sampling_args, output_type);
        
    otherwise
        error(['image type ', sampling_args.image_type, ' not recognised']);
end


%% Test script
function test_script()
clc;

args = default_args();
args.output_type = 'detection';

args.image_type = 'real';
sampling_args = get_sampling_args_from(args);
sampling_args = get_output_args_from(args, sampling_args);
image_lists = create_image_lists(sampling_args);
image_lists(1)

get_next_image('reset');
[img, img_label, fg_map, bg_map, parameters] = ...
    get_next_image(sampling_args, image_lists);
while ~isempty(img)
    figure(1); clf; colormap(gray(256));
    subplot(2,2,1); imagesc(img); axis('image','ij');
    subplot(2,2,2); imagesc(abs(img_label)); axis('image','ij');
    subplot(2,2,3); image(fg_map*255); axis('image','ij');
    if ~isempty(bg_map)
        subplot(2,2,4); image(bg_map*255); axis('image','ij');
    end
    disp(parameters);
    pause(0.5);
    
    [img, img_label, fg_map, bg_map, parameters] = ...
        get_next_image(sampling_args, image_lists);
end

return

args.image_type = 'line';
args.bg_dir = [asymmetryroot,'\data\synthetic_backgrounds\smooth512\train/'];
args.bg_zeros = 5;
args.num_bgs = 10;

[img, img_label, fg_map, bg_map, parameters] = get_next_image(args);
figure(1); clf; colormap(gray(256));
subplot(2,2,1); imagesc(img); axis('image','ij');
subplot(2,2,2); imagesc(abs(img_label)); axis('image','ij');
subplot(2,2,3); image(fg_map*255); axis('image','ij');
if ~isempty(bg_map)
    subplot(2,2,4); image(bg_map*255); axis('image','ij');
end
disp(parameters);
