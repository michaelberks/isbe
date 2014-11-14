function [img, line_map, centre_map, orientation_map, width_map, params] = ...
    create_synthetic_image(sampling_args)

f_debug = (nargin == 0 && nargout == 0);
if f_debug, test_script(); return; end

[img, line_map, centre_map, orientation_map, width_map, params] = func(sampling_args);


%% The function
function [img, line_map, centre_map, orientation_map, width_map, params] = ...
    func(sampling_args)

switch sampling_args.image_type
    case 'line',
        % generate the image
        [img, line_map, centre_map, orientation_map, params] = ...
            create_line_image(sampling_args);

        % Width defined only at the centre of the line
        width_map = params.width * double(centre_map);
        
    case 'grain',
        % generate a 'grain' image
        [img, line_map, centre_map, orientation_map,  params] = ...
            create_grain_image(sampling_args);
        
        % Grain images have no defined width
        width_map = nan(size(line_map));

    otherwise,
        error(['Synthetic image type ', sampling_args.image_type, ' not recognized']);
end

if ~isempty(sampling_args.noise_type)
    img = add_image_noise(img, sampling_args.noise_type, sampling_args.noise_params);
end

% Ensure that orientation is given as a complex valued unit vector whose
% angle is twice the orientation
orientation_map = complex(cosd(2*orientation_map), sind(2*orientation_map));


%% Test script
function test_script()
clc;

rootpath = [asymmetryroot,'\data\synthetic_backgrounds\'];
sampling_args.bg_dir = [rootpath,'smooth512\train/'];
sampling_args.bg_fmt = 'mat';
sampling_args.bg_stem = 'bg';
sampling_args.bg_zeros = 5;
sampling_args.num_bgs = 10;
sampling_args.image_type = 'line';
sampling_args.orientation_range = [0,180];
sampling_args.width_range = [4 16]*4;
sampling_args.contrast_range = [4 8];
sampling_args.decay_rate = 4;
sampling_args.line_type = 'sin';

[img, line_map, centre_map, orientation_map, width_map, params] = ...
    create_synthetic_image(sampling_args);
figure(1); clf; colormap(gray(256));
subplot(2,2,1); imagesc(img); axis('image','ij');
subplot(2,2,2); image(line_map*255); axis('image','ij');
subplot(2,2,3); image(centre_map*255); axis('image','ij');
disp(params);

return

sampling_args.image_type = 'grain';
[img, line_map, centre_map, orientation_map, width_map, params] = ...
    create_synthetic_image(sampling_args);
figure(2); clf; hold on; colormap(gray(256));
imagesc(img); axis('image','ij');
disp(params);
