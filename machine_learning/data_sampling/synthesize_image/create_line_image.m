function [image_out, label, label_centre, label_orientation, parameters] = ...
	create_line_image(args)
%GENERATE_LINE_IMAGE Generate an image with a single line on a false
%background
%   [image_out, label, label_centre, label_orientation, parameters] = generate_line_image(bg, args)
%
% Inputs:
%      bg - Background on which to superimpose line
%
%      args - Parameters of distributions (e.g. width, contrast, etc)
%
%
% Outputs:
%      image_out - The resulting image
%
%      label - Binary map that indicates which pixels are part of a line
%
%      label_centre - Binary map indicating the centre of the line
%
%      label_orientation - The line orientation at each pixel
%
%      parameters - The specific parameters of this image
%
%
% Example:
%
% Notes:
%
% See also:
%		GENERATE_GRAIN_IMAGE
%
% Created: 02-Feb-2011
% Author: Phil Tresadern 
% Email : philip.tresadern@manchester.ac.uk 
% Phone : +44 (0)161 275 5114 
% Copyright: (C) University of Manchester 

%Make background
switch args.bg_type
    
    case 'load'
        % Load in a preprocessed mammographic background
        [bg, bg_index] = get_random_background(args);
        parameters.bg_index = bg_index;
        
        %get size of background
        [row col] = size(bg);
        
    case 'edge'
        row = args.bg_size(1);
        col = args.bg_size(2);
        %We'll actually make the background later
        
    case 'flat'
        row = args.bg_size(1);
        col = args.bg_size(2);
        bg = ones(row, col);
        
    otherwise
        % ?
end
        
    

%randomly select orientation
orientation = sample_uniform(args.orientation_range);

if args.decay_rate
    %randomly select width
    width = sample_rand_exp(args.width_range, args.decay_rate);

    %randomly select contrast on exponential distribution
    contrast = sample_rand_exp(args.contrast_range, args.decay_rate);

    %randomly select squashiness of the line
    if isfield(args,'squash_range')
        squash = sample_rand_exp(args.squash_range, args.decay_rate);
    else
        squash = rand;
    end
else
    %randomly select width
    width = sample_uniform(args.width_range);

    %randomly select contrast
    contrast = sample_uniform(args.contrast_range);

    %randomly select squashiness of the line
    if isfield(args,'squash_range')
        squash = sample_uniform(args.squash_range);
    else
        squash = rand;
    end
end

%Store bar parameters
parameters.width = width;
parameters.contrast = contrast;
parameters.orientation = orientation;
parameters.squash = squash;    

%make  bar
switch args.line_type
    case 'ellipse'
        [bar_image, label, label_centre, label_orientation] = ...
            create_ellipse_bar(width/2, contrast, orientation, ...
                               row, col, col/2, row/2);

    case 'sin'
        [bar_image, label, label_centre, label_orientation] = ...
            create_sin_bar(width/2, contrast, orientation, ...
                           row, col, squash, col/2, row/2);

    case 'curve'
        %randomly select the radial curvature of the line - this really
        %should be parameterised as an args option
        radius = 2^(8*rand + 8);
        parameters.radius = radius;
        [bar_image, label, label_centre, label_orientation] = ...
            create_sin_curve(width/2, contrast, radius, orientation, squash, ...
                             row, col, col/2, row/2);
    otherwise
        error(['Bar type: ', args.line_type, ' not recognised']);
end

%For edge backgrounds we make it now
if strcmp(args.bg_type, 'edge')
    %randomly select orientation
    orientation = sample_from_normal(orientation + 90, 30^2, 1);
    
    if args.decay_rate
        %randomly select width
        width = sample_rand_exp(args.width_range, args.decay_rate);

        %randomly select contrast on exponential distribution
        contrast = sample_rand_exp(args.contrast_range, args.decay_rate);      
    else
        %randomly select width
        width = sample_uniform(args.width_range);

        %randomly select contrast
        contrast = sample_uniform(args.contrast_range);
    end
    
    %Create edge
    bg = 1 + create_sin_step(width, contrast, orientation, row, col, col/2, row/2);
end
    


%Add bar to background
if isfield(args,'polarity') && strcmp(args.polarity,'neg')
	image_out = bar_image - bg;
else
	image_out = bar_image + bg;
end	
