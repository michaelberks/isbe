function [] = generate_bar_images(varargin)
%GENERATE_BAR_IMAGES *Insert a one line summary here*
%   [] = generate_bar_images(varargin)
%
% SAMPLE_BAR_TRAINING_DATA uses the U_PACKARGS interface function
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
% Created: 28-Jan-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
			 '0', ... % non-strict mode
			 {'num_images',...
             'image_dir'},... % the mandatory arguments
             'bg_dir', 'C:\isbe\dev\classification\data\normal_smooth128\',...
             'spicule_ratio', 0.5,...
             'width_range', [2 16],...
             'orientation_range', [0 180],...
             'contrast_range', [8 16],...
             'compute_dt', 0,...
             'plot', 0);

if ~isempty(args.image_dir) && ~strcmp(args.image_dir(end), filesep)
    args.image_dir = [args.image_dir filesep];
end
mkdir(args.image_dir);

%get directory listing of backgrounds
bg_list = dir([args.bg_dir, '*.mat']);


%Create images 
for ii = 1:args.num_images

    %randomly select background - will load a matrix 'bg'
    bg_idx = ceil(length(bg_list)*rand);
    bg = u_load([args.bg_dir, bg_list(bg_idx).name]);

    %get size of background
    [row col] = size(bg);

    %randomly select orientation between 0 and 180 degrees
    orientation = args.orientation_range(1) +...
        (args.orientation_range(2)-args.orientation_range(1))*rand;

    %randomly select width between 4 and 32
    width = args.width_range(1) + (args.width_range(2)-args.width_range(1))*rand;

    %randomly select contrast between 8 and 16
    contrast = args.contrast_range(1) + (args.contrast_range(2)-args.contrast_range(1))*rand;

    %randomly select bar type
    spicule = rand < args.spicule_ratio;

    if spicule    
        %make gaussian bar
        [bar_image, label, label_centre] = create_ellipse_spicule_mb(width/2, contrast, orientation, row, col);
    else       
        %make rectangular bar
        [bar_image, label, label_centre] = create_ellipse_bar_mb(width/2, contrast, orientation, row, col);
    end
    clear dummy;
    
    %Store bar parameters
    parameters.bg_idx = bg_idx; %#ok
    parameters.width = width; %#ok
    parameters.contrast = contrast; %#ok
    parameters.orientation = orientation; %#ok
    parameters.spicule = spicule; %#ok
    
    %Add bar to background
    bar_image = bar_image + bg;
    
    %if we've been asked to, plot image and labels
    if args.plot
        figure;
        subplot(1,3,1); image(bar_image); axis image; colormap(gray(256));
        subplot(1,3,2); imagesc(label); axis image;
        subplot(1,3,3); imagesc(label_centre); axis image;
    end
    
    %if we've been asked to, compute DT-CWT representation
    max_dt = []; %#ok
    if args.compute_dt
        % Create DT-CWT of image
        dt = dtwavexfm2(bar_image, args.num_levels, 'near_sym_b', 'qshift_b');

        %interpolation DT-CWT to full pixel grid
        dt_inter = dt_to_full_image(dt); clear dt;

        %get the maximum response across orientations
        max_dt = squeeze(max(dt_inter, [], 3)); clear dt; %#ok
    end  
    
    %save the image data
    save([args.image_dir 'bar', zerostr(ii,3)], 'bar_image', 'label', 'label_centre', 'parameters', 'max_dt');
    display(['Saved image ', num2str(ii)]);
end
        