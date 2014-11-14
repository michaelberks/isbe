function [] = generate_synthetic_curve_images(varargin)
%SAMPLE_BAR_TRAINING_DATA *Insert a one line summary here*
%   [] = sample_bar_training_data(varargin)
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
    'save_dir'},... % the mandatory arguments
    'bg_dir', [asymmetryroot 'data\synthetic_backgrounds\real512\'],...
    'num_bgs', [],...
    'bg_stem', [],...
    'bg_zeros', [],...
    'width_range', [4 16],...
    'orientation_range', [0 360],...
    'contrast_range', [2 8],...
    'decay_rate', 4,...
    'max_numlines', 20,...
    'line_type', 'curve',...
    'save_path', [asymmetryroot '\data\synthetic_lines\'],...
    'save_para', 'parameters',...
    'plot', 0);

num_images = args.num_images; %total number of feature vectors to sample
save_path = [args.save_path args.save_dir];
save_para = args.save_para;
max_numlines = args.max_numlines;

%Workout whether we need to search for the listing of background patches
if isempty(args.bg_stem) || isempty(args.num_bgs)
    %get directory listing of backgrounds
    bg_list = dir([args.bg_dir, '*.mat']);
    args.num_bgs = length(bg_list);
elseif isempty(args.bg_zeros) 
    args.bg_zeros = floor(log10(args.num_bgs)) + 1;
end

%Create testing images
for curr_image = 1:num_images
    disp(curr_image);
    
    %randomly select background - will load a matrix 'bg'
    bg_idx = ceil(args.num_bgs*rand);
    if exist('bg_list', 'var')
        bg = u_load([args.bg_dir bg_list(bg_idx).name]);
    else   
        bg = u_load([args.bg_dir args.bg_stem zerostr(bg_idx, args.bg_zeros), '.mat']);
    end
    
    %get size of background
    [row col] = size(bg);
    test_image = zeros(row, col);
    label = zeros(row,col);
    label_centre = zeros(row,col);
    label_orientation = nan(row,col);
    num_lines = ceil(max_numlines*rand);
    curr_line = 1;
    while curr_line < num_lines + 1
        
        %randomly select orientation between 0 and 180 degrees
        orientation = args.orientation_range(1) +...
            (args.orientation_range(2)-args.orientation_range(1))*rand;
        
        %randomly select width
        %width = args.width_range(1) + (args.width_range(2)-args.width_range(1))*rand;
        mu = (args.width_range(2) - args.width_range(1)) / (2*log(args.decay_rate));
        width = args.width_range(1) + exp_rand_sample(mu);
    
        %%randomly select contrast between 8 and 16
        %contrast = args.contrast_range(1) + (args.contrast_range(2)-args.contrast_range(1))*rand;
        
        %randomly select contrast on exponential distribution
        mu = (args.contrast_range(2) - args.contrast_range(1)) / (2*log(args.decay_rate));
        contrast = args.contrast_range(1) + exp_rand_sample(mu);
        
        %randomly select squashiness of the line
        squash = 0.6*rand;
        
        %randomly select the radial curvature of the line
        radius = 2^(8*rand + 8);
        
        %randomly select contrast between 8 and 16
        centre_x = ceil(col * rand);
        centre_y = ceil(row * rand);
        
        %make  bar
        switch args.line_type
            case 'ellipse'
                [bar_image, bar_label, bar_centre, bar_orientation] = create_ellipse_bar(width/2, contrast, orientation, row, col, centre_x, centre_y);
                
            case 'sin'
                [bar_image, bar_label, bar_centre, bar_orientation] = create_sin_bar(width/2, contrast, orientation, row, col, squash, centre_x, centre_y);
                
            case 'curve'
                [bar_image, bar_label, bar_centre, bar_orientation] = create_sin_curve(width/2, contrast, radius, orientation, squash, row, col, centre_x, centre_y);
                
            otherwise
                warning(['Bar type: ', args.line_type, ' not recognised']); %#ok
                [bar_image, bar_label, bar_centre, bar_orientation] = create_ellipse_bar(width/2, contrast, orientation, row, col, centre_x, centre_y);
        end
        
        %Add bar to simulate superimposition of structures
        test_image = test_image + bar_image;
        
        %Add bar labels
        label_centre = label_centre | bar_centre;
        junctions = label & bar_label;
        label = label + bar_label;
        label_orientation(bar_label) = bar_orientation(bar_label);
        label_orientation(junctions) = nan;
        
        %label_centre = label_centre + bar_centre;
        %label = label + bar_label;
        %cross_indxcen = (label_centre == 2);
        %cross_indx = (label == 2);
        %test_image(cross_indx) = max(multibar_image_old(cross_indx), bar_image(cross_indx));
        %label(cross_indx) = 1;
        %label_centre(cross_indxcen) = 1;
        
        curr_para(curr_line).halfwidth = width/2; %#ok
        curr_para(curr_line).contrast = contrast;  %#ok
        curr_para(curr_line).orientation = orientation;  %#ok
        curr_para(curr_line).row = row;  %#ok
        curr_para(curr_line).col = col;  %#ok
        curr_para(curr_line).centre_x = centre_x;  %#ok
        curr_para(curr_line).centre_y = centre_y;  %#ok
        curr_para(curr_line).squash = squash;   %#ok
        curr_para(curr_line).radius = radius;   %#ok
        curr_line = curr_line + 1;
    end
    
    %Add bar to background
    test_image = test_image + bg;
    
    %Store bar parameters
    parameters(curr_image).bg_idx = bg_idx; %#ok
    parameters(curr_image).curr_para = curr_para; %#ok
    
    
    if args.plot
        figure;
        subplot(2,2,1); imagesc(test_image); axis image; colormap(gray(256)); hold on;
        subplot(2,2,2); imagesc(label); axis image; colormap(gray(256)); hold on;
        subplot(2,2,3); imagesc(label_centre); axis image; colormap(gray(256)); hold on;
        subplot(2,2,4); imagesc(label_orientation); axis image; colormap(gray(256)); hold on;
    end
    
    %     save([args.save_path, sprintf('test_image%d.mat', curr_image)], 'test_image');
    
    %increment image count
    if ~isempty(save_path)
        if ~exist(save_path, 'dir')
            mkdir(save_path);
            mkdir([save_path 'labels']);
        end
        save([save_path, 'image', zerostr(curr_image, floor(log10(num_images))+1), '.mat'], 'test_image');
        save([save_path, 'labels/label', zerostr(curr_image, floor(log10(num_images))+1), '.mat'], 'label', 'label_centre', 'label_orientation');
    end
end

save([save_path, save_para, '.mat'], 'parameters');

disp('finished!!');

return
