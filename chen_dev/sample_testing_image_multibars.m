function [] = sample_testing_image_multibars(varargin)
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
    {'num_images'},... % the mandatory arguments
    'bg_dir', 'E:\DTCWTmini\data\smooth512x512_patches\test\',...
    'width_range', [4 16],...
    'orientation_range', [0 180],...
    'contrast_range', [1 8],...
    'max_numlines', 20,...
    'bar_type', 'ellipse',...
    'save_path', 'E:\DTCWTmini\data\testimage_contrast1to8_exprnd_sin\',...
    'save_para', 'parameters',...
    'plot', 0);

num_images = args.num_images; %total number of feature vectors to sample
save_path = args.save_path;
save_para = args.save_para;
max_numlines = args.max_numlines;

curr_image = 1;

%Create testing images
while curr_image < num_images+1
    disp(curr_image);
    %randomly select background - will load a matrix 'bg'
    bg_idx = ceil(4903*rand); %% 6130 bg pathces
    bg = u_load([args.bg_dir, 'bg', zerostr(bg_idx,5), '.mat']);
    
    %get size of background
    [row col] = size(bg);
    multibar_image = zeros(row, col);
    multibar_label = zeros(row,col);
    multibar_labelcentre = zeros(row,col);
    num_lines = ceil(max_numlines*rand);
    curr_line = 1;
    while curr_line < num_lines + 1
        
        %randomly select orientation between 0 and 180 degrees
        orientation = args.orientation_range(1) +...
            (args.orientation_range(2)-args.orientation_range(1))*rand;
        
        %randomly select width between 4 and 32
        width = args.width_range(1) + (args.width_range(2)-args.width_range(1))*rand;
        
        %%randomly select contrast between 8 and 16
        %contrast = args.contrast_range(1) + (args.contrast_range(2)-args.contrast_range(1))*rand;
        
        %randomly select contrast on exponential distribution
        decay_rate = 4;
        mu = (args.contrast_range(2) - args.contrast_range(1)) / (2*log(decay_rate));
        contrast = args.contrast_range(1) + exp_rand_sample(mu);
        
        %randomly select squashiness of the line
        squash = rand;
        
        %randomly select contrast between 8 and 16
        centre_x = ceil(col * rand);
        centre_y = ceil(row * rand);
        
        %make  bar
        switch args.bar_type
            case 'ellipse'
                [bar_image, label, label_centre] = create_ellipse_bar(width/2, contrast, orientation, row, col, centre_x, centre_y);
                
            case 'sin'
                [bar_image, label, label_centre] = create_sin_bar(width/2, contrast, orientation, row, col, squash, centre_x, centre_y);
                
            otherwise
                warning(['Bar type: ', args.bar_type, ' not recognised']); %#ok
                [bar_image, label, label_centre] = create_ellipse_bar(width/2, contrast, orientation, row, col, centre_x, centre_y);
        end
        
        %Add bar to simulate superimposition of structures
        multibar_image = multibar_image + bar_image;
        
        %Add bar labels
        multibar_labelcentre = multibar_labelcentre | label_centre;
        multibar_label = multibar_label | label;
        
        
        %multibar_labelcentre = multibar_labelcentre + label_centre;
        %multibar_label = multibar_label + label;
        %cross_indxcen = (multibar_labelcentre == 2);
        %cross_indx = (multibar_label == 2);
        %multibar_image(cross_indx) = max(multibar_image_old(cross_indx), bar_image(cross_indx));
        %multibar_label(cross_indx) = 1;
        %multibar_labelcentre(cross_indxcen) = 1;
        
        curr_para(curr_line). halfwidth = width/2; %#ok
        curr_para(curr_line). contrast = contrast;  %#ok
        curr_para(curr_line). orientation = orientation;  %#ok
        curr_para(curr_line). row = row;  %#ok
        curr_para(curr_line). col = col;  %#ok
        curr_para(curr_line). centre_x = centre_x;  %#ok
        curr_para(curr_line). centre_y = centre_y;  %#ok
        curr_para(curr_line). squash = squash;
        curr_line = curr_line + 1;
    end
    
    %Add bar to background
    image_in = multibar_image + bg;
    
    %Store bar parameters
    parameters(curr_image).bg_idx = bg_idx; %#ok
    parameters(curr_image).curr_para = curr_para; %#ok
    
    
    if args.plot
        figure;
        subplot(1,2,1); imagesc(image_in); axis image; colormap(gray(256)); hold on;
        subplot(1,2,2); imagesc(multibar_labelcentre); axis image; colormap(gray(256)); hold on;
    end
    
    %     save([args.save_path, sprintf('test_image%d.mat', curr_image)], 'test_image');
    
    %increment image count
    if ~isempty(save_path)
        if ~exist(save_path, 'dir')
            mkdir(save_path);
        end
        save([save_path, ['image', int2str(curr_image), '.mat']], 'image_in');
        save([save_path, ['label', int2str(curr_image), '.mat']], 'multibar_label', 'multibar_labelcentre');
    end
    
    curr_image = curr_image + 1;
end

save([args.save_path, save_para, '.mat'], 'parameters');

disp('finished!!');

return
