function [image_out, label, parameters] = create_bar_image(varargin)
%SAMPLE_IMAGE_TRAINING_DATA *Insert a one line summary here*
%   [] = sample_image_training_data(varargin)
%
% SAMPLE_IMAGE_TRAINING_DATA uses the U_PACKARGS interface function
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
			 '0', ... % non-strict mode, no mandatory arguments
			 'bg_dir', 'C:\isbe\dev\classification\data\normal_smooth128\',...
             'width_range', [4 16],...
             'orientation_range', [0 180],...
             'contrast_range', [8 16],...
             'plot', 0);

%get directory listing of backgrounds
bg_list = dir([args.bg_dir, '*.mat']);

%randomly select background - will load a matrix 'bg'
bg_idx = ceil(length(bg_list)*rand);
bg = u_load([args.bg_dir, bg_list(bg_idx).name]);

%get size of background
[row col] = size(bg);

parameters.bg_idx = bg_idx;

%randomly select orientation between 0 and 180 degrees
orientation = args.orientation_range(1) +...
    (args.orientation_range(2)-args.orientation_range(1))*rand;

%randomly select width between 4 and 32
    width = args.width_range(1) + (args.width_range(2)-args.width_range(1))*rand;
    
%randomly select contrast between 8 and 16
    contrast = args.contrast_range(1) + (args.contrast_range(2)-args.contrast_range(1))*rand;
    
%randomly select bar type
g_bar = rand > .5;

if g_bar    
    %make gaussian bar
    [bar_image, label] = create_gauss_bar(width/2, contrast, orientation, row, col);
    label = uint8(label);
else       
    %make rectangular bar
    [bar_image, label] = create_rect_bar(width/2, contrast, orientation, row, col);
    label = uint8(label); %label both bars the same
end

%Add bar to background
image_out = bar_image + bg;

if args.plot
    figure;
    subplot(1,3,1); imagesc(bg); axis image; colormap(gray(256));
    subplot(1,3,2); imagesc(bar_image); axis image; colormap(gray(256));
    subplot(1,3,3); imagesc(image_out); axis image; colormap(gray(256));
end

%Save output parameters
parameters.width = width;
parameters.contrast = contrast;
parameters.orientation = orientation;
parameters.g_bar = g_bar;

%--------------------------------------------------------------------------
return