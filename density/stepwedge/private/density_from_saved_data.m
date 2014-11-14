function [gland_thickness breast_thickness taper_region_idx] =...
    density_from_saved_data(filename, data_path, image_path, resize_factor, debug_mode)
%DENSITY_FROM_SAVED_DATA *Insert a one line summary here*
%   [] = density_from_saved_data(filename,filepath)
%
% Inputs:
%      filename- The name of the mammogram to be processed
%
%      save_directory_name- The filepath to the directory containing the existing user
%       data for this mammogram
%
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 08-May-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%--------------------------------------------------------------------------
%Set all the constants as in the main stepwedge script
leftright = 0; 
%-------------------------------------------------------------------------

%Set default for debug_mode
if nargin < 5
    debug_mode = 0;
end

%Make sure directory paths are file separated
if ~strcmp(data_path(end), filesep);
    data_path = [data_path filesep];
end
if ~strcmp(data_path(end), filesep);
    data_path = [data_path filesep];
end

%Generate the data filename
mammo_filename = [data_path, filename(1:end-4), '_data'];

%load in the existing data
density_data = u_load(mammo_filename);

% define film size
% 18x24 and 24x30 are digitised at different orientations
% 18x24 films will need to be flipped
if ~isempty(strfind(filename, '1824'))
    filmsizes = 1;
elseif ~isempty(strfind(filename, '2430'))
    filmsizes = 2;
else
    %define error action
end

% read in the image
[IMAGE] = imread([image_path, filename],'tif');

if filmsizes == 1
    % this is only needed for the 1824 film sizes
    IMAGE = imrotate(IMAGE,90);
end

% reduce image sizes here (don't need such high resolution now magnification markers have been located)
IMAGE = imresize(IMAGE,resize_factor);
dimensions = size(IMAGE);

%Get thickness profile across the whole mammogram
x_b = thickness_from_polyfit(density_data.x_b_info, dimensions, resize_factor, debug_mode);

% Now reduce the image to 4096 grey levels and invert to make it compatible
% with sw_lookup_table.m
IMAGE = IMAGE./16;
IMAGE = 4095-IMAGE;

% create lookup table to convert pixel value to x_sw
sw_lookup = sw_lookup_table_12bit(density_data.wedgevals);

% get mask of breast area
[background_mask, edgex, edgey] = breast_edge(IMAGE,...
    density_data.coarse_edgex, density_data.coarse_edgey, leftright, debug_mode);
background_mask = double(background_mask > 0);

%now call edge_profile to correct the breast edge mask, applying a
%semicircular profile and returning a thickness map
[breast_thickness taper_region_idx] =...
            mb_edge_profile(x_b,background_mask,edgex,edgey,leftright,resize_factor, debug_mode);

% then produce x_gland image
gland_thickness = sw_vs_gland_thick(sw_lookup(IMAGE+1), breast_thickness);
figure; 
subplot(1,2,1); 
imagesc(-double(IMAGE)); axis image; colormap(gray(256)); 
title('Original mammogram');

subplot(1,2,2);        
imagesc(gland_thickness); axis image; colormap(gray(256)); 
colorbar('west', 'xcolor', 'r', 'ycolor', 'r');
title('Gland thickness map');

% gland_thickness array contains NaN - so convert NaN to 0 or 
% volume calculation will give NaN
gland_thickness(isnan(gland_thickness)) = 0;

