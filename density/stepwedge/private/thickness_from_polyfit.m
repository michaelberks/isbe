function [x_b_matrix] = thickness_from_polyfit(x_b_info,dimensions,resize_factor,debug_mode)
%THICKNESS_FROM_POLYFIT create breat thickness matrix from marker pairs
%assuming paddle can't flex and thicknesses should therefore be a straight
%line
%   [x_b_matrix] = thickness_from_polyfit(x_b_info,dimensions,resize_factor)
%
% Inputs:
%      x_b_info - 2xN array, where N is the number of marker pairs
%                    1st column = x-coordinates of the marker pairs
%                    2nd column = associated thickness
%
%      dimensions - The size of the thickness matrix to create
%
%      resize_factor - The scaling factor between x-coordinates in x_b_info
%                       and dimensions
%
%      debug_mode - Switch debugging images on/off (default off)
%
% Outputs:
%      x_b_matrix - Matrix of breast thickness at each pixel
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 16-Apr-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

if nargin < 4
    debug_mode = 0;
end

% Use polyfit to fit a straight line to the thickness data (scaling the
% x-coordinates in x_b_info)
tilt_line = polyfit(x_b_info(:,1)*resize_factor,x_b_info(:,2),1);

%Now compute the breast thickness for a single row
row_thickness = tilt_line(1) * (1:dimensions(2)) + tilt_line(2);

%Now copy the sibgle row to all rows in the output matrix
x_b_matrix = repmat(row_thickness, dimensions(1), 1);

if debug_mode
    
    %Display the gradient of the paddle tilt
    disp('tilt gradient   offset    width');
    format;
    disp([tilt_line(1) tilt_line(2) dimensions(2)]);
    format bank;
    
    %Plot the input thickness against the fitted line
    figure('Name','Thickness profile');
    hold on
    plot(x_b_info(:,1)*resize_factor,x_b_info(:,2), 'b');
    plot(x_b_info(:,1)*resize_factor,x_b_info(:,2),'or');
    plot(1:dimensions(2), x_b_matrix(1,:), 'g')
    hold off
    
    %Display the resulting thickness matrix
    figure('Name','x_b image');
    imagesc(x_b_matrix); axis image; colormap(gray(256));
    disp(['min of x_b returned from thickness = ', num2str(min(x_b_matrix(:)))]);
    disp(['max of x_b returned from thickness = ', num2str(max(x_b_matrix(:)))]);
    
end
