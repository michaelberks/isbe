function [filtered_image] = median_filter_circle3(image_in, circle_radius)
%MEDIAN_FILTER_CIRCLE1 *Insert a one line summary here*
%   [filtered_image] = median_filter_circle1(image_in, circle_radius)
%
% Inputs:
%      image_in - image to be filtered
%
%      circle_radius - radius about each pixel to take the median value
%
%
% Outputs:
%      filtered_image - median filtered image
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 19-Mar-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

%Get image size
[rows cols] = size(image_in);

%Pre-allocate filtered image
filtered_image = zeros(rows, cols);

%work out size circle matrix needs to be
circle_radius_int = ceil(circle_radius);

%Make circle mask
x = -circle_radius_int:circle_radius_int;
xy = repmat(x, 2*circle_radius_int+1,1);
circle_bw = (xy.^2 + xy'.^2) < circle_radius^2;
        
%Pad the image first
image_in = padarray(double(image_in), [circle_radius_int circle_radius_int], nan);

for i_row = 1:rows
    for i_col = 1:cols
        
        %Get patch from image at this point
        image_patch = image_in(i_row+x+circle_radius_int, i_col+x+circle_radius_int);
        
        %Make sure we deal with the edge of the image
        circle_bwi = circle_bw & ~isnan(image_patch);
                
        %Get values using the mask
        med_vals = image_patch(circle_bwi);

        %Save the median value in our output variable
        filtered_image(i_row, i_col) = median(med_vals);
        
    end
end
        
        


