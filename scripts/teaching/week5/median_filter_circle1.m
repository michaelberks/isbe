function [filtered_image] = median_filter_circle1(image_in, circle_radius)
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

for i_row = 1:rows
    for i_col = 1:cols
        
        %work out size circle matrix needs to be
        circle_radius_int = ceil(circle_radius);
        med_vals = [];
        for i_x = -circle_radius_int:circle_radius_int
            for i_y = -circle_radius_int:circle_radius_int
                
                %First check we're within the radius
                if (i_x^2 + i_y^2) < circle_radius^2
                    
                    %Get the row column subscripts
                    i_xd = i_col - i_x;
                    i_yd = i_row - i_y;
                
                    %Check we're still in the image
                    if (i_xd >= 1) && (i_xd <= cols) && (i_yd >= 1) && (i_yd <= rows)
                        %Add the image value to the med_vals
                        med_vals(end+1) = image_in(i_yd, i_xd);
                    end
                end
            end
        end
        %Save the median value in our output variable
        filtered_image(i_row, i_col) = median(med_vals);
        
    end
end
        
        


