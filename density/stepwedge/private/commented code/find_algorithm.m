function [markers] = find_algorithm(image,x_s,y_s,zero_region,number_of_peaks,allow_map,mask,image_display)
%FIND_ALGORITHM : Given image strip this function finds the marker positions via normalised cross correlation. The code outputs the highest peaks
%				  in the image with the number specified by number_of_peaks.
% 
%
% Inputs:
%			image				image for search [2D array]
%			x_s 				gives the x-split of the edge strips taken for analysis within a single image (varies with image size) [1D array]
%			y_s 				gives the y-split of the edge strips taken for analysis within a single image (varies with image size) [1D array]
%			zero_region			radius of zero region generated around found peaks [#]
%			number_of_peaks		number of markers to output in the image [#]
%			allow_map			map of where markers may exist (regions) [2D array]
%			mask				mask of disallowed regions due to anon region
%			image_display		set to 1 if you want to display results as an image, 0 if not
%
% Outputs:
%			markers 			output of the markers found. Each marker has the data:
%								x-position/y-position/strip number/region ID number output[2D array] 
%
% Example:
%
% Notes:   This function also requires access to 'patch3.mat' (template).
%			
% See also:
%
% Created: 18-06-2013
% Author: Euan Allen
% Email : euan.allen@gmail.com

count1 = 1; count2 =1; %used for marker information storage

%Perform template matching on image strip (defined by x_s/y_s) to find top marker peaks
[markers1, mValues1] = template_matching(image(y_s(1):y_s(2), x_s(1):x_s(2)), ...
    'patch3.mat', zero_region, number_of_peaks, logical(mask(y_s(1):y_s(2),x_s(1):x_s(2))), image_display);

for i=1:length(markers1)%Go through markers 
    if  allow_map(markers1(i,2)+y_s(1)-1,markers1(i,1)+x_s(1)-1)>0; %Do they lie in the allowed regions
        g_markers1(count1,1) = markers1(i,1)+x_s(1); %input position (in original image frame, hence x_s/y_s) into output
        g_markers1(count1,2) = markers1(i,2)+y_s(1);
        g_markers1(count1,3) = 1; %Because its using the first set of coordinates in x_s/y_s (i.e. is strip 1 of the image edge containing markers)
        g_markers1(count1,4) = allow_map(markers1(i,2)+y_s(1)-1,markers1(i,1)+x_s(1)-1); %input specific region number marker lies in (see region_definitions)
        
		%Set allow_map to zero around as marker has now been found in that region (so no more markers can be found in that region)
        allow_map=zero_region_square(allow_map,markers1(i,1)+x_s(1)-1,markers1(i,2)+y_s(1)-1,500);
        count1 = count1+1;
    end
end

[markers2, mValues2] = template_matching(image(y_s(3):y_s(4),x_s(3):x_s(4)),'patch3.mat',zero_region,number_of_peaks,logical(mask(y_s(3):y_s(4),x_s(3):x_s(4))),image_display);
for i=1:length(markers1) %Are the markers in a region
    if  allow_map(markers2(i,2)+y_s(3)-1,markers2(i,1)+x_s(3)-1)>0; %if marker is in allowed region
        g_markers2(count2,1) = markers2(i,1)+x_s(3); %moves into whole image frame
        g_markers2(count2,2) = markers2(i,2)+y_s(3);
        g_markers2(count2,3) = 2; %Because its using the second set of coordinates in x_s/y_s
        g_markers2(count2,4) = allow_map(markers2(i,2)+y_s(3)-1,markers2(i,1)+x_s(3)-1);
        %Set allow_map to zero around as marker has now been found in that region.
        allow_map=zero_region_square(allow_map,markers2(i,1)+x_s(3)-1,markers2(i,2)+y_s(3)-1,500);
        count2 = count2+1;
    end
end

%if image_display; figure; imshow(allow_map);figure; imshow(mask);end

if exist('g_markers1','var')&&exist('g_markers2','var') %if markers have been found in both strips
    markers = [g_markers1;g_markers2]; %output them both
elseif exist('g_markers1','var') 
    markers = g_markers1;
elseif exist('g_markers2','var')
    markers = g_markers2;
else %if no markers found at all
    markers = [0,0,0,0]; %output zero result
end

end
