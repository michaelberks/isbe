function [ f_image, f_r_map,delete_range ] = remove_white_frame( image,r_map )
%REMOVE_WHITE_FRAME : Removes any white frame caused by the scanning in of mammograms reads in image and if one of the all edges has more than 
%					  "fraction" filled with pure white (or nearly - defined by GS_threshold), then remove edge from image. Remove points are 
% 					  also removed from region map (r_map). 
% 
%
% Inputs:
%			image				image that you want to find edges in [2D array]
%			r_map		 		region map [2D array]
%
% Outputs:
%			f_image				cut down image [2D array]
%			f_r_map 			cut down region map [2D array]
%			delete_range 		number of pixel columns or rows removed from the image edges during function
%
% Example:
%
% Notes:	I think it would be best to make the white frame around the image a mask rather than removing it. As the code stands at the moment
%			the marker co-ordinates are in the frame of the cut image and not the original. 
%			
% See also:
%
% Created: 18-06-2013
% Author: Euan Allen
% Email : euan.allen@gmail.com

GS_threshold = 7536; %If GS value is >58000
delete_numL = 0; %Components of "delete_range"
delete_numR = 0;
delete_numT = 0;
delete_numB = 0;
N = 99;%how many pixels into image to check
[Y X] = size(image);
fraction = 0.3; %Fraction of border that needs to be white.

%left side
for i = 1:N %search N edges in
    temp = image(:,i); %generate 1D temp array of strip (single row or column) to search in
    count = temp > 65536 - GS_threshold; %count the number of pixels in temp with a value > 65536 - GS_threshold
    if sum(count) > fraction*Y; delete_numL=delete_numL+1; end %if count is >fraction then delete it
end
clear temp count

%right side
for i = X-N-1:X %last N edges in
    temp = image(:,i);
    count = temp > 65536 - GS_threshold;
    if sum(count) > fraction*Y;  delete_numR=delete_numR+1; end
end
clear temp count

%top side
for i = 1:N %search N edges in
    temp = image(i,:);
    count = temp > 65536 - GS_threshold;
    if sum(count) > fraction*X; delete_numT=delete_numT+1; end
end
clear temp count

%bottom side
for i = Y-N-1:Y %last N edges in
    temp = image(i,:);
    count = temp > 65536 - GS_threshold;
    if sum(count) > fraction*X;  delete_numB=delete_numB+1; end
end

f_image = image(delete_numT+1:Y-delete_numB,delete_numL+1:X-delete_numR);
f_r_map = r_map(delete_numT+1:Y-delete_numB,delete_numL+1:X-delete_numR);
delete_range = [delete_numL,delete_numR,delete_numT,delete_numB];
end

