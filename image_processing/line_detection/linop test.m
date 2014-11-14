%Linop_test
% Script to test 2D line operator.
%Initially set up to threshold the chromosome image, chromos1.tif.
%Make sure it is in the path.
% Edit the source to change the threshold or use other images 

close all; % get rid of figures from previous runs
iptsetpref('ImshowBorder','tight'); % display preferences
src = imread('102_2.tif'); % load the image; 
figure, imshow(src); % display the image
num_angles = 24; % number of orientations for detecting lines
min_scale = 5; % smallest scale for detecting linear segments;
max_scale = 15; % largest scale for detecting linear segments;
[response, orientations] = line_operator_conv(src, num_angles,min_scale, max_scale, 'degrees');
% response is an image containing the line operator response, maximised across scales. 
% orientations is an image with the orientations of the maximum responses.
%           bin = (src > 50);% binary image; change the threshold here
figure, imshow(response);% display the response image
line_im = non_maximal_supp(response, orientations)

figure, imshow(orientations);
figure, imshow(line_im);