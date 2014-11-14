function [coarse_edgex, coarse_edgey] = breast_edge_user(IMAGE)

% this function is to try and identify breast periphery areas

% reduce image a bit further anyway
further_resize = 0.5;
IMAGE_smaller = imresize(IMAGE,further_resize);

% select a region containing only the breast - avoid markers, wedge, etc...

f1 = figure('Name', 'Select approximate breast boundary');
imagesc(IMAGE_smaller); colormap(gray(256)); axis image;
title('Select ROI containing breast region, but no confusing markers, etc..');
ylabel('On nipple side define region loosely - program will autofind here');
xlabel('On chest-wall side define breast region tightly');
set(gca, 'xticklabel', [], 'yticklabel', []);
[mask, coarse_edgex, coarse_edgey] = roipoly; %#ok xi yi are saved
close(f1);