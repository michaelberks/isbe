function [pectoral_mask pectoral_x pectoral_y] = pectoral_user(IMAGE)

% manually mark the pectoral muscle

f1 = figure('Name', 'Select approximate breast boundary');
imagesc(padarray(IMAGE,[20 20],max(IMAGE(:)))); colormap(gray(256)); axis image;

title('Mark the pectoral muscle');
set(gca, 'xticklabel', [], 'yticklabel', []);
[mask, pectoral_x, pectoral_y] = roipoly; %#ok xi yi are saved
close(f1);
pectoral_x = pectoral_x(1:end-1) - 20;
pectoral_y = pectoral_y(1:end-1) - 20;

pectoral_mask = roipoly(IMAGE, pectoral_x, pectoral_y);