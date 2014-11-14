function [image_out, skin_air] = segment_ferrari(image_in)
% segment_breast
% created: 07/03/2007 16:34
% author: Michael Berks
% function: segment the breast border using Ferrari et al (2004)
%           method

% 0) Check image_in is uint8 and resize
%%%%

[init_row init_col] = size(image_in);

if init_row == 1024
    row = 1024;
    col = init_col*row / init_row;
    image_in = imresize(image_in, [row col]);
else
    [row col] = size(image_in);
end

if ~isa(image_in, 'uint8')
    image_in = uint8(image_in);
end
    
% 1) logarithmic contrast enhancement
%%%%
 image_in = double(image_in);
 
 image_in = log(1 + image_in);
 image_in = 255 * ...
     (image_in - min(image_in(:)))/(max(image_in(:)) - min(image_in(:)));

 f1 = figure; image(image_in), colormap(gray(256)); axis image;
 
% 2) make mask using Lloyd-Max quantiser
%%%%
image_mask = lloyd_max(image_in);

% 3) Opening then select largest continous area
%%%%
image_mask = imopen(image_mask, strel('disk', 7, 0));

[ll, no_objects] = bwlabel(image_mask, 4);
[nn xx] = hist(ll(:), no_objects + 1); clear xx no_objects;
[mm ind] = max(nn(2:end)); clear max;
image_mask = not(ll - ind); clear ind ll;

f2 = figure; imagesc(image_mask), axis image;

% 4) Extract breast boundary
%%%%
rp = regionprops(bwlabel(image_mask, 4), 'Centroid');
[start_r start_c] = find(image_mask, 1);
breast_border = bwtraceboundary(image_mask, [start_r, start_c], 'E');
breast_border = fliplr(breast_border);

N1_x = min(breast_border(1)); N2_x = N1_x;
yy = breast_border(:,2);
N1_y = min(yy(breast_border(:,1) == N1_x));
N2_y = max(yy(breast_border(:,1) == N2_x));
idx2 = find((breast_border(:,1) == N2_x) & (breast_border(:,2) == N2_y)); 
d = (breast_border(:,1) - N2_x).^2 + (breast_border(:,2) - N2_y).^2;
idx3 = find(d == max(d));
N3_x = breast_border(idx3,1); N3_y = breast_border(idx3,2);

d = (breast_border(:,1) - N1_x).^2 + (breast_border(:,2) - N1_y).^2;
idx4 = find(d == max(d));
N4_x = breast_border(idx4,1); N4_y = breast_border(idx4,2);

skin_air = breast_border(idx3:idx2, :);

figure(f1); hold on;
plot(breast_border(:,1), breast_border(:,2), 'y', 'LineWidth', 1);
plot(skin_air(:,1), skin_air(:,2), 'r', 'LineWidth', 2);
plot(rp.Centroid(1), rp.Centroid(2), 'yx');
plot(N1_x, N1_y, 'mx', N2_x, N2_y, 'gx')
plot(N3_x, N3_y, 'mo', N4_x, N4_y, 'go')

% 5) Adaptive contour model
%%%%

skin_air2 = active_contour(skin_air(1:10:end,:), image_in);
plot(skin_air2(:,1), skin_air2(:,2), 'b', 'LineWidth', 1.5);
image_out = image_in;
