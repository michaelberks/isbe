centre_im = ones(128);
theta = 45;
[row col] = size(centre_im);
pad_size = ceil(0.5*(sqrt(2)-1)*max([row col]));

pad_im = padarray(centre_im, [pad_size pad_size], 0);

rot_im = imrotate(pad_im, theta, 'bilinear', 'crop');

sum_im = repmat(sum(rot_im), size(rot_im,1), 1);
sum_im2 = conv2(rot_im, ones(4*pad_size+2*row,1), 'same');

inv_im = imrotate(sum_im, -theta, 'bilinear', 'crop');
final_im = inv_im(pad_size+1:pad_size+row, pad_size+1:pad_size+col);

inv_im2 = imrotate(sum_im2, -theta, 'bilinear', 'crop');
final_im2 = inv_im2(pad_size+1:pad_size+row, pad_size+1:pad_size+col);

figure; imagesc(centre_im); axis image;
figure; imagesc(pad_im); axis image;
figure; imagesc(rot_im); axis image;
figure; subplot(1,2,1); imagesc(sum_im); axis image; subplot(1,2,2); imagesc(sum_im2); axis image;
figure; subplot(1,2,1); imagesc(final_im); axis image; subplot(1,2,2); imagesc(final_im2); axis image;
%%
mam_list = dir('C:\isbe\mammograms\new_CAD\BMP_2004_half\*.mat');
for ii = 1:540
    %load mammo 001LCC is good...
    load(['C:\isbe\mammograms\new_CAD\BMP_2004_half\o04_' mam_list(ii).name(5:10) '.mat']);

    %load segmentation
    load(['C:\isbe\dev\segmentation\breast_borders\o04_' mam_list(ii).name(5:10) '_segmentation.mat']);

    %resize segmentation and create mask
    [breast_border] = segment_breast_resize(size(mammogram), segmentation);
    mask = roipoly(mammogram, breast_border(:,1), breast_border(:,2));

    %Compute line strength
    class_forest = u_load('M:\chen\data\DTCWT_rf_fulltrees_W3L5\tree_combine\random_forest.mat');
    line_prob = classify_image(...
        'image_in', mammogram,...
        'forest', class_forest,...
        'forest_type', 'isbe',...
        'num_levels', 5,...
        'mask', mask);

    %Compute orientation
    reg_forest = u_load(['C:\isbe\dev\classification\rf\rf_reg_ori_02.mat']);
    line_ori = classify_image(...
        'image_in', mammogram,...
        'forest', reg_forest,...
        'forest_type', 'regression',...
        'num_levels', 5,...
        'mask', mask);

    save(['C:\isbe\dev\ad\' mam_list(ii).name(5:10) '_data.mat'], 'line*');

    line_combined = line_prob .* exp(i*pi*line_ori/180);
    figure; imagesc(line_prob); axis image; colormap(gray(256));
    figure; image(complex2rgb(line_combined(1:4:end, 1:4:end).^2)); axis image;
    clear; pack;
    load(['C:\isbe\dev\ad\' mam_list(ii).name(5:10) '_data.mat']);
    [angle_bands] = radial_line_projection(line_prob, line_ori, [36 1]);
    save(['C:\isbe\dev\ad\' mam_list(ii).name(5:10) '_data.mat'], 'line*', 'angle_bands');
end
%%
load C:\isbe\dev\ad\001LCC_data.mat
line_combined = line_prob .* exp(i*pi*line_ori/180);
imwrite(complex2rgb(line_combined(1:4:end, 1:4:end).^2), 'C:\isbe\dev\ad\figures\001LCC_orientation_map.bmp');
write_im_from_colormap(angle_bands, 'C:\isbe\dev\ad\figures\001LCC_radial_map_1.bmp', jet(256), [0 1100]);
write_im_from_colormap(imfilter(angle_bands, fspecial('disk', 8)), 'C:\isbe\dev\ad\figures\001LCC_radial_map_8.bmp', jet(256), [0 1100]);
write_im_from_colormap(imfilter(angle_bands, fspecial('disk', 16)), 'C:\isbe\dev\ad\figures\001LCC_radial_map_16.bmp', jet(256), [0 1100]);
write_im_from_colormap(imfilter(angle_bands, fspecial('disk', 32)), 'C:\isbe\dev\ad\figures\001LCC_radial_map_32.bmp', jet(256), [0 1100]);
write_im_from_colormap(imfilter(angle_bands, fspecial('disk', 64)), 'C:\isbe\dev\ad\figures\001LCC_radial_map_64.bmp', jet(256), [0 1100]);