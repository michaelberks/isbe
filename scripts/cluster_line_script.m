%%
%make line masks

mask{1,1} = zeros(64);

mask{1,1}(9:56, 31:34) = 1;
mask{1,1} = logical(mask{1,1});

mask{2,1} = imrotate(mask{1,1}, 45, 'nearest', 'crop');
mask{3,1} = imrotate(mask{1,1}, 90, 'nearest', 'crop');
mask{4,1} = imrotate(mask{1,1}, 135, 'nearest', 'crop');

for ii = 1:4
    figure; imagesc(mask{ii}); colormap(gray(256)); axis image;
end

%% 
% make and save 10 images
mkdir('C:\isbe\dev\background\images\lines');
for ii = 1:10
    bg_int = uint8(100+round(30*rand(1)));
    
    line_ints = bg_int + 20 + uint8(10*randn(4,1));
    
    A = uint8(zeros(64));
    A(mask{1}) = line_ints(1);
    A(~mask{1}) = bg_int;
    
    B = uint8(zeros(64));
    B(mask{2}) = line_ints(2);
    B(~mask{2}) = bg_int;
    
    C = uint8(zeros(64));
    C(mask{3}) = line_ints(3);
    C(~mask{3}) = bg_int;
    
    D = uint8(zeros(64));
    D(mask{4}) = line_ints(4);
    D(~mask{4}) = bg_int;
    
    im = [A B; C D];
    im = padarray(im, [64 64], bg_int);
    im = im + uint8(5*randn(256));
    figure; image(im); colormap(gray(256)); axis image;
    imwrite(im, ['C:\isbe\dev\background\images\lines\line', zerostr(ii,3), '.bmp']);
    clear im;
end

%%
%Construct pyramids from images
build_args.ImageDir = 'C:\isbe\dev\background\images\lines\';
build_args.OutputDir = 'C:\isbe\dev\background\pyramid\lines\';
build_args.NumLevels = 5;
build_args.NumOrientations = 2;

mb_build_pyramids(build_args);

%%
%
% Clustering performed on hydra
%
%%

% Assign clusters from models to pyramid bands
input_list = dir('C:\isbe\dev\background\pyramid\lines\*pyramid*');

for ii = 1:length(input_list)
    pyramid = u_load(['C:\isbe\dev\background\pyramid\lines\', input_list(ii).name]);
    cluster_image = cell(7,2);
    
    for level = 2:5
        for ori = 1:2
            pyr_model = u_load(['C:\isbe\dev\background\results\lines\lines_model_',...
                num2str(level), '_', num2str(ori)]);

            [cluster_image{level, ori}] = assign_cluster_to_image(pyr_model, pyramid{level, ori});
            clear pyr_model;
        end   
    end
    c_name = ['C:\isbe\dev\background\misc\lines\', input_list(ii).name, 'cluster_image'];
    save(c_name, 'cluster_image');
    
    clear cluster_image pyramid
        
end
%%

i1 = imread('C:\isbe\dev\background\images\lines\line001.bmp');
i1 = double(i1);

i2 = i1;
i2(65:128, 65:128) = i1(129:192, 65:128);
i2(65:128, 129:192) = i1(65:128, 65:128);
i2(129:192, 129:192) = i1(65:128, 129:192);
i2(129:192, 65:128) = i1(129:192, 129:192);

i3 = i2;
i3(65:128, 65:128) = i2(129:192, 65:128);
i3(65:128, 129:192) = i2(65:128, 65:128);
i3(129:192, 129:192) = i2(65:128, 129:192);
i3(129:192, 65:128) = i2(129:192, 129:192);

i4 = i3;
i4(65:128, 65:128) = i3(129:192, 65:128);
i4(65:128, 129:192) = i3(65:128, 65:128);
i4(129:192, 129:192) = i3(65:128, 129:192);
i4(129:192, 65:128) = i3(129:192, 129:192);

figure; image(i1); colormap(gray(256)); axis image;
figure; image(i2); colormap(gray(256)); axis image;
figure; image(i3); colormap(gray(256)); axis image;
figure; image(i4); colormap(gray(256)); axis image;

%

[pyr1 p_sizes] = buildSFpyr(double(i1), 5, 3);
[pyr2 p_sizes] = buildSFpyr(double(i2), 5, 3);
[pyr3 p_sizes] = buildSFpyr(double(i3), 5, 3);
[pyr4 p_sizes] = buildSFpyr(double(i4), 5, 3);

pyramid1 = mb_change_pyramid_form(pyr1, p_sizes);
pyramid2 = mb_change_pyramid_form(pyr2, p_sizes);
pyramid3 = mb_change_pyramid_form(pyr3, p_sizes);
pyramid4 = mb_change_pyramid_form(pyr4, p_sizes);

%
pyramid_ori = pyramid1;
pyramid_ori(2:5, 2) = pyramid2(2:5, 2);
pyramid_ori(2:5, 3) = pyramid3(2:5, 3);
pyramid_ori(2:5, 4) = pyramid4(2:5, 4);

pyramid_scale = pyramid1;
pyramid_scale(3, 1:4) = pyramid2(3, 1:4);
pyramid_scale(4, 1:4) = pyramid3(4, 1:4);
pyramid_scale(5, 1:4) = pyramid4(5, 1:4);

i_ori = reconSFpyr(mb_change_pyramid_form(pyramid_ori), p_sizes);
i_scale = reconSFpyr(mb_change_pyramid_form(pyramid_scale), p_sizes);

figure; image(i_ori); colormap(gray(256)); axis image;
figure; image(i_scale); colormap(gray(256)); axis image;
%%

imwrite(uint8(i1), 'C:\isbe\dev\background\report\figures\lines1.bmp');
imwrite(uint8(i2), 'C:\isbe\dev\background\report\figures\lines2.bmp');
imwrite(uint8(i3), 'C:\isbe\dev\background\report\figures\lines3.bmp');
imwrite(uint8(i4), 'C:\isbe\dev\background\report\figures\lines4.bmp');

imwrite(uint8(i_ori), 'C:\isbe\dev\background\report\figures\lines_ori.bmp');
imwrite(uint8(i_scale), 'C:\isbe\dev\background\report\figures\lines_scale.bmp');
%%
for ii = 2:5
    for jj = 1:4
        figure; imagesc(pyramid1{ii, jj}); axis image; colormap(jet(256));
    end
end

%%

i_blank = i1;
i_blank(65:128, 65:192) = i1(1:64, 65:192);
i_blank(129:192, 65:192) = i1(193:256, 65:192);
figure; imagesc(i_blank); colormap(gray(256)); axis image;

%
[pyr_blank p_sizes] = buildSFpyr(double(i_blank), 5, 3);
pyramid_blank = mb_change_pyramid_form(pyr_blank, p_sizes);
%
pyramido1 = pyramid_blank;
pyramido1(2:5, 1) = pyramid1(2:5, 1);

pyramido2 = pyramid_blank;
pyramido2(2:5, 2) = pyramid1(2:5, 2);

pyramido3 = pyramid_blank;
pyramido3(2:5, 3) = pyramid1(2:5, 3);

pyramido4 = pyramid_blank;
pyramido4(2:5, 4) = pyramid1(2:5, 4);

i_o1 = reconSFpyr(mb_change_pyramid_form(pyramido1), p_sizes);
i_o2 = reconSFpyr(mb_change_pyramid_form(pyramido2), p_sizes);
i_o3 = reconSFpyr(mb_change_pyramid_form(pyramido3), p_sizes);
i_o4 = reconSFpyr(mb_change_pyramid_form(pyramido4), p_sizes);

figure; image(i_o1); colormap(gray(256)); axis image;
figure; image(i_o2); colormap(gray(256)); axis image;
figure; image(i_o3); colormap(gray(256)); axis image;
figure; image(i_o4); colormap(gray(256)); axis image;
%%
imwrite(uint8(i_o1), 'C:\isbe\dev\background\report\figures\lines_o1.bmp');
imwrite(uint8(i_o2), 'C:\isbe\dev\background\report\figures\lines_o2.bmp');
imwrite(uint8(i_o3), 'C:\isbe\dev\background\report\figures\lines_o3.bmp');
imwrite(uint8(i_o4), 'C:\isbe\dev\background\report\figures\lines_o4.bmp');

%%
pyramids1 = pyramid_blank;
pyramids1(2, 1:4) = pyramid1(2, 1:4);

pyramids2 = pyramid_blank;
pyramids2(3, 1:4) = pyramid1(3, 1:4);

pyramids3 = pyramid_blank;
pyramids3(4, 1:4) = pyramid1(4, 1:4);

pyramids4 = pyramid_blank;
pyramids4(5, 1:4) = pyramid1(5, 1:4);

i_s1 = reconSFpyr(mb_change_pyramid_form(pyramids1), p_sizes);
i_s2 = reconSFpyr(mb_change_pyramid_form(pyramids2), p_sizes);
i_s3 = reconSFpyr(mb_change_pyramid_form(pyramids3), p_sizes);
i_s4 = reconSFpyr(mb_change_pyramid_form(pyramids4), p_sizes);

figure; image(i_s1); colormap(gray(256)); axis image;
figure; image(i_s2); colormap(gray(256)); axis image;
figure; image(i_s3); colormap(gray(256)); axis image;
figure; image(i_s4); colormap(gray(256)); axis image;
%%
imwrite(uint8(i_s1), 'C:\isbe\dev\background\report\figures\lines_s1.bmp');
imwrite(uint8(i_s2), 'C:\isbe\dev\background\report\figures\lines_s2.bmp');
imwrite(uint8(i_s3), 'C:\isbe\dev\background\report\figures\lines_s3.bmp');
imwrite(uint8(i_s4), 'C:\isbe\dev\background\report\figures\lines_s4.bmp');
%%
real_line = double(rgb2gray(imread('C:\isbe\dev\background\images\normal_2\normal016.bmp')));
[pyr_real p_sizes] = buildSFpyr(double(real_line), 5, 4);
pyramid_real = mb_change_pyramid_form(pyr_real, p_sizes);
imwrite(uint8(real_line), 'C:\isbe\dev\background\report\figures\real_line.bmp');

real_plain = double(rgb2gray(imread('C:\isbe\dev\background\images\normal_2\normal017.bmp')));
[pyr_plain p_sizes] = buildSFpyr(double(real_plain), 5, 4);
pyramid_plain = mb_change_pyramid_form(pyr_plain, p_sizes);
imwrite(uint8(real_plain), 'C:\isbe\dev\background\report\figures\real_plain.bmp');
%%

for ii = 2:2
    for jj = 1:5
        figure;
        cmin = min(pyramid_real{ii, jj}(:));
        cmax = max(pyramid_real{ii, jj}(:));
        subplot(1,2,1); imagesc(pyramid_real{ii, jj}); axis image; colormap(jet(256)); caxis([cmin cmax]);
        subplot(1,2,2); imagesc(new_pyr{ii, jj}); axis image; colormap(jet(256)); caxis([cmin cmax]);
    end
end

%%
roi = roipoly(uint8(real_line));
roi = find(roi);
%%
new_pyr = mb_swap_pyramid_region(pyramid_plain, pyramid_real, roi, [], [2 5], [1 5]);
test_line = reconSFpyr(mb_change_pyramid_form(new_pyr), p_sizes);
figure; imagesc(test_line); axis image; colormap(gray(256));
imwrite(uint8(test_line), 'C:\isbe\dev\background\report\figures\real_line_1to5.bmp');

new_pyr = mb_swap_pyramid_region(pyramid_plain, pyramid_real, roi, [], [2 5], [1 1]);
test_line = reconSFpyr(mb_change_pyramid_form(new_pyr), p_sizes);
figure; imagesc(test_line); axis image; colormap(gray(256));
imwrite(uint8(test_line), 'C:\isbe\dev\background\report\figures\real_line_1.bmp');

new_pyr = mb_swap_pyramid_region(pyramid_plain, pyramid_real, roi, [], [2 5], [2 2]);
test_line = reconSFpyr(mb_change_pyramid_form(new_pyr), p_sizes);
figure; imagesc(test_line); axis image; colormap(gray(256));
imwrite(uint8(test_line), 'C:\isbe\dev\background\report\figures\real_line_2.bmp');

new_pyr = mb_swap_pyramid_region(pyramid_plain, pyramid_real, roi, [], [2 5], [3 3]);
test_line = reconSFpyr(mb_change_pyramid_form(new_pyr), p_sizes);
figure; imagesc(test_line); axis image; colormap(gray(256));
imwrite(uint8(test_line), 'C:\isbe\dev\background\report\figures\real_line_3.bmp');

new_pyr = mb_swap_pyramid_region(pyramid_plain, pyramid_real, roi, [], [2 5], [4 4]);
test_line = reconSFpyr(mb_change_pyramid_form(new_pyr), p_sizes);
figure; imagesc(test_line); axis image; colormap(gray(256));
imwrite(uint8(test_line), 'C:\isbe\dev\background\report\figures\real_line_4.bmp');

new_pyr = mb_swap_pyramid_region(pyramid_plain, pyramid_real, roi, [], [2 5], [5 5]);
test_line = reconSFpyr(mb_change_pyramid_form(new_pyr), p_sizes);
figure; imagesc(test_line); axis image; colormap(gray(256));
imwrite(uint8(test_line), 'C:\isbe\dev\background\report\figures\real_line_5.bmp');
