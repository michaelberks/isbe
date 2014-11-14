%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script producing figures for the IWDM 2008 mass model/pyramid paper

for ii = 2:4
    for jj = 1:5
        cmin = min(orig_pyr{ii,jj}(:));
        cmax = max(orig_pyr{ii,jj}(:));
        figure; 
        subplot(1,3,1); imagesc(orig_pyr{ii, jj}); axis image; colormap(jet(256)); caxis([cmin cmax]);
        subplot(1,3,2); imagesc(pyramid{ii, jj}); axis image; colormap(jet(256)); caxis([cmin cmax]);
        subplot(1,3,3); imagesc(pyramid{ii, jj} - orig_pyr{ii, jj}); axis image; colormap(jet(256));
    end
end
%%
for ii = 2:4
    for jj = 1:5
        cmin = min(orig_pyr{ii,jj}(:));
        cmax = max(orig_pyr{ii,jj}(:));
        figure; 
        subplot(1,2,1); imagesc(mc_pyr{ii, jj}); axis image; colormap(jet(256)); caxis([cmin cmax]);
        subplot(1,2,2); imagesc(pyramid{ii, jj}); axis image; colormap(jet(256)); caxis([cmin cmax]);
    end
end

%%
%orig_im = imread('C:/isbe/dev/background/images/normal_2/normal006.bmp');
orig_im = uint8(pyramid);
[rows, cols] = size(orig_im);

phi = linspace(0, 2*pi, 200);
x = 100*cos(phi) + round(cols/2);
y = 100*sin(phi) + round(rows/2);

xx = ceil(x);
yy = ceil(y);
clear phi x y
%contrast enhance image
orig_im = double(orig_im);
orig_im = (orig_im - min(orig_im(:))) / (max(orig_im(:)) - min(orig_im(:)));

%scale to 0-255 and round
orig_im = uint8(255*orig_im);

mask = zeros(size(orig_im));
mask(sub2ind(size(orig_im), yy, xx)) = 1;

clear anno*
anno_r = orig_im;
anno_r(logical(mask)) = 255;
anno_gb = orig_im;
anno_gb(logical(mask)) = 0;

anno_rgb(:,:,1) = anno_r;
anno_rgb(:,:,2) = anno_gb;
anno_rgb(:,:,3) = anno_gb;

figure; image(anno_rgb); axis image;
imwrite(anno_rgb, 'K:\isbe\iwdm2008\figures\norm_006_syn.bmp');
%%
%clear
orig_im = imread('K:\isbe\iwdm2008\figures\norm_006_3_1_syn.bmp');
[rows, cols, dummy] = size(orig_im);

phi = linspace(0, 2*pi, 100);
x = 50*cos(phi) + round(cols/2);
y = 50*sin(phi) + round(rows/2);

xx = ceil(x);
yy = ceil(y);
clear phi x y

mask = zeros(rows, cols);
mask(sub2ind([rows cols], yy, xx)) = 1;
%
anno_r = orig_im(:,:,1);
anno_r(logical(mask)) = 255;
anno_g = orig_im(:,:,2);
anno_g(logical(mask)) = 0;
anno_b = orig_im(:,:,3);
anno_b(logical(mask)) = 0;

anno_rgb(:,:,1) = anno_r;
anno_rgb(:,:,2) = anno_g;
anno_rgb(:,:,3) = anno_b;

figure; image(anno_rgb); axis image;
imwrite(anno_rgb, 'K:\isbe\iwdm2008\figures\norm_006_3_1_syn.bmp');
%%
%orig_im = imread('C:/isbe/dev/background/images/pectoral_2/pectoral001.bmp');
%orig_pyr = u_load('C:/isbe/dev/background/pyramid/pectoral_2/pectoral001_pyramid.mat');
% make mask for circle
%
r = 100;
[rows cols] = size(orig_im);
c_x = round(cols/2); c_y = round(rows/2);
[x y] = meshgrid(1:cols, 1:rows);
circ = (x - c_x).^2 + (y - c_y).^2 > r.^2;
clear r c_x c_y x y
%
%set up synthesis args
syn_args.FilledImage = circ;   
%make sure uint8 mode is turned off since we're synthesising doubles
syn_args.Uint8Mode = 1;

syn_args.PathToTextureGMM = 'C:\isbe\dev\background\results\pectoral_model';
syn_args.SeededImage = orig_im;
[syn_image cluster_image] = mb_gmm_tex_synthesis(syn_args);
save C:\isbe\dev\background\syn\pectoral001_circ100 syn_image cluster_image
imwrite(syn_image, 'K:\isbe\iwdm2008\figures\normal006_syn2_1.bmp');
%%
clear
orig_pyr = u_load('C:/isbe/dev/background/pyramid/normal_2/normal006.bmp_pyramid.mat');
r = 100;
[rows cols] = size(orig_pyr{1,1});
c_x = round(cols/2); c_y = round(rows/2);
[x y] = meshgrid(1:cols, 1:rows);
circ = (x - c_x).^2 + (y - c_y).^2 > r.^2;
clear r c_x c_y x y

syn_args.SamplePyramid = orig_pyr;
syn_args.FilledImage = circ;
syn_args.ModelDir = 'C:\isbe\dev\background\results\normal_2\level_5';
syn_args.ModelName = 'normal_2_model_a';
syn_args.CutOffLevel = 5;
syn_args.Plot = 0;

syn_args.SaveFile = [mberksroot, 'background\syn\normal006_normal_2_model_a'];
[synthesised_image, new_pyr] = mb_gmm_pyr_synthesis(syn_args);
save K:\isbe\iwdm2008\normal006_circ100_pyr synthesised_image new_pyr
%%
orig_pyr = u_load('C:/isbe/dev/background/pyramid/normal_2/normal006.bmp_pyramid.mat');
for ii = 2:4
    for jj = 1:5
        write_im_from_colormap(orig_pyr{ii, jj},...
            ['C:\isbe\dev\background\misc\pyr_ims\normal006\normal006_',...
            num2str(ii), '_', num2str(jj), '.bmp'], colormap(jet(256)));
    end
end

%%
%look at some synthesised textures:
for ii = 1:15
    load(['C:\isbe\dev\background\syn\normal_2_2levels\conditioned_synthesis_circle_6_', zerostr(ii,3)]);
    orig_image = double(imread(['C:\isbe\dev\background\images\normal_2\normal', zerostr(ii,3), '.bmp']));
    figure; 
    subplot(1,2,1); imagesc(synthesised_image); axis image; colormap(gray(256));
    subplot(1,2,2); imagesc(orig_image); axis image; colormap(gray(256));
    clear synthesised_image pyramid cluster_image
end
%%
%look at some synthesised textures:
for ii = 1:15
    load(['C:\isbe\dev\background\syn\mass_2_2levels\conditioned_synthesis_circle_6_', zerostr(ii,3)]);
    orig_image = double(imread(['C:\isbe\dev\background\images\normal_2\normal', zerostr(ii,3), '.bmp']));
    figure; 
    subplot(1,2,1); image(synthesised_image); axis image; colormap(gray(256));
    subplot(1,2,2); image(orig_image); axis image; colormap(gray(256));
    clear synthesised_image pyramid cluster_image
end
%%
% Of these I like 6, 7 and 11 best, so save them, and the difference maps
% to figures;
for ii = [6 7 11]
    load(['C:\isbe\dev\background\syn\mass_2_2levels\conditioned_synthesis_circle_6_', zerostr(ii,3)]);
    orig_image = double(imread(['C:\isbe\dev\background\images\normal_2\normal', zerostr(ii,3), '.bmp']));
    %write_im_from_colormap(uint8(synthesised_image), ['C:\isbe\iwdm2008\figures\syn_normal', zerostr(ii,3), '.bmp'], gray(256), [0 255]);
    %write_im_from_colormap(uint8(orig_image), ['C:\isbe\iwdm2008\figures\syn_normal_orig_', zerostr(ii,3), '.bmp'], gray(256), [0 255]);
    write_im_from_colormap(synthesised_image-orig_image, ['C:\isbe\iwdm2008\figures\syn_normal_dif_', zerostr(ii,3), '.bmp'], jet(256), [-20 20]);
    
    clear synthesised_image pyramid cluster_image orig_im
end
%%
%also print some pyramid levels
load('C:\isbe\dev\background\syn\mass_2_2levels\conditioned_synthesis_circle_6_006');
pyr_orig = u_load('C:\isbe\dev\background\pyramid\normal_2\normal006.bmp_pyramid.mat');
%%
figure; imagesc(pyramid{3,4}); axis image; caxis([-15 15]);
figure; imagesc(pyr_orig{3,4}); axis image; caxis([-15 15]);
figure; imagesc(pyr_orig{3,4} - pyramid{3,4}); axis image; caxis([-15 15]);
%%
write_im_from_colormap(pyramid{3,4}, 'C:\isbe\iwdm2008\figures\syn_normal006_3_4.bmp', jet(256), [-15 15]);
write_im_from_colormap(pyr_orig{3,4}, 'C:\isbe\iwdm2008\figures\syn_normal006_3_4_real.bmp', jet(256), [-15 15]);
write_im_from_colormap(pyr_orig{3,4} - pyramid{3,4}, 'C:\isbe\iwdm2008\figures\syn_normal006_3_4_dif.bmp', jet(256), [-15 15]);
%%
%look at some synthesised textures:
for ii = 10:15
    load(['C:\isbe\dev\background\syn\mass_2_2levels\conditioned_synthesis_circle', zerostr(ii,3)]);
    orig_image = double(imread(['C:\isbe\dev\background\images\normal_2\normal', zerostr(ii,3), '.bmp']));
    [p pp] = mb_change_pyramid_form(pyramid);
    synthesised_image = reconSFpyr(p, pp);
    save(['C:\isbe\dev\background\syn\mass_2_2levels\conditioned_synthesis_circle', zerostr(ii,3)], 'synthesised_image', 'pyramid', 'cluster_image');
    figure; 
    subplot(1,2,1); image(synthesised_image); axis image; colormap(gray(256));
    subplot(1,2,2); image(orig_image); axis image; colormap(gray(256));
    clear synthesised_image pyramid cluster_image
end
%%
clear; pack;
%%
pyr_orig = u_load('C:\isbe\dev\background\pyramid\normal_2\normal008.bmp_pyramid.mat');
%%
[rows cols] = size(pyr_orig{2,1});
row_centre = round(rows / 2);
col_centre = round(cols / 2);

% Make a the biggest circular mask
m = min([rows cols]);
rad = floor((m - 128) / 2);
[x y] = meshgrid(1:cols, 1:rows);
filled_image = (x - col_centre).^2 + (y - row_centre).^2 > rad.^2;
temp = logical(ones(size(pyr_orig{2,1})));
temp(1:row_centre, 1:col_centre) = 0;
filled_image = filled_image | temp;
clear x y m rad row_centre col_centre temp;

%%
%Draw diagram of synthesis
pyramid = u_load('C:\isbe\dev\background\pyramid\normal_2\normal016_pyramid.mat');

figure; imagesc(pyramid{3, 4}); axis image; colormap(jet(256));
figure; imagesc(pyramid{4, 4}); axis image; colormap(jet(256));
%%
%nice area centred at r=80 c= 100 in {3,4} - take windows from each level
win1 = sample_window(pyramid{3, 4}, 32, 80, 100);
win2 = sample_window(pyramid{4, 4}, 16, 40, 50);

figure; imagesc(win1); axis image; colormap(jet(256)); caxis([-5 5]);
figure; imagesc(win2); axis image; colormap(jet(256)); caxis([-5 5]);
%%
figure; mask1 = roipoly(win1);
figure; imagesc(mask1); colormap(gray); axis image;
win1(mask1) = min(win1(~mask1)) - 1;
%
% mask2 = imresize(mask1, 0.5);
% figure; imagesc(mask2); colormap(gray); axis image;
% win2(mask2) = min(win2(~mask2)) - 1;

my_map = [1 1 1; jet(256)];
my_map(end,:) = [];
%%
figure; imagesc(win1); axis image; colormap(my_map);
figure; imagesc(win2); axis image; colormap(jet);
figure; imagesc(sample_window(pyramid{5, 4}, 16, 20, 25)); axis image; colormap(jet); caxis([-4 4]);
draw_grid(8,8,[4.5 4.5], 8, 'k')
draw_grid(5,5,[6.5 6.5], 1, 'k')

figure; imagesc(win2); axis image; colormap(jet); caxis(clims);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Code for mass model paper
%%
% Look at mass + mass regenerations to pick best image
load('C:\isbe\dev\files\u_files.mat');
%%
for ii = 81:101
    load(['C:\isbe\dev\masses\', u_files1(ii).name]);
    load(['C:\isbe\dev\recon\recon_', zerostr(ii, 3), '.mat']);
    figure;
    subplot(1,3,1); imagesc(mass.subtract_ROI); axis image; colormap(gray(256));
    subplot(1,3,2); imagesc(old_shape_ROI); axis image; colormap(gray(256));
    subplot(1,3,3); imagesc(new_shape_ROI); axis image; colormap(gray(256));
    clear mass *ROI
end
%%
gd = [11 28 45 51 56 61 89];
for ii = 1:length(gd)
    load(['C:\isbe\dev\masses\', u_files1(gd(ii)).name]);
    load(['C:\isbe\dev\recon\recon_', zerostr(gd(ii), 3), '.mat']);
    figure;
    subplot(1,3,1); imagesc(mass.subtract_ROI); axis image; colormap(gray(256));
    subplot(1,3,2); imagesc(old_shape_ROI); axis image; colormap(gray(256));
    subplot(1,3,3); imagesc(new_shape_ROI); axis image; colormap(gray(256));
    clear mass *ROI
end
%%
gd = [11 28 45 51 56 61 89];
for ii = 1:length(gd)
    load(['C:\isbe\dev\masses\', u_files1(gd(ii)).name]);
    load(['C:\isbe\dev\recon\recon_', zerostr(gd(ii), 3), '.mat']);
    write_im_from_colormap([mass.subtract_ROI old_shape_ROI new_shape_ROI], ...
        ['K:\isbe\iwdm2008\figures\recon_', zerostr(gd(ii), 3), '.bmp']);
    clear mass *ROI
end
