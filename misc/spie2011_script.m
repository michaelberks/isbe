%--------------------------------------------------------------------------
%% Load in the images we have saved

% load in the 1st retinogram
ret = imread('G:\Retina Drive\DRIVE\test\images\01_test.tif');
ret = rgb2gray(ret);

%load in the vessel map
vessel_map = u_load('G:\1000Trees-In-paintingmethod4m4m5\Votes\V01.mat');
vessel_map = reshape(vessel_map(:,1),584,566);

% load in vessel centreline map for image
centre_map = u_load('G:\1000TreesWithin-painting-centerlines\votes\V01.mat');
centre_map = reshape(centre_map(:,1),584,566);

%load in ground truth
gt = logical(imread('G:\Retina Drive\DRIVE\test\1st_manual\01_manual1.gif'));
gt(:,566) = 0;

%skelentonise ground truth
gts = bwmorph(gt ,'skel',Inf);

%load in the out-of-image mask
mask = logical(imread('G:\mask\01_test_mask.gif'));
mask(:,566) = 0;
%%
%--------------------------------------------------------------------------
% %load in the regression forest for orientation
% load 'G:\Regression\Results\rf_forest.mat';
% 
% %Use regression forest to compute orientation of retinogram and save
% [orientation_map] = classify_image(...
%         'image_in', image_in,...
%         'forest', random_forest,...
%         'forest_type', 'regression',...
%         'num_levels', 6,...
%         'win_size', 1);
% orientation_map(:,566)=orientation_map(:,565);   
% save G:\orientation_map01.mat orientation_map
%%
%load orientation_map
******** LOAD ORIENTATION MAP HERE ********************
%%

%Use the orientation map to NMS the vessel centreline map
smooth_centre = imfilter(centre_map, fspecial('Gaussian'));
centre_nms = mb_non_maximal_supp(smooth_centre, orientation_map);
%%
% Method 1: threshold vessel map
bw_method1 = vessel_map > 700 & mask;

% Method 2: threshold vessel centreline map
bw_method2 = centre_map > 700 & mask;

% Method 1: threshold vessel centreline map
bw_method3 = centre_map > 700 & mask;

% Method 4: threshold NMS and dilate
bw_method4 = imdilate(centre_nms> 700, strel('disk', 1))  & mask;

% figure; imagesc(bw_method1); axis image;
% figure; imagesc(bw_method2); axis image;
% figure; imagesc(bw_method3); axis image;
% figure; imagesc(bw_method4); axis image;
%
%Creat output images of TP and FP for each method

% Method 1: TP = GT & V; FP = ~GT & V
[tp_fp_image1] = make_tp_fp_image(gt, ~gt, bw_method1);

% Method 2: TP = GTS & C; FP = ~GTS & C
[tp_fp_image2] = make_tp_fp_image(gts, ~gts, bw_method2);

% Method 3: TP = GTS & C; FP = ~GTS & C
[tp_fp_image3] = make_tp_fp_image(gts, ~gt, bw_method3);

% Method 4: TP = GTS & C; FP = ~GTS & C
[tp_fp_image4] = make_tp_fp_image(gts, ~gt, bw_method4);

figure; imagesc(tp_fp_image1); axis image; a1 = gca;
figure; imagesc(tp_fp_image2); axis image; a2 = gca;
figure; imagesc(tp_fp_image3); axis image; a3 = gca;
figure; imagesc(tp_fp_image4); axis image; a4 = gca;
linkaxes([a1 a2 a3 a4]);
%
imwrite(tp_fp_image1, 'G:\tp_fp_image1.bmp');
imwrite(tp_fp_image2, 'G:\tp_fp_image2.bmp');
imwrite(tp_fp_image3, 'G:\tp_fp_image3.bmp');
imwrite(tp_fp_image4, 'G:\tp_fp_image4.bmp');

imwrite(tp_fp_image1(230:429,230:429,:), 'G:\tp_fp_image1_roi.bmp');
imwrite(tp_fp_image2(230:429,230:429,:), 'G:\tp_fp_image2_roi.bmp');
imwrite(tp_fp_image3(230:429,230:429,:), 'G:\tp_fp_image3_roi.bmp');
imwrite(tp_fp_image4(230:429,230:429,:), 'G:\tp_fp_image4_roi.bmp');
%%
