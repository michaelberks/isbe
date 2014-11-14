i1 = imread('C:\isbe\dev\background\images\normal512\o04_001LCC_1024_3427_865.bmp');

args.image1 = double(i1(1:128,1:128));
args.image2 = double(i1(129:256,129:256));
args.num_train = 1000;
args.d = 20;
args.tree_dir = 'C:\isbe\dev\misc';
args.do_test1 = 1;
args.do_test2 = 2;
args.max_test_size = 64;
args.n_trees = 2;
args.mask1 = [false(128,28) true(128,100)];
args.mask2 = [false(28,128); true(100,128)];
args.n_trees = 2;
clear i1;
%%
[random_forest] =...
    mb_random_forest_two_images(args);
%%
i1 = double(imresize(imread('C:\isbe\mammograms\new_CAD\bMP_2004_Normals\024LML.bmp'), [1024 NaN], 'bilinear'));
i2 = double(imresize(imread('C:\isbe\mammograms\new_CAD\bMP_2004_Normals\024RML.bmp'), [1024 NaN], 'bilinear'));
s1 = u_load('C:\isbe\mammograms\new_CAD\BMP_2004_Normals\segmentations\024LML_segmentation.mat');
s2 = u_load('C:\isbe\mammograms\new_CAD\BMP_2004_Normals\segmentations\024RML_segmentation.mat');
mask1 = roipoly(i1, s1.breast_border(:,1), s1.breast_border(:,2));
mask2 = roipoly(i2, s2.breast_border(:,1), s2.breast_border(:,2));
%
args.image1 = i1;
args.image2 = fliplr(i2);
args.num_train = 10000;
args.d = 5;
args.tree_dir = 'C:\isbe\dev\misc';
args.do_test1 = 1;
args.do_test2 = 0;
args.do_max = 0;
args.win_size = 1;
args.max_test_size = 64;
args.n_trees = 2;
args.mask1 = mask1;
args.mask2 = fliplr(mask2);

clear i1 i2 s1 s2 mask1 mask2;
%
args.n_trees = 100;
[random_forest] =...
    mb_random_forest_two_images(args);
save M:\asymmetry_project\results\contralateral_data\rfs\rf_024ML.mat random_forest
%%
i1 = double(imresize(u_load('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\o04_024LML.mat'), [1024 NaN], 'bilinear'));
i2 = double(imresize(u_load('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\o04_024RML.mat'), [1024 NaN], 'bilinear'));
s1 = u_load('C:\isbe\asymmetry_project\data\segmentations\2004_screening\abnormals\o04_024LML_segmentation.mat');
s2 = u_load('C:\isbe\asymmetry_project\data\segmentations\2004_screening\abnormals\o04_024RML_segmentation.mat');
mask1 = roipoly(i1, s1.breast_border(:,1), s1.breast_border(:,2));
mask2 = roipoly(i2, s2.breast_border(:,1), s2.breast_border(:,2));

% sub_mask1 = false(size(mask1));
% sub_mask1(2:2:end, 2:2:end) = true;
% sub_mask2 = false(size(mask2));
% sub_mask2(2:2:end, 2:2:end) = true;
%
args.image2 = i1;
args.image1 = fliplr(i2);
args.num_train = 10000;
args.num_levels = 3;
args.d = 5;
args.tree_dir = 'C:\isbe\asymmetry_project\data\misc';
args.do_test1 = 1;
args.do_test2 = 2;
args.use_probs = 0;
args.do_max = 0;
args.win_size = 1;
args.max_test_size = 64;
% args.train_mask2 = mask1;
% args.train_mask1 = fliplr(mask2);
% args.test_mask2 = mask1 & sub_mask1;
% args.test_mask1 = fliplr(mask2 & sub_mask2);
args.mask2 = mask1;
args.mask1 = fliplr(mask2);

args.n_trees = 10;
args.save_path = ['C:\isbe\asymmetry_project\data\misc\rf_024ML_W' num2str(args.win_size) 'L' num2str(args.num_levels) 'P.mat'];
clear i1 i2 s1 s2 mask1 mask2 sub_mask1 sub_mask2;
%%
profile on    
[random_forest] =...
    mb_random_forest_two_images(args);
profile off
profsave(profile('info'),'C:\isbe\asymmetry_project\data\contralateral_rfs\misc\profile_results_orig')
clear random_forest;

profile on    
[random_forest] =...
    mb_random_forest_two_images_a(args);
profile off
profsave(profile('info'),'C:\isbe\asymmetry_project\data\contralateral_rfs\misc\profile_results_nag')
%%
im1_votes1 = random_forest.image1_votes1 ./ random_forest.image1_total_votes;
im1_votes1(~random_forest.image1_total_votes) = 0;

im2_votes1 = random_forest.image2_votes1 ./ random_forest.image2_total_votes;
im2_votes1(~random_forest.image2_total_votes) = 0;

figure; imagesc(im1_votes1); axis image;
figure; imagesc(im2_votes1); axis image;
figure; imagesc(args.image1); axis image; colormap(gray(256));
figure; imagesc(args.image2); axis image; colormap(gray(256));
%%
write_im_from_colormap(im1_votes1, 'M:\asymmetry_project\data\contralateral_rfs\2004_screening\figures\024RML_cancer_prob_map_W1L3P.jpg', jet(256));
write_im_from_colormap(im2_votes1, 'M:\asymmetry_project\data\contralateral_rfs\2004_screening\figures\024LML_cancer_prob_map_W1L3P.jpg', jet(256));
%write_im_from_colormap(args.image1, 'M:\asymmetry_project\data\contralateral_data\rfs\figures\024RML_cancer.jpg', jet(256));
%write_im_from_colormap(args.image2, 'M:\asymmetry_project\data\contralateral_data\rfs\figures\024LML_cancer.jpg', jet(256));

f1 = figure; hist(im1_votes1(random_forest.image1_total_votes>0), 101);
saveas(f1, 'M:\asymmetry_project\data\contralateral_rfs\2004_screening\figures\024RML_cancer_prob_hist_W1L3P.jpg');

f2 = figure; hist(im2_votes1(random_forest.image1_total_votes>0), 101);
saveas(f2, 'M:\asymmetry_project\data\contralateral_rfs\2004_screening\figures\02LRML_cancer_prob_hist_W1L3P.jpg');