%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make like an Egyptation and build me a pyramid
mb_build_pyramids('ImageDir', 'C:\isbe\dev\background\images\mc_2\',...
                  'OutputDir', 'C:\isbe\dev\background\pyramid\mc_2\',...
                  'NumLevels', 5,...
                  'NumOrientations', 5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to swap 3rd/4th levels of pyrmaids between 2 regions - It's good fun
% this!!

%%
%load images
im5 = imread('C:\isbe\dev\background\images\normal_2\normal010.bmp');
im10 = imread('C:\isbe\dev\background\images\mc_2\mc010.bmp');
%

% make mask for circles
r = 100;
c_x = 180; c_y = 180;

[x5 y5] = meshgrid(1:size(im5, 1), 1:size(im5, 2));
[x10 y10] = meshgrid(1:size(im10, 1), 1:size(im10, 2));

circ5 = (x5 - c_x).^2 + (y5 - c_y).^2 < r.^2;
circ10 = (x10 - c_x).^2 + (y10 - c_y).^2 < r.^2;
%
% Construct pyramids
[pyramid5 p_sizes5] = buildSFpyr(double(im5), 5, 4);
[pyramid10 p_sizes10] = buildSFpyr(double(im10), 5, 4);

% Change pyramid form
pyramid5 = mb_change_pyramid_form(pyramid5, p_sizes5);
pyramid10 = mb_change_pyramid_form(pyramid10, p_sizes10);

% Compute indices to swap at levels 3 and 4;
% idx5 = find(circ5); %level 2 indices
% idx10 = find(circ10);

[r5 c5] = find(circ5);
[r10 c10] = find(circ10);

idx5_2 = sub2ind(size(pyramid5{2,1}), r5, c5);
idx10_2 = sub2ind(size(pyramid10{2,1}), r10, c10);

r5_3 = ceil(r5/2); c5_3 = ceil(c5/2);
r10_3 = ceil(r10/2); c10_3 = ceil(c10/2);

idx5_3 = unique(sub2ind(p_sizes5(7,:), r5_3, c5_3)); %level 3 indices
idx10_3 = unique(sub2ind(p_sizes10(7,:), r10_3, c10_3));

r5_4 = ceil(r5/4); c5_4 = ceil(c5/4);
r10_4 = ceil(r10/4); c10_4 = ceil(c10/4);

idx5_4 = unique(sub2ind(p_sizes5(12,:), r5_4, c5_4)); %level 4 indices
idx10_4 = unique(sub2ind(p_sizes10(12,:), r10_4, c10_4));

%
%swap pyramid values at given indices
pyramid5a = pyramid5;
pyramid10a = pyramid10;

%change lev 4
for ori = 1:5
    pyramid5a{4,ori}(idx5_4) = pyramid10{4,ori}(idx10_4);
    pyramid10a{4,ori}(idx10_4) = pyramid5{4,ori}(idx5_4);
    
    %plot results
    figure;
    subplot(2,2,1); imagesc(pyramid5{4,ori}); colormap(gray(256)); axis image;
    subplot(2,2,2); imagesc(pyramid5a{4,ori}); colormap(gray(256)); axis image;
    subplot(2,2,3); imagesc(pyramid10{4,ori}); colormap(gray(256)); axis image;
    subplot(2,2,4); imagesc(pyramid10a{4,ori}); colormap(gray(256)); axis image;
    
end

%change lev 3
for ori = 1:5
    pyramid5a{3,ori}(idx5_3) = pyramid10{3,ori}(idx10_3);
    pyramid10a{3,ori}(idx10_3) = pyramid5{3,ori}(idx5_3);
    
    %plot results
    figure;
    subplot(2,2,1); imagesc(pyramid5{3,ori}); colormap(gray(256)); axis image;
    subplot(2,2,2); imagesc(pyramid5a{3,ori}); colormap(gray(256)); axis image;
    subplot(2,2,3); imagesc(pyramid10{3,ori}); colormap(gray(256)); axis image;
    subplot(2,2,4); imagesc(pyramid10a{3,ori}); colormap(gray(256)); axis image;
end

%change lev 3
for ori = 1:5
    pyramid5a{2,ori}(idx5_2) = pyramid10{2,ori}(idx10_2);
    pyramid10a{2,ori}(idx10_2) = pyramid5{2,ori}(idx5_2);
    
    %plot results
    figure;
    subplot(2,2,1); imagesc(pyramid5{2,ori}); colormap(gray(256)); axis image;
    subplot(2,2,2); imagesc(pyramid5a{2,ori}); colormap(gray(256)); axis image;
    subplot(2,2,3); imagesc(pyramid10{2,ori}); colormap(gray(256)); axis image;
    subplot(2,2,4); imagesc(pyramid10a{2,ori}); colormap(gray(256)); axis image;
end

%
%reconstruct images
pyramid5a = mb_change_pyramid_form(pyramid5a);
pyramid10a = mb_change_pyramid_form(pyramid10a);

im5_r = reconSFpyr(pyramid5a, p_sizes5);
im10_r = reconSFpyr(pyramid10a, p_sizes10);
%
figure; image(im5); colormap(gray(256)); axis image;
figure; image(im5_r); colormap(gray(256)); axis image;
figure; image(im10); colormap(gray(256)); axis image;
figure; image(im10_r); colormap(gray(256)); axis image;
%%
save C:\isbe\dev\background\syn\n10_with_n5 im5_r im5 im10 im10_r;
%%
temp2 = zeros(size(pyramid5{2,1}));
temp2(idx5_2) = 1;
temp3 = zeros(size(pyramid5{3,1}));
temp3(idx5_3) = 1;
temp4 = zeros(size(pyramid5{4,1}));
temp4(idx5_4) = 1;
%
figure; imagesc(temp2); colormap(gray(256)); axis image;
figure; imagesc(temp3); colormap(gray(256)); axis image;
figure; imagesc(temp4); colormap(gray(256)); axis image;
%%
hist_args.PyramidDir = 'C:\isbe\dev\background\pyramid\normal_2';

hist_args.Level = 2;
hist_args.Orientation = 1;
hist_args.PyramidModel = u_load('C:\isbe\dev\background\results\pyramid\normal_2_model_2_1');
hist_args.Border = 20;
[clust_hist_2_1] = mb_find_cluster_support(hist_args);

save C:\isbe\dev\background\misc\clust_hist_2_1 clust_hist_2_1

ch = clust_hist_2_1;
model = u_load('C:/isbe/dev/background/results/pyramid/normal_2_model_2_1.mat');

pch = 100*ch / sum(ch);
scrap_idx = find(pch < 0.1);
display(['num cluster to scrap = ', num2str(length(scrap_idx))]);

model.Means(scrap_idx, :) = [];
model.CovMats(scrap_idx) = [];
model.NumClusters = model.NumClusters - length(scrap_idx);
model.ClusterProbs(scrap_idx) = [];
save C:/isbe/dev/background/results/pyramid/normal_2_model_a_2_1.mat model
%%
hist_args.PyramidDir = 'C:\isbe\dev\background\pyramid\normal_2';

hist_args.Level = 4;
hist_args.Orientation = 4;
hist_args.PyramidModel = u_load('C:\isbe\dev\background\results\pyramid\normal_2_model_4_4');
hist_args.Border = 20;
[clust_hist_4_4] = mb_find_cluster_support(hist_args);

save C:\isbe\dev\background\misc\clust_hist_4_4 clust_hist_4_4

ch = clust_hist_4_4;
model = u_load('C:/isbe/dev/background/results/pyramid/normal_2_model_4_4.mat');

pch = 100*ch / sum(ch);
scrap_idx = find(pch < 0.1);
display(['num cluster to scrap = ', num2str(length(scrap_idx))]);

model.Means(scrap_idx, :) = [];
model.CovMats(scrap_idx) = [];
model.NumClusters = model.NumClusters - length(scrap_idx);
model.ClusterProbs(scrap_idx) = [];
save C:/isbe/dev/background/results/pyramid/normal_2_model_a_4_4.mat model
%
hist_args.PyramidDir = 'C:\isbe\dev\background\pyramid\normal_2';

hist_args.Level = 4;
hist_args.Orientation = 5;
hist_args.PyramidModel = u_load('C:\isbe\dev\background\results\pyramid\normal_2_model_4_5');
hist_args.Border = 20;
[clust_hist_4_5] = mb_find_cluster_support(hist_args);

save C:\isbe\dev\background\misc\clust_hist_4_5 clust_hist_4_5

ch = clust_hist_4_5;
model = u_load('C:/isbe/dev/background/results/pyramid/normal_2_model_4_5.mat');

pch = 100*ch / sum(ch);
scrap_idx = find(pch < 0.1);
display(['num cluster to scrap = ', num2str(length(scrap_idx))]);

model.Means(scrap_idx, :) = [];
model.CovMats(scrap_idx) = [];
model.NumClusters = model.NumClusters - length(scrap_idx);
model.ClusterProbs(scrap_idx) = [];
save C:/isbe/dev/background/results/pyramid/normal_2_model_a_4_5.mat model
%%
hist_args.Level = 3;
hist_args.Orientation = 4;
hist_args.PyramidModel = u_load('C:\isbe\dev\background\results\pyramid\normal_2_model_3_4');
hist_args.Border = 20;
[sum_elements] = mb_find_cluster_support(hist_args);
clust_hist_3_4(end+1) = sum_elements - sum(clust_hist_3_4);
save C:\isbe\dev\background\misc\clust_hist_3_4 clust_hist_3_4
%
hist_args.Level = 3;
hist_args.Orientation = 5;
hist_args.PyramidModel = u_load('C:\isbe\dev\background\results\pyramid\normal_2_model_3_5');
hist_args.Border = 20;
[sum_elements] = mb_find_cluster_support(hist_args);
clust_hist_3_5(end+1) = sum_elements - sum(clust_hist_3_5);
save C:\isbe\dev\background\misc\clust_hist_3_5 clust_hist_3_5
%
hist_args.Level = 4;
hist_args.Orientation = 2;
hist_args.PyramidModel = u_load('C:\isbe\dev\background\results\pyramid\normal_2_model_4_2');
hist_args.Border = 20;
[sum_elements] = mb_find_cluster_support(hist_args);
clust_hist_4_2(end+1) = sum_elements - sum(clust_hist_4_2);
save C:\isbe\dev\background\misc\clust_hist_4_2 clust_hist_4_2
%
hist_args.Level = 4;
hist_args.Orientation = 3;
hist_args.PyramidModel = u_load('C:\isbe\dev\background\results\pyramid\normal_2_model_4_3');
hist_args.Border = 20;
[sum_elements] = mb_find_cluster_support(hist_args);
clust_hist_4_3(end+1) = sum_elements - sum(clust_hist_4_3);
save C:\isbe\dev\background\misc\clust_hist_4_3 clust_hist_4_3
%
hist_args.Level = 2;
hist_args.Orientation = 2;
hist_args.PyramidModel = u_load('C:\isbe\dev\background\results\pyramid\normal_2_model_2_2');
hist_args.Border = 20;
[sum_elements] = mb_find_cluster_support(hist_args);
clust_hist_2_2(end+1) = sum_elements - sum(clust_hist_2_2);
save C:\isbe\dev\background\misc\clust_hist_2_2 clust_hist_2_2
%
hist_args.Level = 2;
hist_args.Orientation = 4;
hist_args.PyramidModel = u_load('C:\isbe\dev\background\results\pyramid\normal_2_model_2_4');
hist_args.Border = 20;
[sum_elements] = mb_find_cluster_support(hist_args);
clust_hist_2_4(end+1) = sum_elements - sum(clust_hist_2_4);
save C:\isbe\dev\background\misc\clust_hist_2_4 clust_hist_2_4
%
hist_args.Level = 2;
hist_args.Orientation = 5;
hist_args.PyramidModel = u_load('C:\isbe\dev\background\results\pyramid\normal_2_model_2_5');
hist_args.Border = 20;
[sum_elements] = mb_find_cluster_support(hist_args);
clust_hist_2_5(end+1) = sum_elements - sum(clust_hist_2_5);
save C:\isbe\dev\background\misc\clust_hist_2_5 clust_hist_2_5
%
hist_args.Level = 5;
hist_args.Orientation = 1;
hist_args.PyramidModel = u_load('C:\isbe\dev\background\results\pyramid\normal_2_model_5_1');
hist_args.Border = 20;
[sum_elements] = mb_find_cluster_support(hist_args);
clust_hist_5_1(end+1) = sum_elements - sum(clust_hist_5_1);
save C:\isbe\dev\background\misc\clust_hist_5_1 clust_hist_5_1
%
hist_args.Level = 5;
hist_args.Orientation = 2;
hist_args.PyramidModel = u_load('C:\isbe\dev\background\results\pyramid\normal_2_model_5_2');
hist_args.Border = 20;
[sum_elements] = mb_find_cluster_support(hist_args);
clust_hist_5_2(end+1) = sum_elements - sum(clust_hist_5_2);
save C:\isbe\dev\background\misc\clust_hist_5_2 clust_hist_5_2
%
hist_args.Level = 5;
hist_args.Orientation = 3;
hist_args.PyramidModel = u_load('C:\isbe\dev\background\results\pyramid\normal_2_model_5_3');
hist_args.Border = 20;
[sum_elements] = mb_find_cluster_support(hist_args);
clust_hist_5_3(end+1) = sum_elements - sum(clust_hist_5_3);
save C:\isbe\dev\background\misc\clust_hist_5_3 clust_hist_5_3
%
hist_args.Level = 5;
hist_args.Orientation = 4;
hist_args.PyramidModel = u_load('C:\isbe\dev\background\results\pyramid\normal_2_model_5_4');
hist_args.Border = 20;
[sum_elements] = mb_find_cluster_support(hist_args);
clust_hist_5_4(end+1) = sum_elements - sum(clust_hist_5_4);
save C:\isbe\dev\background\misc\clust_hist_5_4 clust_hist_5_4
%%
%
% Construct pyramids
[pyramid5 p_sizes5] = buildSFpyr(double(im5), 5, 4);
[pyramid10 p_sizes10] = buildSFpyr(double(im10), 5, 4);

% Change pyramid form
pyramid5 = mb_change_pyramid_form(pyramid5, p_sizes5);
pyramid10 = mb_change_pyramid_form(pyramid10, p_sizes10);

%[r5 c5] = find(circ5);
[r10 c10] = find(tight_bw);
r5 = r10; c5 = c10;

idx5_2 = sub2ind(size(pyramid5{2,1}), r5, c5);
idx10_2 = sub2ind(size(pyramid10{2,1}), r10, c10);

r5_3 = ceil(r5/2); c5_3 = ceil(c5/2);
r10_3 = ceil(r10/2); c10_3 = ceil(c10/2);

idx5_3 = unique(sub2ind(p_sizes5(7,:), r5_3, c5_3)); %level 3 indices
idx10_3 = unique(sub2ind(p_sizes10(7,:), r10_3, c10_3));

r5_4 = ceil(r5/4); c5_4 = ceil(c5/4);
r10_4 = ceil(r10/4); c10_4 = ceil(c10/4);

idx5_4 = unique(sub2ind(p_sizes5(12,:), r5_4, c5_4)); %level 4 indices
idx10_4 = unique(sub2ind(p_sizes10(12,:), r10_4, c10_4));

%
%swap pyramid values at given indices
pyramid5a = pyramid5;
pyramid10a = pyramid10;

%change lev 4
for ori = 1:5
    pyramid5a{4,ori}(idx5_4) = pyramid10{4,ori}(idx10_4);
    pyramid10a{4,ori}(idx10_4) = pyramid5{4,ori}(idx5_4);
end

%change lev 3
for ori = 1:5
    pyramid5a{3,ori}(idx5_3) = pyramid10{3,ori}(idx10_3);
    pyramid10a{3,ori}(idx10_3) = pyramid5{3,ori}(idx5_3);
end

%change lev 3
for ori = 1:5
    pyramid5a{2,ori}(idx5_2) = pyramid10{2,ori}(idx10_2);
    pyramid10a{2,ori}(idx10_2) = pyramid5{2,ori}(idx5_2);
end

%
%reconstruct images
pyramid5a = mb_change_pyramid_form(pyramid5a);
pyramid10a = mb_change_pyramid_form(pyramid10a);

im5_r = reconSFpyr(pyramid5a, p_sizes5);
im10_r = reconSFpyr(pyramid10a, p_sizes10);
%
figure; image(im5); colormap(gray(256)); axis image;
figure; image(im5_r); colormap(gray(256)); axis image;
figure; image(im10); colormap(gray(256)); axis image;
figure; image(im10_r); colormap(gray(256)); axis image;
%%
orig_im = imread('C:/isbe/dev/background/images/pectoral_2/pectoral001.bmp');
orig_pyr = u_load('C:/isbe/dev/background/pyramid/pectoral_2/pectoral001_pyramid.mat');
% make mask for circle
%
r = 25;
[rows cols] = size(orig_im);
c_x = round(cols/2); c_y = round(rows/2);
[x y] = meshgrid(1:cols, 1:rows);
circ = (x - c_x).^2 + (y - c_y).^2 > r.^2;
clear r c_x c_y x y
%%
%set up synthesis args
syn_args.FilledImage = circ;   
%make sure uint8 mode is turned off since we're synthesising doubles
syn_args.Uint8Mode = 0;
%%
ori = 5;

%Set arguments specific to current band
syn_args.PathToTextureGMM = ...
    ['C:\isbe\dev\background\results\mc\mc_2_model_2_',...
    num2str(ori), '.mat'];

syn_args.SeededImage = orig_pyr{2, ori};

%Synthesise the new texture for the band and save back to pyramid
[new_pyr{2, ori} cluster_image{2,ori}] = mb_gmm_tex_synthesis(syn_args);

%display new texture
figure; imagesc(new_pyr{2, ori}); axis image; colormap(jet(256));
figure; imagesc(cluster_image{2, ori}); axis image; colormap(jet(256));

%%
clear syn_args;
syn_args.SamplePyramid = orig_pyr;
syn_args.FilledImage = circ;
syn_args.ModelDir = 'C:\isbe\dev\background\results\mc\';
syn_args.ModelName = 'mc_2_model';
syn_args.CutOffLevel = 4;
syn_args.Plot = 0;

syn_args.SaveFile = [mberksroot, 'background\syn\pectoral_mc_model'];
[synthesised_image, new_pyr] = mb_gmm_pyr_synthesis(syn_args);