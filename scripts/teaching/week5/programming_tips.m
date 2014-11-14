%--------------------------------------------------------------------------
% ------------------------- Displaying Images -----------------------------
%--------------------------------------------------------------------------

%Script showing examples of good code practice
data_dir = 'P:\MatlabTutorial\data\';
load([data_dir '024RCC.mat']); %loads in a data structure mammogram
load([data_dir '024RCC_cancer_map.mat'])
%%
figure; imagesc(mammogram); axis image; colormap(gray(256));
figure; imagesc(cancer_map); axis image; colormap(gray(256));

%%
%Now try and use the function: localimagemaximanocomments to find all the
%points in the cancer map with a score greater than 25, larger than all
%their neighbours (i.e. the local maxima)
localimagemaximanocomments

%%
%Now try and use the function local_image_maxima_comments
[maxima_pos, maxima_vals] = local_image_maxima_comments(cancer_map, 50, [], 25);

plot(maxima_pos(:,1), maxima_pos(:,2), 'kx', 'markersize', 10);
plot(maxima_pos(1,1), maxima_pos(1,2), 'ko', 'markersize', 10);

%%
mammo_small = imresize(mammogram, 0.25, 'bilinear');
figure; imagesc(mammo_small); axis image; colormap(gray(256));
%
tic;
[filtered_mammogram1] = median_filter_circle1(mammo_small, 5);
toc;
figure; imagesc(filtered_mammogram1); axis image; colormap(gray(256));
%
tic;
[filtered_mammogram2] = median_filter_circle2(mammo_small, 5);
toc;
figure; imagesc(filtered_mammogram1); axis image; colormap(gray(256));
%
tic;
[filtered_mammogram3] = median_filter_circle3(mammo_small, 5);
toc;
figure; imagesc(filtered_mammogram1); axis image; colormap(gray(256));
%%
profile clear;
profile on;
[filtered_mammogram1] = median_filter_circle1(mammo_small, 5);
profile report;
%
profile clear;
profile on;
[filtered_mammogram2] = median_filter_circle2(mammo_small, 5);
profile report;
%
profile clear;
profile on;
[filtered_mammogram3] = median_filter_circle3(mammo_small, 5);
profile report;
