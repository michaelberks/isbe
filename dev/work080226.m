clear

pyramid1 = u_load('C:\isbe\dev\background\pyramid\normal_2\normal001.bmp_pyramid.mat');
filled_image = logical(ones(size(pyramid1{2,1})));

row_centre = round(size(pyramid1{2,1}, 1) / 2);
col_centre = round(size(pyramid1{2,1}, 2) / 2);
filled_image(row_centre-64:row_centre+64-1, col_centre-64:col_centre+64-1) = 0;

wei_args.FilledImage = filled_image;

wei_args.SampleData = u_load('C:\isbe\dev\background\trees\normal_2\normal004.bmp_pyramid_tree\tsvq_names.mat');
wei_args.TargetPyramid = pyramid1;
wei_args.WindowSize1 = 5;
wei_args.WindowSize2 = 9;
wei_args.Plot = true;
wei_args.CutOffLevel = 3;
wei_args.K_method = 'tree';

for ii = 2:5
    for jj = 1:5
        for kk = 1:2
            wei_args.SampleData{ii,jj,kk} =...
                [mberksroot wei_args.SampleData{ii,jj,kk}(31:end)];
        end
    end
end
%profile on
[syn_im_tree, new_pyr_tree] =  mb_wei_pyr_synthesis(wei_args);
%profile off

%%
gp_1 = u_load('C:\isbe\dev\background\g_pyramid\normal001_pyramid.mat');
%%
for ii = 1:5
    figure; image(gp_1{ii}); colormap(gray(256)); axis image;
end

%%
clear

pyramid1 = u_load('C:\isbe\dev\background\g_pyramid\lines\line001_pyramid.mat');
filled_image = logical(ones(size(pyramid1{1})));

row_centre = round(size(pyramid1{1}, 1) / 2);
col_centre = round(size(pyramid1{1}, 2) / 2);
filled_image(row_centre-64:row_centre+64-1, col_centre-64:col_centre+64-1) = 0;

wei_args.FilledImage = filled_image;
%
wei_args.SampleData = u_load('C:\isbe\dev\background\g_pyramid\lines\tsvq\line003_pyramid_tree\tsvq_names.mat');
wei_args.TargetPyramid = pyramid1;
wei_args.WindowSize1 = 5;
wei_args.WindowSize2 = 9;
wei_args.Plot = true;
wei_args.CutOffLevel = 3;
wei_args.K_method = 'tree';

[syn_im_tree, new_pyr_tree] =  mb_wei_g_pyr_synthesis(wei_args);
%%
pyr_args.TargetPyramid = u_load('C:\isbe\dev\background\pyramid\normal_2\normal001.bmp_pyramid.mat');
pyr_args.FilledImage = logical(ones(size(pyr_args.TargetPyramid{2,1})));
row_centre = round(size(pyr_args.TargetPyramid{1}, 1) / 2);
col_centre = round(size(pyr_args.TargetPyramid{1}, 2) / 2);
pyr_args.FilledImage(row_centre-64:row_centre+64-1, col_centre-64:col_centre+64-1) = 0;

pyr_args.ModelDir = 'C:\isbe\dev\background\results\2levels\';
pyr_args.ModelName = 'normal_2_model';
%%
profile on;
[synthesised_image, pyramid, cluster_image] = mb_gmm_pyr_synthesis(pyr_args);
profile viewer;