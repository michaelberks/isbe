%Script for mucking around with CLS synthesis
pyr_args.TargetPyramid = u_load([mberksroot,...
    'background/g_pyramid/normal_2/normal007_pyramid.mat']);
[rows cols] = size(pyr_args.TargetPyramid{1,1});
row_centre = round(rows / 2);
col_centre = round(cols / 2);
%     pyr_args.FilledImage = logical(ones(rows, cols)); %#ok
%     pyr_args.FilledImage(row_centre-64:row_centre+63, col_centre-64:col_centre+63) = 0;

% Make a the biggest circular mask leaving an overlap of 48 at all edges
m = min([rows cols]);
rad = floor((m - 96) / 2);
[x y] = meshgrid(1:cols, 1:rows);
pyr_args.FilledImage = (x - col_centre).^2 + (y - row_centre).^2 > rad.^2;
clear x y m rad row_centre col_centre;

pyr_args.ModelName = 'mass_2_g_20_model';
pyr_args.ModelNameCls = 'mass_2_g_cls_20_model';

% cls_result = u_load(['C:/isbe/dev/background/cls/mass_2/mass_cls_', zerostr(ii, 3), '.mat']);

cls_007_3 = double(rgb2gray(imread(('C:/isbe/dev/background/misc/cls_temp007_3.bmp'))));
cls_result{1,1}.CLS = zeros(size(pyr_args.TargetPyramid{1,1}));
cls_result{2, 1}.CLS = zeros(size(pyr_args.TargetPyramid{2,1}));
cls_result{3, 1}.CLS = ~cls_007_3;

pyr_args.ClsMap = mb_cls_map_merge(cls_result);

%only synthesis CLS pixels
pyr_args.FilledImage = pyr_args.FilledImage | ~pyr_args.ClsMap{1};
figure; imagesc(pyr_args.FilledImage); axis image; colormap(gray(256));
%%

pyr_args.ModelDir = [mberksroot, 'background/results/g_pyramid/'];
pyr_args.CutOffLevel = 4;
pyr_args.CutOffLevelCls = 3;
pyr_args.ConditionLevels = 1;

[synthesised_image, pyramid, cluster_image] = mb_gmm_g_pyr_synthesis_cls(pyr_args);

%%
model_3_1 = u_load('C:/isbe/dev/background/results/g_pyramid/mass_2_g_cls_20_model_3_1.mat');

for ii = 1:20
    mean_window = reshape(model_3_1.Means(ii, 1:121), 11, 11);
    mean_window = (mean_window * model_3_1.Stats.SD) + model_3_1.Stats.Mean;
    figure; imagesc(mean_window); axis image; colormap(gray(256)); colorbar
end
%%
cls_results_syn = cell(5,5);
cls_results_syn{1,1} = CLS_result_syn_2_1;
cls_results_syn{1,2} = CLS_result_syn_2_2;
cls_results_syn{1,3} = CLS_result_syn_2_3;
cls_results_syn{1,4} = CLS_result_syn_2_4;
cls_results_syn{1,5} = CLS_result_syn_2_5;
%%
for level = 3:4
    for ori = 1:5
        cls_results_ori{level-1, ori} = ...
            mb_cls_selection('ImageIn', pyramid{level,ori}, 'GaborFilterArgs', GaborFilterArgs);
        cls_results_syn{level-1, ori} = ...
            mb_cls_selection('ImageIn', syn.pyramid{level,ori}, 'GaborFilterArgs', GaborFilterArgs);
    end
end

%%
for level = 1:3
    for ori = 1:5
        figure;
        subplot(1,2,1); imagesc(cls_results_ori{level, ori}.CLS); axis image; colormap(gray(256));
        subplot(1,2,2); imagesc(cls_results_syn{level, ori}.CLS); axis image; colormap(gray(256));
    end
end
        

