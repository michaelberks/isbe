
cluster_args.ClusteringFunctionArgs.k = 20;
cluster_args.ClusteringFunctionArgs.MaxIter = 1000;

cluster_args.NextDataFunctionArgs.LoadIdx = 1;
cluster_args.NextDataFunctionArgs.Idx = 'C:\isbe\dev\background\misc\mc_idx1';
cluster_args.NextDataFunctionArgs.ImageDir = 'C:\isbe\dev\background\pyramid\mc_2\';
cluster_args.NextDataFunctionArgs.Level = 2;
cluster_args.NextDataFunction = 'get_band_data_for_clustering';

for ii = 1:5
    display(['Clustering orientation ', num2str(ii), ' at ', datestr(now)]);
    cluster_args.NextDataFunctionArgs.Orientation = ii; 
    cluster_args.ResultFile = ['C:\isbe\dev\background\results\mc_pyramid_',...
        num2str(cluster_args.NextDataFunctionArgs.Level), '_',...
        num2str(cluster_args.NextDataFunctionArgs.Orientation)];
    mb_cluster_once(cluster_args);
end
%%
ch = clust_hist_5_4;
model = u_load('C:/isbe/dev/background/results/pyramid/normal_2_model_5_4.mat');

ch(1) = [];
pch = 100*ch / sum(ch);
scrap_idx = find(pch < 0.1);
display(['num cluster to scrap = ', num2str(length(scrap_idx))]);

model.Means(scrap_idx, :) = [];
model.CovMats(scrap_idx) = [];
model.NumClusters = model.NumClusters - length(scrap_idx);
model.ClusterProbs(scrap_idx) = [];
save C:/isbe/dev/background/results/pyramid/normal_2_model_a_5_4.mat model

%%
syn_args.SamplePyramid = u_load('C:/isbe/dev/background/pyramid/normal_2/normal006.bmp_pyramid.mat');
%syn_args.SampleImage = orig_image;
syn_args.FilledImage = ones(size(syn_args.SamplePyramid{1,1}));
syn_args.FilledImage(100:249, 100:249) = 0;
syn_args.ModelDir = 'C:/isbe/dev/background/results/pyramid/';
syn_args.ModelName = 'normal_2_model_a';
syn_args.CutOffLevel = 4;
syn_args.Plot = 0;

syn_args.SaveFile = [mberksroot, 'background\syn\syn_image_14of15'];
[synthesised_image, pyramid] = mb_gmm_pyr_synthesis(syn_args);
%%
NextDataFunctionArgs.LoadIdx = 1;
NextDataFunctionArgs.Idx = 'C:\isbe\dev\background\misc\mc_idx1';
NextDataFunctionArgs.ImageDir = 'C:\isbe\dev\background\pyramid\mc_2\';
NextDataFunctionArgs.Level = 2;
NextDataFunctionArgs.Orientation = 1;
data = get_band_data_for_clustering(NextDataFunctionArgs);
%%
args.LoadIdx = 1;
args.Idx = 'C:\isbe\dev\background\idx\normal_2\cluster_idx001';
args.ImageDir = 'C:\isbe\dev\background\pyramid\normal_2\';
args.Level = 3;
args.MaxMemory = 1;
data = get_complete_data_for_clustering(args);
%%
args.LoadIdx = 1;
args.Idx = 'C:\isbe\dev\background\idx\normal_2\cluster_idx001';
args.ImageDir = 'C:\isbe\dev\background\pyramid\normal_2\';
args.CutOffLevel = 4;
args.MaxMemory = 1;
[model] = mb_build_gaussian_pyramid_model(args);
%%
syn_args.PathToPyramidGM = 'C:\isbe\dev\background\misc\test_model';
syn_args.SampleImage = imread('C:\isbe\dev\background\images\normal_2\normal005.bmp');
syn_args.FilledImage = ones(size(syn_args.SampleImage));
syn_args.FilledImage(100:149, 100:149) = 0;

profile on;
[synthesised_image, pyramid] = mb_gmm_pyr_synthesis3(syn_args);
profile viewer;
%%
for ii = 1:10
    p_sample = (randn(length(L_x), 1) .* sqrt(L_x'));
    new_win = mean_x + (P_x*p_sample)';
    figure; imagesc(reshape(new_win, 15, 15)); axis image; colormap(jet(256));
end
%%
for ii = 2:4
    for jj = 1:5
        try
            load(['C:\isbe\dev\background\results\pyramid\normal_2_model_a_',...
                num2str(ii), '_', num2str(jj), '.mat']);
        catch
            continue
        end
        [cluster_image] = assign_cluster_to_image(model, pyramid{ii, jj}); 
        figure; imagesc(cluster_image); colormap(jet(256)); axis image; caxis([0 40]);
    end
end