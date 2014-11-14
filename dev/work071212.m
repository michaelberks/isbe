m_list = dir('C:\isbe\dev\background\images\mc\mc*');

for ii = 1:10
    load(['C:\isbe\dev\background\images\mc\', m_list(ii).name]);
    figure; imagesc(mass.mass_ROI); colormap(gray(256)); axis image;
    imwrite(mass.mass_ROI, ['C:\isbe\dev\background\images\mc\mc',...
        zerostr(ii,2), '.bmp']);
    clear mass;
end
%%
%
for ii = 1:10
    im1 = imread(['C:\isbe\dev\background\images\mc\mc',...
        zerostr(ii,2), '.bmp']);
    p_small = imresize(im1, 0.5, 'bicubic');
    imwrite(p_small, ['C:\isbe\dev\background\images_small\mc\mc',...
        zerostr(ii,2), '.bmp']);
    clear mass;
end
%
total_pts = 0;
for ii = 1:10
    im1 = imread(['C:\isbe\dev\background\images_small\mc\mc',...
        zerostr(ii,2), '.bmp']);
    total_pts = numel(im1) + total_pts;
    clear mass;
end
%%
idx_file =...
    ['C:\isbe\dev\background\idx\pectoral\cluster_idx', zerostr(1, 3)];
result_file =...
    ['C:\isbe\dev\background\results\pectoral\cluster_result', zerostr(1, 3)];

display(['Idx file is', idx_file]);
display(['Result file is', result_file]);

load('C:\isbe\dev\background\args\divide_cluster_args');
divide_cluster_args.IdxFile = idx_file;
divide_cluster_args.ResultFile = result_file;
divide_cluster_args.ImageDir = 'C:\isbe\dev\background\images_small\pectoral\';

mb_cluster_once(divide_cluster_args);
%%
clear;

final_cluster_args.ClusteringFunctionArgs.k = 40;
final_cluster_args.ClusteringFunctionArgs.MaxIter = 2;
final_cluster_args.FinalFile = [mberksroot, 'background/results/pectoral_model'];
final_cluster_args.ResultsDir = [mberksroot, 'background/results/pectoral/'];
final_cluster_args.MaxFinalMemory = 256;

%%
save C:\isbe\dev\background\args\final_cluster_args.mat final_cluster_args
%%
divide_cluster_args.IdxFile = [mberksroot, 'background/idx/pectoral/cluster_idx001'];
divide_cluster_args.ResultFile = [mberksroot, 'background/results/del_me'];
divide_cluster_args.ClusteringFunctionArgs.k = 40;
divide_cluster_args.ClusteringFunctionArgs.MaxIter = 2;
divide_cluster_args.ClusteringFunctionArgs.PropRepresentativePoints = 0.2;
divide_cluster_args.ImageDir = [mberksroot, '/background/images_small/pectoral/'];
divide_cluster_args.WindowSize = 15;
%%
save C:\isbe\dev\background\args\divide_cluster_args.mat divide_cluster_args
%%
cluster_idx_args.NumberSets = 12;
cluster_idx_args.IdxDir = 'C:\isbe\dev\background\idx\mc\';
cluster_idx_args.NextIdxFunctionArgs.ImageDir = 'C:\isbe\dev\background\images_small\mc\';
cluster_idx_args.NextIdxFunctionArgs.TempStorageDir = 'C:\isbe\dev\background\temp\';
cluster_idx_args.NextIdxFunctionArgs.WindowSize = 15;
cluster_idx_args.NextIdxFunctionArgs.MaxMemory = 256;
cluster_idx_args.NextIdxFunction = 'mb_next_idx_for_clustering';
