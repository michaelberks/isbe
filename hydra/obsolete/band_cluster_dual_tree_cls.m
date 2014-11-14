function band_cluster_dual_tree_cls(id_idx, data_type, size, w1, w2, k, r)
%final cluster wrapper script
if nargin < 7
    r = 100;
end

%id_idx = 1,...,4 and comes from $SGE_IDX batch job variable
[levels clss] = meshgrid(1:5, 0:1);
level = levels(id_idx);
cls = clss(id_idx);

%Change the random seed to we don't select the same data points every time
rand('twister', sum(100*clock));

band_cluster_args.ClusteringFunctionArgs.Replicates = r;
band_cluster_args.ClusteringFunctionArgs.k = k;
band_cluster_args.ClusteringFunctionArgs.MaxIter = 10000;
band_cluster_args.MaxFinalMemory = size;
band_cluster_args.DualTreeDir = [mberksroot, 'background/dual_tree/', data_type, '/'];
band_cluster_args.CLSDir = [mberksroot, 'background/cls/', data_type, '/'];
band_cluster_args.CLSMode = cls;
band_cluster_args.WindowSize = w1;
band_cluster_args.WindowSize2 = w2;
band_cluster_args.Standardise = 0;

band_cluster_args.FinalFile = [mberksroot, 'background/models/dual_tree/',...
    data_type, '_k', num2str(k),...
    '_w1_', num2str(w1),...
    '_w2_', num2str(w2),...
    '_', num2str(level),...
    '_cls_', num2str(cls)];

band_cluster_args.Level = level;
band_cluster_args.Overlap = 64 / 2^level;

display(['Results file is ', band_cluster_args.FinalFile]);
mb_cluster_dual_tree_cls(band_cluster_args);
clear;

if ~ispc, exit; end


