function band_cluster_dual_tree_all(id_idx, data_type, size, w1, w2, r)
%final cluster wrapper script

%id_idx = 1,...,4 and comes from $SGE_IDX batch job variable
[levels ks] = meshgrid(1:5, 10:10:50);
level = levels(id_idx);
k = ks(id_idx);



%Change the random seed to we don't select the same data points every time
rand('twister', sum(100*clock));

band_cluster_args.ClusteringFunctionArgs.k = k;
band_cluster_args.ClusteringFunctionArgs.Replicates = r;
band_cluster_args.ClusteringFunctionArgs.MaxIter = 10000;
band_cluster_args.MaxFinalMemory = size;
band_cluster_args.DualTreeDir = [mberksroot, 'background/dual_tree/', data_type, '/'];
band_cluster_args.WindowSize = w1;
band_cluster_args.WindowSize2 = w2;
band_cluster_args.Standardise = 0;

band_cluster_args.FinalFile = [mberksroot, 'background/models/dual_tree/',...
    data_type, '_k', num2str(k),...
    '_rep_', num2str(r),...
    '_w1_', num2str(w1),...
    '_w2_', num2str(w2),...
    '_', num2str(level)];

band_cluster_args.Level = level;
band_cluster_args.Overlap = 64 / 2^level;

display(['Results file is ', band_cluster_args.FinalFile]);
mb_cluster_dual_tree_final(band_cluster_args);
clear;

if ~ispc, exit; end


