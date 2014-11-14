function band_cluster_dual_tree_oris(id_idx, data_type, size, w1, w2)
%final cluster wrapper script

%id_idx = 1,...,20 and comes from $SGE_IDX batch job variable
[levels orientations reals] = meshgrid(1:4, [1 4], 0:1);

k = 10;
band_cluster_args.ClusteringFunctionArgs.k = k;
band_cluster_args.ClusteringFunctionArgs.MaxIter = 10000;
band_cluster_args.MaxFinalMemory = size;
band_cluster_args.DualTreeDir = [mberksroot, 'background/dual_tree/', data_type, '/'];
band_cluster_args.WindowSize = w1;
band_cluster_args.WindowSize2 = w2;
band_cluster_args.Standardise = 0;
    
band_cluster_args.FinalFile = [mberksroot, 'background/models/dual_tree/',...
    data_type, '_k', num2str(k),...
    '_w1_', num2str(w1),...
    '_w2_', num2str(w2),...
    '_', num2str(levels(id_idx)),...
    '_', num2str(orientations(id_idx)),...
    '_', num2str(reals(id_idx))];

band_cluster_args.Level = levels(id_idx);
band_cluster_args.Orientation = orientations(id_idx);
band_cluster_args.Real = reals(id_idx);
band_cluster_args.Overlap = 64 / 2^levels(id_idx);

display(['Results file is ', band_cluster_args.FinalFile]);
mb_cluster_dual_tree_oris(band_cluster_args);
clear;

if ~ispc, exit; end


