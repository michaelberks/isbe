function depth_cluster_dual_tree(data_type, size, k, id_idx)
%final cluster wrapper script

%set rand seed to unique id_idx
rand('twister', id_idx);

band_cluster_args.ClusteringFunctionArgs.k = k;
band_cluster_args.ClusteringFunctionArgs.MaxIter = 10000;
band_cluster_args.MaxFinalMemory = size;
band_cluster_args.DualTreeDir = [mberksroot, 'background/dual_tree/', data_type, '/'];
    
band_cluster_args.FinalFile = [mberksroot, 'background/models/dual_tree/',...
    data_type, '_k', num2str(k),...
    '_size_', num2str(size),...
    '_uid_', num2str(id_idx), '_std'];

display(['Results file is ', band_cluster_args.FinalFile]);
mb_cluster_dual_tree(band_cluster_args);
clear;

if ~ispc, exit; end


