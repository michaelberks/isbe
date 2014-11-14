function band_cluster_pyramid_window_size(id_idx, data_type)
%final cluster wrapper script

%id_idx = 1,...,30 and comes from $SGE_IDX batch job variable
[levels window_sizes] = meshgrid(2:5, [5 7 9 11]);

k = 10;
band_cluster_args.ClusteringFunctionArgs.k = k;
band_cluster_args.ClusteringFunctionArgs.MaxIter = 10000;
band_cluster_args.MaxFinalMemory = 512;
band_cluster_args.PyramidsDir = [mberksroot, 'background/pyramid/', data_type, '/'];
band_cluster_args.Orientation = 1;
band_cluster_args.WindowSize2 = 0;
band_cluster_args.Standardise = 0;
  
band_cluster_args.FinalFile = [mberksroot, 'background/models/window_size/',...
    data_type, '_level_',...
    num2str(levels(id_idx)) '_win_',...
    num2str(window_sizes(id_idx))];

band_cluster_args.WindowSize = window_sizes(id_idx);
band_cluster_args.Level = levels(id_idx);
band_cluster_args.Overlap = 2^(8 - levels(id_idx));

display(band_cluster_args);
display(['Results file is ', band_cluster_args.FinalFile]);
mb_cluster(band_cluster_args);
clear;
exit;


