function band_cluster_g_pyramid(id_idx, data_type)
%final cluster wrapper script

%id_idx = 1,...,4 and comes from $SGE_IDX batch job variable

band_cluster_args.ClusteringFunctionArgs.k = 40;
band_cluster_args.ClusteringFunctionArgs.MaxIter = 10000;
band_cluster_args.MaxFinalMemory = 512;
band_cluster_args.PyramidsDir = [mberksroot, 'background/g_pyramid/', data_type, '/'];
band_cluster_args.WindowSize = 11;
band_cluster_args.WindowSize2 = 5;

   
band_cluster_args.FinalFile = [mberksroot, 'background/results/g_pyramid/',...
    data_type, '_g_20_model_', num2str(id_idx)];
band_cluster_args.Level = id_idx;
band_cluster_args.Overlap = 64 / 2^(id_idx - 1);

display(['Results file is ', band_cluster_args.FinalFile]);
mb_cluster_gp(band_cluster_args);
clear;
exit;


