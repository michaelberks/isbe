function band_cluster_pyramid_lower(id_idx, data_type)
%final cluster wrapper script

%id_idx = 1,...,5 and comes from $SGE_IDX batch job variable

band_cluster_args.ClusteringFunctionArgs.k = 10;
band_cluster_args.ClusteringFunctionArgs.MaxIter = 10000;
band_cluster_args.MaxFinalMemory = 256;
band_cluster_args.PyramidsDir = [mberksroot, 'background/pyramid/', data_type, '/'];
band_cluster_args.WindowSize = 5;
band_cluster_args.WindowSize2 = 3;

   
band_cluster_args.FinalFile = [mberksroot, 'background/models/',...
    data_type, '_20_model_6_', num2str(id_idx)];
band_cluster_args.Level = 6;
band_cluster_args.Orientation = id_idx;
band_cluster_args.Overlap = 4;

display(['Results file is ', band_cluster_args.FinalFile]);
mb_cluster(band_cluster_args);
clear;
exit;


