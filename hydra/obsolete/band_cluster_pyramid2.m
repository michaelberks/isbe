function band_cluster_pyramid2(data_type)
%final cluster wrapper script

band_cluster_args.ClusteringFunctionArgs.k = 40;
band_cluster_args.ClusteringFunctionArgs.MaxIter = 1000;
band_cluster_args.FinalFile = [mberksroot, 'background/results/',...
    data_type, '_model_1_2'];
band_cluster_args.MaxFinalMemory = 512;
band_cluster_args.PyramidsDir = [mberksroot, 'background/pyramid/', data_type, '/'];
band_cluster_args.Level = 1;
band_cluster_args.Orientation = 1;

display(['Results file is ', band_cluster_args.FinalFile]);

mb_cluster(band_cluster_args);
clear;
exit;


