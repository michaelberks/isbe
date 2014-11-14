function band_cluster_g_pyramid_cls(id_idx, data_type)
%final cluster wrapper script

%id_idx = 1,...,8 and comes from $SGE_IDX batch job variable

[levels cls_mode] = meshgrid(1:4, 0:1);

band_cluster_args.ClusteringFunctionArgs.k = 40;
band_cluster_args.ClusteringFunctionArgs.MaxIter = 10000;
band_cluster_args.MaxFinalMemory = 512;
band_cluster_args.PyramidsDir = [mberksroot, 'background/g_pyramid/', data_type, '/'];
band_cluster_args.CLSDir = [mberksroot, 'background/cls/', data_type, '/'];
band_cluster_args.WindowSize = 11;
band_cluster_args.WindowSize2 = 5;
band_cluster_args.CLSMode = cls_mode(id_idx);
   
band_cluster_args.FinalFile = [mberksroot, 'background/results/g_pyramid/',...
    data_type, '_g_cls_20_model_', num2str(levels(id_idx)) '_', num2str(cls_mode(id_idx))];
band_cluster_args.Level = levels(id_idx);
band_cluster_args.Overlap = 64 / 2^(levels(id_idx) - 1);

display(band_cluster_args);
display(['Results file is ', band_cluster_args.FinalFile]);
mb_cluster_gp_cls(band_cluster_args);
clear;
exit;


