function band_cluster_pyramid_cls(id_idx, data_type)
%final cluster wrapper script

%id_idx = 1,...,30 and comes from $SGE_IDX batch job variable
[levels orientations cls_mode] = meshgrid(2:4, 1:5, 0:1);

k = 10;
band_cluster_args.ClusteringFunctionArgs.k = k;
band_cluster_args.ClusteringFunctionArgs.MaxIter = 10000;
band_cluster_args.MaxFinalMemory = 512;
band_cluster_args.PyramidsDir = [mberksroot, 'background/pyramid/', data_type, '/'];
band_cluster_args.CLSDir = [mberksroot, 'background/cls/', data_type, '/'];
band_cluster_args.WindowSize = 11;
band_cluster_args.WindowSize2 = 0;
band_cluster_args.Standardise = 0;
  
band_cluster_args.FinalFile = [mberksroot, 'background/models/cls/',...
    data_type, '_cls_k', num2str(k), '_c_model_',...
    num2str(levels(id_idx)) '_',...
    num2str(orientations(id_idx)) '_',...
    num2str(cls_mode(id_idx))];

band_cluster_args.Level = levels(id_idx);
band_cluster_args.Orientation = orientations(id_idx);
band_cluster_args.CLSMode = cls_mode(id_idx);

band_cluster_args.Overlap = 64 / 2^(levels(id_idx) - 1);

display(band_cluster_args);
display(['Results file is ', band_cluster_args.FinalFile]);
mb_cluster_cls(band_cluster_args);
clear;
exit;


