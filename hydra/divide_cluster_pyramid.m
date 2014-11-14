function divide_cluster_pyramid(idx_id, data_name)
%wrapper function to call the divide stage of clustering algorithm

divide_cluster_args.ResultFile = ...
    [mberksroot, 'background/results/pyramid/',data_name,'/clustering_result', zerostr(idx_id, 3)];

divide_cluster_args.NextDataFunction = 'get_pyramid_data_for_clustering';
divide_cluster_args.NextDataFunctionArgs.IdxFile = ...
    [mberksroot, 'background/idx/',data_name,'/cluster_idx', zerostr(idx_id, 3)];
divide_cluster_args.NextDataFunctionArgs.ImageDir = ...
    [mberksroot, '/background/pyramid/',data_name,'/'];

divide_cluster_args.ClusteringFunctionArgs.k = 40;
divide_cluster_args.ClusteringFunctionArgs.MaxIter = 1000;
divide_cluster_args.ClusteringFunctionArgs.PropRepresentativePoints = 0.4;

display(['Idx file is', divide_cluster_args.NextDataFunctionArgs.IdxFile]);
display(['Result file is', divide_cluster_args.ResultFile]);
display(['Image dir is', divide_cluster_args.NextDataFunctionArgs.ImageDir]);

mb_cluster_once(divide_cluster_args);

clear;
exit;