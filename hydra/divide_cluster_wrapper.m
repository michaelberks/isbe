function divide_cluster_wrapper(idx_id, data_name)
%wrapper function to call the divide stage of clustering algorithm

divide_cluster_args.IdxFile = [mberksroot, 'background/idx/',data_name,'/cluster_idx', zerostr(idx_id, 3)];
divide_cluster_args.ResultFile = [mberksroot, 'background/results/',data_name,'/clustering_result', zerostr(idx_id, 3)];
divide_cluster_args.ClusteringFunctionArgs.k = 40;
divide_cluster_args.ClusteringFunctionArgs.MaxIter = 1000;
divide_cluster_args.ClusteringFunctionArgs.PropRepresentativePoints = 0.2;
divide_cluster_args.ImageDir = [mberksroot, '/background/images_small/',data_name,'/'];
divide_cluster_args.WindowSize = 15;


display(['Idx file is', divide_cluster_args.IdxFile]);
display(['Result file is', divide_cluster_args.ResultFile]);
display(['Image dir is', divide_cluster_args.ImageDir]);

mb_cluster_once(divide_cluster_args);

clear;
exit;