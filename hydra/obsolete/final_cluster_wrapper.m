function final_cluster_wrapper(data_name)
%final cluster wrapper script

final_cluster_args.ClusteringFunctionArgs.k = 40;
final_cluster_args.ClusteringFunctionArgs.MaxIter = 1000;
final_cluster_args.FinalFile = [mberksroot, 'background/results/',data_name,'_model'];
final_cluster_args.ResultsDir = [mberksroot, 'background/results/',data_name,'/'];
final_cluster_args.MaxFinalMemory = 256;
display(['Results dir is ', final_cluster_args.ResultsDir]);

mb_cluster_final(final_cluster_args);
clear;
exit;


