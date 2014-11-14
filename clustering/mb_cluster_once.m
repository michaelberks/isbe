function mb_cluster_once(varargin)

% Clustering method based on the divide and conquer method implemented by
% CJR. This function implements inital clustering on one block. This is so 

args = u_packargs(varargin,... % the user's input
			 '0', ... % strict mode
			 {'ResultFile',...
             'ClusteringFunctionArgs',...
             'NextDataFunctionArgs'}, ... % the mandatory arguments
			 'ClusteringFunction', 'mb_k_means_clustering2',...
             'NextDataFunction', 'get_window_data_for_clustering');
clear varargin;

%get data using specified function
args.ClusteringFunctionArgs.Data =...
    feval(args.NextDataFunction, args.NextDataFunctionArgs);

%display how points we're going to cluster
display(['number of data points to cluster = ', num2str(size(args.ClusteringFunctionArgs.Data,1))]);

%cluster the data
clustering_result = feval(args.ClusteringFunction,... 
    args.ClusteringFunctionArgs); %#ok

%save the result
save(args.ResultFile, 'clustering_result');




            
        
