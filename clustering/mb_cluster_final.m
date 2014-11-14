function mb_cluster_final(varargin)

% Clustering method based on the divide and conquer method implemented by
% CJR. This function implements the final gathering of the divided cluster
% results and performs the final clustering

% get a dir listing in the results dir of files matching the naming
% convention of our cluster results

args = u_packargs(varargin,... % the user's input
			 '0', ... % strict mode
			 {'ClusteringFunctionArgs',...
             'ResultsDir',...
             'FinalFile'}, ... % the mandatory arguments
			 'ClusteringFunction', 'mb_k_means_clustering2',...
             'ResultFile', 'clustering_result',...
             'MaxFinalMemory', 512,...
             'Dimensions', 15^2);
clear varargin;

if args.ResultsDir(end) ~= filesep
	args.ResultsDir = [args.ResultsDir filesep];
end
% get a dir listing in the results dir of files matching the naming
% convention of our cluster results        
results_list = dir([args.ResultsDir, args.ResultFile, '*.mat']);
num_results = length(results_list);
display(['Number of results to group = ', num2str(num_results)]);

% randomly select one of the cluster results to supply the initial cluster
% means
% init_mean_idx = randsample(num_results, 1);

% work out how many data points we can fit in memory for the final clustering
% each element occupies 2^-17 MB (i.e. 8 bytes), each point has args.Dimensions elements 
total_points = floor(args.MaxFinalMemory*2^17 / args.Dimensions);

% work out how many data points to take from each file -- in case we need to sample from each file
num_points_from_each_file = floor(total_points / num_results);

% add each set of rep points to the central pool; pre-allocate memory -
% better to be too big now, and free up spare later
args.ClusteringFunctionArgs.Data = zeros(total_points, args.Dimensions);
curr_point = 1; %maintain counter to current point
for ii = 1 : num_results
	load([args.ResultsDir   results_list(ii).name]); % variable is called clustering_results
	
	% we have to select a subset of the representative points to use in the final clustering
    rp_size = size(clustering_result.RepPoints, 1);
    if num_points_from_each_file < rp_size %#ok
        random_indices = randsample(rp_size, num_points_from_each_file);
        args.ClusteringFunctionArgs.Data(curr_point:curr_point + num_points_from_each_file - 1,:) =...
            clustering_result.RepPoints(random_indices, :);
        curr_point = curr_point + num_points_from_each_file;
    else
        args.ClusteringFunctionArgs.Data(curr_point:curr_point + rp_size - 1,:) =...
            clustering_result.RepPoints;
        curr_point = curr_point + rp_size;
	end
	
    % set initial cluster means
%     if ii == init_mean_idx;
%         args.ClusteringFunctionArgs.Start = clustering_result.Means;
%     end
    clear rep_points clustering_results;
end
% discard unused rows of data
args.ClusteringFunctionArgs.Data(curr_point:end, :) = [];

disp(datestr(now))
disp(['Large Data Clustering: About to perform the final clustering on ' num2str(size(args.ClusteringFunctionArgs.Data,1)) ' data points'])

% Do the final clustering
final_clustering_result = feval(args.ClusteringFunction, args.ClusteringFunctionArgs); %#ok
% Get the data used in the final clustering, in case we got all this way
% and produced a bad result
final_clustering_result.RepresentativeDataPoints = args.ClusteringFunctionArgs.Data;

% Save the clustering_result
save(args.FinalFile, 'final_clustering_result');




