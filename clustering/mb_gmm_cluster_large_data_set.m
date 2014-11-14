function [clustering_result] = mb_gmm_cluster_large_data_set(varargin)
%
% MB_CLUSTER_LARGE_DATA_SET Clusters large data sets.
%
% This function uses a "divide and conquer"  approach to clustering a large
% data set. This function allows the user to provide a clustering function
% that is suitable to their problem, and a function that delivers chunks of
% data, to allow this function to be used in a wide range of applications.
%
% MB_CLUSTER_LARGE_DATA_SET uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% 'k'
%			- the number of clusters to find. Note that the actual number
%			of functions may differ, depending upon the clustering function
%			used (see the 'ClusteringFunction' argument).
% 'NextDataFunction' 
%			- a string which stores the name of a function (M file) that returns the
%			next chunk 	of unseen data to be clustered. So, if the function you
%			wrote to get a new chunk of data is called get_next_chunk.m, then
%			the string should be 'get_next_chunk' This function must return the
%			data in a matrix where each row is a data point and the columns
%			are the dimensions. It must also return a boolean which indicates
%			whether the data sampling has finished. The prototype of the
%			next data function is:
%				[data, finished_sampling_flag] = mynextdatafunc(args_struct)
%					(see the 'NextDataFunctionArgs' argument for details on
%					args_struct)
%			The data returned by each call to the next data function should be
%			as statistically representative of the whole data set as possible; so
%			the function must sample the data and record what has been sampled.
%			The next data function must be able to clean up its temporary files
%			when instructed. This is achieved using the CleanUp member of the
%			NextDataFunctionArgs struct or cell array, which MUST be supported.
%			When CleanUp is set to true (1==1), the function must delete its
%			temporary files and then return (it must not return any data when asked
%			to clean up).
% 'NextDataFunctionArgs'
%			- a struct or cell array containing any arguments required by the
%			'NextDataFunction', i.e. use u_packargs to write the function interface.
% 'ClusteringFunction'
%			- a string which stores the name of a function (M file) that performs
%			clustering on a matrix of data (as returned by the 'NextDataFunction');
%			an example would be a k-means algorithm (see MB_K_MEANS_CLUSTERING).
%			The clustering function should return a cell array where cell i contains a
%			matrix of a 	sample of the data points assigned to cluster i. The number of
%			points 	to be returned should be specified by the 'ClusteringFunctionArgs'
%			struct (see below). The clustering function assumes that the arguments
%			to the clustering function are passed by a struct that has the following
%			fields: 'k', the number of clusters to find and 'Data', a matrix of the data
%			points to clusters (where each row is a data point).These arguments are
%			set by MB_CLUSTER_LARGE_DATA_SET, other arguments must be passed
%			via the 'ClusteringFunctionArgs' argument (see below). This clustering
%			function will be called many times, so it is best if the function is fast.
%			See the 'FinalPassClusteringFunction' argument for details on how to
%			use a more accurate (and probably slower) clustering function to perform
%			a final clustering.
% 'TempDirPath'
%			- a string containing the path to a directory that can be used for temporary
%			storage.
%
%'MaxMemoryFinalClusteringInMB'
%			- the maximum amount of memory to use to store that data for the final
%			clustering, specified in MB.
%
% Optional Arguments:
%
% 'DataTransformationFunction'
%			- a string which stores the name of a function (M file) that performs a
%			transformation on the data returned by the NextDataFunction and before
%			it undergoes the first round of clustering (ClusteringFunction);
%			this will normally be some kind of dimensionality-reduction such as projecting
%			the data into a PCA space, but is can be any useful transformation. See the
%			DataTransformationFunctionArgs documentation (below) to determine how the
%			function must be written. The function must return a matrix of transformed data
%			where each row is a vector and each column is a dimension.
% 'DataTransformationFunctionArgs'
%			- a struct containing any arguments required by the data transformation function.
%			One of the struct fields must be called 'Data', and will contain the data returned by
%			the NextDataFunction (i.e. the data to be transformed).
%
% 'ClusteringFunctionArgs'
%			- a struct containing any arguments required by the clustering function.
% 'FinalPassClusteringFunction'
%			- a string which stores the name of a function (M file) that is used in the
%			same way as 'ClusteringFunction', however this clustering function is only
%			called once to perform the final clustering of the set of representative
%			points. Since the function is only 	called once, it can be slower (and hopefully
%			more accurate) than 'ClusteringFunction'. For example, the non-iterative
%			k-means algorithm is 	less-accurate (but faster) than the iterative version.
%			Hence the non-iterative version could be specified for the main clustering
%			function, and the iterative version used for the final pass. If this function is
%			not specified, it is assumed that the user wants to use the same clustering
%			function during the final pass as during the main part of the algorithm.
% 'FinalPassClustFunctionArgs'
%			- a struct containing any arguments required by the clustering function used in
%			the final pass.
% 'MaxRunTime'
%			- the maximum time (in hours) to run for before moving on to perform the
%			final clustering. Defaults to 12.
%
% Return Value:
%
% MB_CLUSTER_LARGE_DATA_SET returns a final clustering of the data in the same form
% as the clustering function used in the final pass (see the 'FinalPassClusteringFunction'
% argument for more details). Therefore, at least a cell array where cell i contains a matrix of a
% sample of the data points assigned to cluster i will be returned, but other information may be
% useful (for example cluster means and covariances, the number of points assigned to each
% cluster etc.) -- if so, the clustering function used in the final pass of the data should be edited
% to return any additional information (as fields in the returned struct). The struct also contains
% a field called RepresentativeDataPoints, which contains the set of data points that was used
% to produce the final clustering. This is included in case a poor final clustering occurs.
%
% References:
%
% A. K. Jain, M. N. Murty and P. J. Flynn. Data Clustering: A Review. ACM Computing Surveys,
% vol 1, no. 3, Sept 1999.
%
% See Also:
%
% MB_K_MEANS_CLUSTERING

args = u_packargs(varargin,... % the user's input
			 'strict', ... % strict mode
			 {'k', 'NextDataFunction', 'NextDataFunctionArgs', 'ClusteringFunction',...
             'TempDirPath', 'MaxMemoryFinalClusteringInMB'}, ... % the mandatory arguments
			 'DataTransformationFunction', [],... % the optional arguments
			 'DataTransformationFunctionArgs', [],...
			 'ClusteringFunctionArgs', [],...
			 'FinalPassClusteringFunction', [],...
			 'FinalPassClustFunctionArgs', [],...
			 'MaxRunTime', 12);

% get the start time
start_time = clock;
			 
% check to see if a final pass clustering function was provided
if isempty(args.FinalPassClusteringFunction)
	% no, so assume they want to use the same clustering function for the final pass as the other passes. 
	args.FinalPassClusteringFunction = args.ClusteringFunction;
	args.FinalPassClustFunctionArgs = args.ClusteringFunctionArgs;
end

%Make sure temp dir exists
if ~isdir(args.TempDirPath)
    mkdir(args.TempDirPath)
end
% ensure that the TempDirPath ends in a filesep
if ~(args.TempDirPath(end) == filesep)
    args.TempDirPath = [args.TempDirPath filesep];
end

% initialisation/pre-allocation
% args.ClusteringFunctionArgs.k = args.k;
% more_data_flag = [];
% this_clustered_samples = [];
% cluster_sizes = [];
% clustered_samples = [];

% do the clustering
total_num_rep_points = 0;
rep_points_batch_counter = 0; % this counts the number of batches of representative points -- used in saving
data_dimensionality = []; % this will store then number of dims in our data
while 1
	% get the next chunk of data
    disp(datestr(now))
    disp('Large Data Clustering: Getting next chunk of data')
	[args.ClusteringFunctionArgs.Data, finished_sampling_flag] = feval(args.NextDataFunction, args.NextDataFunctionArgs);
	if finished_sampling_flag % if there is no more data to sample, then quit this loop
		break;
	end
    disp(datestr(now))
    disp(['Large Data Clustering: Sampled ' num2str(size(args.ClusteringFunctionArgs.Data,1)) ' data points'])
	
	% apply any transformation specified
	if ~isempty(args.DataTransformationFunction)
		args.DataTransformationFunctionArgs.Data = args.ClusteringFunctionArgs.Data;
		args.ClusteringFunctionArgs.Data = feval(args.DataTransformationFunction, args.DataTransformationFunctionArgs);
	end
	
    % determine if we have enough data point to do the clustering (at end of clustering, this might not be the case)
    % if we have pretty much finished sampling, then we won't cluster again until
    % the next data function reports to us that it has finished sampling -- and then we will skip this clustering loop.
    if (size(args.ClusteringFunctionArgs.Data,1) / args.ClusteringFunctionArgs.k) > (args.ClusteringFunctionArgs.SmallClusterSize * 1.3)
        % the 1.3 is a safety margin
        
        % cluster this data
        disp(datestr(now))
        disp('Large Data Clustering: About to cluster the data')
        clustering_result = feval(args.ClusteringFunction, args.ClusteringFunctionArgs);
        disp(datestr(now))
        disp(['Large Data Clustering: Clustering algorithm returned ' num2str(length(clustering_result.RepPoints)) ' clusters'])
        
        % make a note of the data dimensionality
        if isempty(data_dimensionality)
            data_dimensionality = size(clustering_result.RepPoints{1}, 2); % the number of dims in the data
        end
        
        % now save out the representative points from the clustering to the caller's specified directory
        this_batch_rep_points = cell2mat(clustering_result.RepPoints(:)); % get the cell array of rep points as a matrix, each row a vector
        total_num_rep_points = total_num_rep_points + size(this_batch_rep_points, 1); % update the number of representative points
        filename = [args.TempDirPath 'rep_points_batch' num2str(rep_points_batch_counter) '.mat'];
        save(filename, 'this_batch_rep_points');
        clear('this_batch_rep_points'); % remove the temp batch result from memory
        clear('clustering_result'); % remove the whole clustering result from memory
        
        disp(datestr(now))
        disp(['Large Data Clustering: The total number of representative points is now ' num2str(total_num_rep_points)])
        
        % increment the number of batches of rep points we've seen
        rep_points_batch_counter = rep_points_batch_counter + 1;
    end
	
	% see if we have been doing this long enough yet
	if etime(clock, start_time) > (args.MaxRunTime * 60 * 60) % convert the num of hours to num of seconds
		% warm the user about finishing early and not deleting the temp files
		disp(['mb_gmm_cluster_large_data_set has been running for ' num2str(etime(clock, start_time) / 60 / 60) ' hours']); % convert secs to hours
		disp('So we''re quitting early, as instructed.')
		
		args.NextDataFunctionArgs.CleanUp = (1==1); % we will want to clean up after ourselves
		feval(args.NextDataFunction, args.NextDataFunctionArgs); % perform the clean up 
		
		disp('Have deleted the temporary files used in clustering')
		
		break; % break out of the clustering loop
	end
end

% now we have a bunch of representative data from clusters found,
% perform the final clustering pass

% get a dir listing in the temp dir of all files matching the naming convention of our representative points
dir_list = dir([args.TempDirPath 'rep_points_batch*.mat']);

% work out how many data points we can fit in memory for the final clustering
mem_for_one_data_point_in_MB = (data_dimensionality * 8) / 1024 / 1024;
num_data_points_in_allocated_memory = args.MaxMemoryFinalClusteringInMB / mem_for_one_data_point_in_MB;

% work out how many data points to take from each file -- in case we need to sample from each file
num_data_points_from_each_file = floor(num_data_points_in_allocated_memory/ length(dir_list));

% add each set of rep points to the central pool
args.FinalPassClustFunctionArgs.Data = []; % ensure this is empty and defined
for i = 1 : length(dir_list)
	load([args.TempDirPath   dir_list(i).name]); % variable is called this_batch_rep_points
	% work out if the all the representative points amount to more memory than we have allocated
	if  total_num_rep_points > num_data_points_in_allocated_memory
		% we have to select a subset of the representative points to use in the final clustering
        if num_data_points_from_each_file > size(this_batch_rep_points,1)
            num_this_time = size(this_batch_rep_points,1);
        else
            num_this_time = num_data_points_from_each_file;
        end
		random_indices = randsample(1:size(this_batch_rep_points,1), num_this_time);
		this_batch_rep_points = this_batch_rep_points(random_indices, :); % select the sample
	end
	args.FinalPassClustFunctionArgs.Data = [args.FinalPassClustFunctionArgs.Data; this_batch_rep_points]; % collect together all the rep points into a matrix
end
clear('this_batch_rep_points'); % can do away with the copy of the last rep points file
% end 

disp(datestr(now))
disp(['Large Data Clustering: About to perform the final clustering on ' num2str(size(args.FinalPassClustFunctionArgs.Data,1)) ' data points'])

% Do the final clustering
args.FinalPassClustFunctionArgs.k = args.k;
clustering_result = feval(args.FinalPassClusteringFunction, args.FinalPassClustFunctionArgs);

% Save the clustering_result, in case we run out of memory in the last few stages
save([args.TempDirPath 'clustering_result.mat'], 'clustering_result');

% Get the data used in the final clustering, in case we got all this way and produced a bad result
clustering_result.RepresentativeDataPoints = args.FinalPassClustFunctionArgs.Data;

% the clustering has now been completed -- clear up after ourselves
delete([args.TempDirPath 'rep_points_batch*.mat']); % delete the files containing the representative points
% delete([args.TempDirPath 'clustering_result.mat']); % delete the file containing the clustering result
