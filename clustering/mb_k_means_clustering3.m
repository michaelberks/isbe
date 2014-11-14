function [result] = mb_k_means_clustering3(varargin)
%MB_K_MEANS_CLUSTERING Perform k-means clustering.
%
% The function allows both the iterative and non-iterative versions of the algorithm to be
% used. The cluster means and covariance matrices are returned, along with other useful information.
% The function allows any similarity function to be used to measure the similarity
% between data points (so a measure meaningful to the application can be used).
% The function discards clusters with few data points assigned to them, assuming
% that such clusters are unrepresentative of the data (see the argument 'SmallClusterSize'
% for more details). The function can return a number of representative points from each
% cluster, which may be useful when using this function as part of a more sophisticated
% clustering scheme.
%
% MB_K_MEANS_CLUSTERING uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%	'Data'
%			- the data to cluster where a row is a vector (i.e. a data point), and
%			the columns are the dimensions of the data.
%	'k' 		
%			- the number of prototypes to find -- although the actual number of
%			clusters may be less (see the 'SmallClusterSize' argument).
%	'Iterative'
%			- a string that specifies whether the iterative or non-iterative version of the k-means
%			algorithm should be used. The string should be 'iterate' for the iterative version
%			or 'noiterate' for the non-iterative version. The iterative version is likely to give a
%			better clustering, but is slower (approx 10 times slower than the non-iterative one,
%			but this is a very rough guide and depands on the stopping criterea; see the
%			'StoppingPercentage' argument).
%
% Optional Arguments:
%	'SimilarityFunction'
%			- a function handle to a function that measure the similarity betwen
%			two data points in a meaninful way to your application. The
%			function must be specified as follows:
%				function [sim] = my_similarity_masure(x1, x2, params_struct)
%			where params_struct is a structure of any parameters that the
%			similarity function may need. This can be passed via
%			'SimilarityFunctionArgs' -- see below for more details.
%			The similarity measure sim must be small for similar vectors and large for
%			dissimilar ones; ideally the similarity function should be a metric, but this
%			may not be too important. If a similarity function is not specified, then the
%			Euclidean metric is used by default.
%	'SimilarityFunctionArgs'
%			- a struct of any arguments required by the similarity function being used.
%	'PropRepresentativePoints'
%			- this algorithm can return a set of representative points from each cluster.
%			Use this argument to specify the proportion to return from each cluster. I.e. it
%			MUST be a number from 0 to 1. If specified, the struct returned by this function
%			will have the field 'RepPoints' which will be a cell array, where RepPoints{i} will
%			be a matrix of the representative points of cluster i (where a row is a data point).
%			The representative points returned in RepPoints{i} are uniformly sampled from the
%			points assigned to cluster i. A proportion of points is chosen so that the density of the
%			representative point space should be approx. the same as the original data space
%			(i.e. if the number of representative points for all clusters was the same, the statistics
%			of the representaqtive point space would not be anywhere near those of the original
%			data space).
%
%	'SmallClusterSize'
%			- a scalar that describes the minimum size for a cluster; small clusters are
%			discarded. The default value for this is 30 data points.
%	'StoppingPercentage'
%			- if the iterative version of the algorithm is selected, then iterating stops when
%			the proportion of reassigned data points drops below this number (expressed as
%			a percentage). If this argument is not specified, then the default is 5 (i.e. 5%).
%	'DebugFlag'
%			- if set to 1, then debug info is shown (currently for the iterative version only)
%
% MB_K_MEANS_CLUSTERING returns a struct with the following fields:
%	'NumClusters'
%			- the number of clusters that were actually formed (this may differ from k
%			since this function discards small clusters)
%	'N_S'
%			- the total number of data points that were used in clustering. This may not be
%			the same as the number of input data points, since small clusters are discarded.
%	'ClusterProbs'
%			- a vector where ClusterProbs(i) is computed from the number of data points
%			assigned to cluster i, and represents the probability of each cluster.
%	'Means'
%			- a matrix where Means(i,:) is the mean vector for cluster i
%	'CovMats'
%			- a cell array where CovMats{i} returns the covariance matrix for cluster i
%	'RepPoints'
%			- only included if the optional argument 'PropRepresentativePoints' was set.
%			Is a cell array where cell i is a matrix containing a number (specified by
%			'PropRepresentativePoints') of representative points from cluster i (the points
%			are sampled uniformly from each cluster). Only non-small clusters are included for
%			this sampling (see the argument 'SmallClusterSize'). 
%
% References:
%	E. Forgey. Cluster analysis of multivariate data; Efficiency vs. interpretability of classifications
%	. Biometrics, 21(768), 1965.
%	
%	J. MacQueen. Some methods for classification and analysis of multivariate observations. In
%	L. M. Le Cam and J. Neyman, editors, Proc. 5th Berkley Symposium on Mathematical Statistics
%	and Probability, vol 1 of Statistics. University of California Press, 1967.
%	
%	A. K. Jain, M. N. Murty and P. J. Flynn. Data Clustering: A Review. ACM Computing Surveys,
%	vol 1, no. 3, Sept 1999.
%
% See Also:
%	MB_CLUSTER_LARGE_DATA_SET

% pack the args
args = u_packargs(varargin,...
				'strict',...
				{'Data', 'k', 'Iterative'}, ... % the mandatory arguments
				'DebugFlag', 0, ... % now the optional arguments and their defaults
				'SimilarityFunction', @euclidean_metric, ...
				'SimilarityFunctionArgs', [], ...
				'SmallClusterSize', 30,...
				'StoppingPercentage', 0,...
				'PropRepresentativePoints', []);

% CHECK DATA:
if size(args.Data,1) / args.k < args.SmallClusterSize
    disp('k is too big (or there are too few data points); unless otherwise specified, this algorithm discards')
    disp('clusters with fewer than 30 data points assigned to it. See the argument ''SmallClusterSize'' for details.')
    disp('Or, the files that keep a record of what we''ve sampled are still lying in the temp dir.')
	error('k is too big, there are too few data points, or temp files are lying around')
end
if ~isempty(args.PropRepresentativePoints)
	if args.PropRepresentativePoints < 0 || args.PropRepresentativePoints > 1
		error('PropRepresentativePoints must be between 0 and 1')
	end
end

% INITIALISATION STEPS:

% choose k initial estimates of the prototypes randomly from the data, without replacement
random_indices = randsample(size(args.Data,1), args.k);

C = args.Data(random_indices, :);

% DO THE CLUSTERING:

if args.Iterative
    result = do_iterative(args.Data, args.k, args.SimilarityFunction, C, ...
    						args.SimilarityFunctionArgs, args.SmallClusterSize,...
                            args.PropRepresentativePoints, args.StoppingPercentage);
else
	result = do_non_iterative(args.Data, args.k, args.SimilarityFunction, C, ...
    						args.SimilarityFunctionArgs, args.SmallClusterSize,...
                            args.PropRepresentativePoints);
end

% normalise the ClusterProbs so that they sum to unity
result.ClusterProbs = result.ClusterProbs ./ sum(result.ClusterProbs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [result] = do_non_iterative(data, k, similarity_function, C_means, ...
								similarity_function_args, small_cluster_size, prop_rep_points)
%
% Perform the non-iterative method of clustering
% loop over all the data, and assign it to it's closest cluster,
% update the cluster estimate.

n_pts = zeros(k); %running count of the points assigned to each cluster
c_idx = zeros(size(data, 1), 1); %cluster assigned to each point
for i = 1 : size(data, 1)
%     % get this data point
    
    % find nearest cluster mean
    % repmat sum is quicker than looping through each cluster mean
    this_data = repmat(data(i,:), k, 1);
    dists = feval(similarity_function, this_data, C_means, similarity_function_args);
    [dummy closest_cluster] = min(dists);
    
    % assign point to its closest cluster and increment cluster's count
    c_idx(i) = closest_cluster;
    n_pts(closest_cluster) = n_pts(closest_cluster) + 1;
    
    % update this cluster's mean
    C_means(closest_cluster,:) = ((n_pts(closest_cluster)-1)*C_means(closest_cluster,:) + ....
        this_data(1,:)) / n_pts(closest_cluster);
    clear this_data;
end

% Go through clusters, check if we're keeping them, caluclate the covs etc.
retained_k = 0;
total_pts = 0;
for i = 1 : k
    curr_idx = c_idx == i;
    % Check to see we have enough points 
    curr_num_pts = n_pts(i);
    if curr_num_pts >= small_cluster_size
        % if we're keeping this cluster...
        retained_k = retained_k+1; % update count of clusters
        total_pts = total_pts + curr_num_pts; % update count of samples
        result.ClusterProbs(retained_k) = curr_num_pts; % record number of points in cluster (will normalise later)
        curr_data = data(curr_idx,:); % extract samples
        result.CovMats{retained_k} = cov(curr_data); % compute covariance for cluster
        result.Means(retained_k, :) = C_means(retained_k, :); % get mean for cluster (pre-computed by Matlab kmeans)
        
        % check to see if we're recording some representative points
        if ~isempty(prop_rep_points)
            % Compute number of points for this cluster - 
            % we know num_rep_points >= 1, since we've checked
            % prop_rep_points*small_cluster_size > 1
            num_rep_points = floor(prop_rep_points * curr_num_pts);
            
            % Take random sample from data 
            random_indices = randsample(curr_num_pts, num_rep_points);
            result.RepPoints{retained_k} = curr_data(random_indices, :);
        end       
    end
end

result.NumClusters = retained_k;
result.N_S = total_pts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [result] = do_iterative(data, k, similarity_function, C_means, ...
							similarity_function_args, small_cluster_size,...
                            prop_rep_points, stopping_percentage)
%
% Perform the iterative method of clustering

% get the number of item of data we have and the number of dimensions in the data
[num_data_items] = size(data, 1);

% initialisation step -- get an initial assignment of data points to clusters:

% loop over all the data, and assign it to it's closest cluster,
% update the cluster estimate.
n_pts = zeros(k); %running count of the points assigned to each cluster
c_idx = cell(k, 1);%cluster assigned to each point
for i = 1 : size(data, 1)
%     % get this data point
    
    % find nearest cluster mean
    % repmat sum is quicker than looping through each cluster mean
    this_data = repmat(data(i,:), k, 1);
    dists = feval(similarity_function, this_data, C_means, similarity_function_args);
    [dummy closest_cluster] = min(dists);
    
    % assign point to its closest cluster and increment cluster's count
    c_idx{i, 1} = [c_idx{i, 1}, closest_cluster];
    n_pts(closest_cluster) = n_pts(closest_cluster) + 1;
    
    % update this cluster's mean
    C_means(closest_cluster,:) = ((n_pts(closest_cluster)-1)*C_means(closest_cluster,:) + ....
        this_data(1,:)) / n_pts(closest_cluster);
    clear this_data;
end

% this is the iterative bit:

% look through all the data and work out which data points to re-assign,
% terminate if the number of data points re-assigned is less than 5% of the total number of
% data points

number_of_iterations = 0;
while 1==1
    % record how many points we've re-assigned
    total_points_moved = 0;
    cluster_points_moved = zeros(k,1);
    for cluster_counter = 1 : k % loop over the clusters
        
        data_point_counter = 1;
        
        while data_point_counter <= n_pts(cluster_counter); % this loops over the data points assigned to cluster number cluster_counter
            % get this data point
            this_data_point = data(c_idx{cluster_counter, 1}(data_point_counter),:);
            
            % find the cluster closest to this_data_point
            % repmat sum is quicker than looping through each cluster mean
            dists = feval(similarity_function, repmat(this_data_point, k, 1), C_means, similarity_function_args);
            [dummy closest_cluster] = min(dists);
 
            % if the closest cluster is not the cluster the data point is already in then remove the data point from where it was
            % and place it in the closest_cluster, and then incrememnt a count of how many data points we've moved
            if closest_cluster ~= cluster_counter
                
                % update the cluster index lists
                c_idx{closest_cluster, 1}(data_point_counter) = ...
                    [c_idx{closest_cluster, 1} c_idx{cluster_counter, 1}(data_point_counter)];
                c_idx{cluster_counter, 1}(data_point_counter) = [];
                
                % update the cluster counts
                n_pts(cluster_counter) = n_pts(cluster_counter) - 1;
                n_pts(closest_cluster) = n_pts(closest_cluster) + 1;
                
                % increment the number of points re-assigned
                cluster_points_moved(cluster_counter) =...
                    cluster_points_moved(cluster_counter) + 1;
                
                %but don't increment the data_point_counter
            else
                % move to next point by incrementing the data_point_counter
                data_point_counter = data_point_counter + 1;
            end
        end
        total_points_moved = total_points_moved +...
            cluster_points_moved(cluster_counter);
    end
    
    % update all the cluster means
    for i = 1 : k
        %No need to update if no points moved 
        if cluster_points_moved(i)
        % only update the mean if there is at least one data point assigned
        % to the cluster, else choose a new value for the mean at random
            if n_pts(i) >= 1
                C_means(i, :) = mean(data(c_idx{i, 1},:), 1);
            else
                C_means = data(randsample(1:num_data_items, 1),:);
            end
        end
    end
    
    % increment the number of iterations
    number_of_iterations = number_of_iterations + 1;
    
    % work out if we should continue iterating or stop
    if number_of_points_moved <= ceil(stopping_percentage*num_data_items /100)
        break; % stop now
    end
end

% Go through clusters, check if we're keeping them, caluclate the covs etc.
retained_k = 0;
total_pts = 0;
for i = 1 : k
    curr_idx = c_idx == i;
    % Check to see we have enough points 
    curr_num_pts = n_pts(i);
    if curr_num_pts >= small_cluster_size
        % if we're keeping this cluster...
        retained_k = retained_k+1; % update count of clusters
        total_pts = total_pts + curr_num_pts; % update count of samples
        result.ClusterProbs(retained_k) = curr_num_pts; % record number of points in cluster (will normalise later)
        curr_data = data(curr_idx,:); % extract samples
        result.CovMats{retained_k} = cov(curr_data); % compute covariance for cluster
        result.Means(retained_k, :) = C_means(retained_k, :); % get mean for cluster (pre-computed by Matlab kmeans)
        
        % check to see if we're recording some representative points
        if ~isempty(prop_rep_points)
            % Compute number of points for this cluster - 
            % we know num_rep_points >= 1, since we've checked
            % prop_rep_points*small_cluster_size > 1
            num_rep_points = floor(prop_rep_points * curr_num_pts);
            
            % Take random sample from data 
            random_indices = randsample(curr_num_pts, num_rep_points);
            result.RepPoints{retained_k} = curr_data(random_indices, :);
        end       
    end
end

% now build all these into a struct:
result.NumClusters = retained_k;
result.N_S = total_pts; % N_S means number of samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dist] = euclidean_metric(x1, x2, varargin)
%
% Return the Euclidean distance between x1 and x2

square_diffs = (x1 - x2) .^ 2;
dist = sqrt(sum(square_diffs, 2));