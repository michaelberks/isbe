function [result] = mb_k_means_clustering2(varargin)
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
				'0',...
				{'Data', 'k'}, ... % the mandatory arguments
				'SmallClusterSize', 30,... % now the optional arguments and their defaults
				'PropRepresentativePoints', [],...
                ... the following parameters to be used by Matlab's kmeans
                'Display', 'Notify', ... 
				'SimilarityFunction', 'sqEuclidean', ...
				'Start', 'sample', ...
				'MaxIter', 1000,...
                'Replicates', 1,...
                'EmptyAction', 'drop');
clear varargin;

% CHECK DATA:
if size(args.Data,1) / args.k < args.SmallClusterSize
    disp('k is too big (or there are too few data points); unless otherwise specified, this algorithm discards')
    disp('clusters with fewer than 30 data points assigned to it. See the argument ''SmallClusterSize'' for details.')
    disp('Or, the files that keep a record of what we''ve sampled are still lying in the temp dir.')
	error('k is too big, there are too few data points, or temp files are lying around')
end
if ~isempty(args.PropRepresentativePoints)
    if args.PropRepresentativePoints < 0 || args.PropRepresentativePoints > 1
		error('PropRepresentativePoints must be between 0 and 1');
    end
    if args.PropRepresentativePoints * args.SmallClusterSize < 1
		error('PropRepresentativePoints * SmallClusterSize must be greater than 1');
    end
end

%Use inbuilt matlab clustering algorithm
[IDX, cluster_means, cluster_dists] = kmeans(args.Data, args.k, 'Display', args.Display,...
    'Distance', args.SimilarityFunction, 'Start', args.Start, ...
    'MaxIter', args.MaxIter, 'Replicates', args.Replicates, ...
    'EmptyAction', args.EmptyAction);
k = 0; %k is count of non-discarded clusters
total_pts = 0; %count of the total samples we have

if ~isempty(args.PropRepresentativePoints)
    result.RepPoints = [];
end
for ii = 1:args.k
    %find out which rows of data belong to current cluster
    curr_idx = IDX == ii;
    % Check to see we have enough points 
    curr_num_pts = sum(curr_idx);
    if curr_num_pts >= args.SmallClusterSize
        % if we're keeping this cluster...
        k = k+1; % update count of clusters
        total_pts = total_pts + curr_num_pts; % update count of samples
        result.ClusterProbs(k) = curr_num_pts; % record number of points in cluster (will normalise later)
        curr_data = args.Data(curr_idx,:); % extract samples
        result.CovMats{k} = cov(curr_data); % compute covariance for cluster
        result.Means(k, :) = cluster_means(k, :); % get mean for cluster (pre-computed by Matlab kmeans)
        
        % check to see if we're recording some representative points
        if ~isempty(args.PropRepresentativePoints)
            % Compute number of points for this cluster - 
            % we know num_rep_points >= 1, since we've checked
            % prop_rep_points*small_cluster_size > 1
            num_rep_points = floor(args.PropRepresentativePoints * curr_num_pts);
            
            % Take random sample from data 
            random_indices = randsample(curr_num_pts, num_rep_points);
            result.RepPoints = [result.RepPoints; curr_data(random_indices, :)];
        end       
    end
end
result.N_S = total_pts; % Record total samples kept
result.NumClusters = k; % Record number of clusters kept
result.ClusterDists = cluster_dists;
% normalise the ClusterProbs so that they sum to unity
result.ClusterProbs = result.ClusterProbs ./ sum(result.ClusterProbs);
