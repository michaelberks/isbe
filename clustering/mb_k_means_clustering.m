function [result] = mb_k_means_clustering(varargin)
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
				{'Data', 'k', 'Iterative'}, ... % the mandatory arguments
				'DebugFlag', 0, ... % now the optional arguments and their defaults
				'SimilarityFunction', @euclidean_metric, ...
				'SimilarityFunctionArgs', [], ...
				'SmallClusterSize', 30,...
				'StoppingPercentage', 5,...
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

% GET SOME USEFUL INFO:

% get the number of item of data we have and the number of dimensions in the data
[num_data_items num_dims] = size(args.Data);

% INITIALISATION STEPS:

% make a cell array to hold the prototypes and the data assigned to each
% prototype. C is a cell array with k rows and 2 columns.
% C{n,1} contains the n-th prototype (a row vector).
% C{n,2} contains a matrix of row vectors, where each row
% is the data item assigned to the n-th prototype.
C = cell(args.k,2);

% initialise the prototypes

% choose k initial estimates of the prototypes randomly from the data, without replacement
random_indices = randsample(num_data_items, args.k);

% get the initial estimates of the prototypes
for i = 1 : args.k
    C{i,1} = args.Data(random_indices(i), :);
end

% DO THE CLUSTERING:

if args.Iterative
    result = do_iterative(args.Data, args.k, args.SimilarityFunction, C, num_data_items, num_dims, args.DebugFlag, ...
    						args.SimilarityFunctionArgs, args.SmallClusterSize, args.PropRepresentativePoints);
else
	result = do_non_iterative(args.Data, args.k, args.SimilarityFunction, C, num_data_items, num_dims, ...
    						args.SimilarityFunctionArgs, args.SmallClusterSize, args.PropRepresentativePoints);
end

% normalise the ClusterProbs so that they sum to unity
result.ClusterProbs = result.ClusterProbs ./ sum(result.ClusterProbs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [result] = do_non_iterative(data, k, similarity_function, C, num_data_items, num_dims, ...
								similarity_function_args, small_cluster_size, prop_rep_points)
%
% Perform the non-iterative method of clustering

% loop over all the data, and assign it to it's closest cluster,
% update the cluster estimate.
for i = 1 : num_data_items
    % get this data point
    this_data = data(i,:);
    
    % find which cluster the data point is closest to
    closest_cluster = 1; % assume it's closest to cluster 1 for the moment
    for j = 2 : k
        if feval(similarity_function, this_data, C{j,1}, similarity_function_args) < feval(similarity_function, this_data, C{closest_cluster,1}, similarity_function_args)
            closest_cluster = j;
        end
    end
    
    % now assign the data point to its closest cluster
    C{closest_cluster, 2} = [C{closest_cluster, 2} ; this_data];
    
    % update this cluster's mean -- the sum of the dims of the data points assigned to the cluster, divided by the number of data points assigned
    C{closest_cluster, 1} = sum(C{closest_cluster, 2},1) / size(C{closest_cluster, 2}, 1);
end

% pre-allocate clusters
clusters = zeros(k, num_dims);
num_dat_points_assigned = zeros(1, k);

for i = 1 : k
    % build clusters
    clusters(i,:) = C{i,1};
    % build the covariance matrices
    CM{i} = cov(C{i,2});
    % build the vector of the number of data points assigned to each
    % cluster
    num_dat_points_assigned(i) = size(C{i,2},1);
end

% now build all these into a struct:
result.NumClusters = k;
result.N_S = num_data_items; % N_S means number of samples
result.ClusterProbs = num_dat_points_assigned;
result.Means = clusters;
result.CovMats = CM;

% discard small clusters (will update NumClusters)
[result, included_clusters] = discard_small_clusters(result, k, small_cluster_size);

% get the representative points for each cluster, if that's required, and assign it to result.RepPoints
if ~isempty(prop_rep_points)
	result.RepPoints = get_representative_points(prop_rep_points, C, included_clusters, result.NumClusters, k); 
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [result] = do_iterative(data, k, similarity_function, C, num_data_items, num_dims,...
							debug_flag, similarity_function_args, small_cluster_size, prop_rep_points)
%
% Perform the iterative method of clustering

% initialisation step -- get an initial assignment of data points to clusters:

% loop over all the data, and assign it to it's closest cluster,
% update the cluster estimate.
for i = 1 : num_data_items
    % get this data point
    this_data = data(i,:);
    
    % find which cluster the data point is closest to
    closest_cluster = 1; % assume it's closest to cluster 1 for the moment
    for j = 2 : k
        if feval(similarity_function, this_data, C{j,1}, similarity_function_args) < feval(similarity_function, this_data, C{closest_cluster,1}, similarity_function_args)
            closest_cluster = j;
        end
    end
    
    % now assign the data point to its closest cluster
    C{closest_cluster, 2} = [C{closest_cluster, 2} ; this_data];
    
    % update this cluster's mean -- the sum of the dims of the data points assigned to the cluster, divided by the number of data points assigned
    C{closest_cluster, 1} = sum(C{closest_cluster, 2},1) / size(C{closest_cluster, 2}, 1);
    
end

%update all the cluster means -- the sum of the dims of the data points assigned to the cluster, divided by the number of data points assigned
for i = 1 : k
    C{i, 1} = sum(C{i, 2},1) / size(C{i, 2}, 1);
end

% pre-allocate clusters
% orig_clusters = zeros(k, num_dims);

% % build clusters
% for i = 1 : k
%     orig_clusters(i,:) = C{i,1};
% end

% this is the iterative bit:

% look through all the data and work out which data points to re-assign,
% terminate if the number of data points re-assigned is less than 5% of the total number of
% data points

number_of_iterations = 0;
while 1==1
    % record how many points we've re-assigned
    number_of_points_moved = 0;
    for cluster_counter = 1 : k % loop over the clusters
        num_data_points_in_cluster = size(C{cluster_counter, 2},1);
        data_point_counter = 1;
        while data_point_counter <= num_data_points_in_cluster % this loops over the data points assigned to cluster number cluster_counter
            % get this data point
            this_clusters_data = C{cluster_counter, 2};
            this_data_point = this_clusters_data(data_point_counter, :);
            
            % find the cluster closest to this_data_point
            closest_cluster = cluster_counter; % assume it's closest to the cluster it's in for the moment
            for j = 1 : k
                % get the cluster mean that the data point belongs to
                cluster_mean = C{cluster_counter,1};
                % get the cluster mean of the next candidate cluster
                candidate_cluster_mean = C{j,1};
                % work out the distance between the data point and the cluster mean
                current_cluster_dist = feval(similarity_function, this_data_point, cluster_mean, similarity_function_args);
                % work out the distance between the data point and the next cluster mean
                potential_new_cluster_dist = feval(similarity_function, this_data_point, candidate_cluster_mean, similarity_function_args);
                % work out if the data point is closer to the cluster
                if potential_new_cluster_dist < current_cluster_dist
                    closest_cluster = j;
                end
            end
            
            
            
            % if the closest cluster is not the cluster the data point is already in then remove the data point from where it was
            % and place it in the closest_cluster, and then incrememnt a count of how many data points we've moved
            if closest_cluster ~= cluster_counter
                % get a copy of the data points assigned to cluster we'll be removing from
                cluster_data_points = C{cluster_counter, 2};
                % remove the data point
                cluster_data_points(data_point_counter, :) = [];
                % put the cluster data points back into the cluster
                C{cluster_counter, 2} = cluster_data_points;
                % update the number of data points in this cluster
                num_data_points_in_cluster = num_data_points_in_cluster - 1;
                % append this_data_point to the cluster it was closest to
                C{closest_cluster, 2} = [C{closest_cluster, 2}; this_data_point];
                % increment the number of points re-assigned
                number_of_points_moved = number_of_points_moved + 1;
            end
            
            % increment the data_point_counter
            data_point_counter = data_point_counter + 1;
        end
    end
    
    % update all the cluster means -- the sum of the dims of the data points assigned to the cluster, divided by the number of data points assigned
    for i = 1 : k
        % only update the mean if there is at least one data point assigned to the cluster, else choose a new value for the mean at random
        if size(C{i, 2}, 1) >= 1
            C{i, 1} = sum(C{i, 2},1) / size(C{i, 2}, 1);
        else
            C{i, 1} = data(randsample(1:num_data_items, 1),:);
        end
    end
    
    % increment the number of iterations
    number_of_iterations = number_of_iterations + 1;
    
    % work out if we should continue iterating or stop
    if number_of_points_moved == 0
    %if number_of_points_moved < ceil((num_data_items /100) * 5)
        break; % stop now
    end
end

% pre-allocate clusters
clusters = zeros(k, num_dims);
num_dat_points_assigned = zeros(1, k);

for i = 1 : k
    % build clusters
    clusters(i,:) = C{i,1};
    % build the covariance matrices
    CM{i} = cov(C{i,2});
    % build the vector of the number of data points assigned to each
    % cluster
    num_dat_points_assigned(i) = size(C{i,2},1);
end

% now build all these into a struct:
result.NumClusters = k;
result.N_S = num_data_items; % N_S means number of samples
result.ClusterProbs = num_dat_points_assigned;
result.Means = clusters;
result.CovMats = CM;

% discard small clusters
[result, included_clusters] = discard_small_clusters(result, k, small_cluster_size);

% get the representative points for each cluster, if that's required, and assign it to result.RepPoints
if ~isempty(prop_rep_points)
	result.RepPoints = get_representative_points(prop_rep_points, C, included_clusters, result.NumClusters, k); 
end

% now do any debug
if debug_flag == 1
    error('This needs to be fixed to work on result instead of clusters')
    % plot the ata points in a colour that indicates their cluster membership
%     figure;
%     hold on;
%     colour_index = 1;
%     for i = 1 : k
%         % get the next cluster's data points
%         this_debug_data = C{i,2};
%         
%         % plot there data points in the current colour
%         switch colour_index
%         case 1
%             plot_string = 'b.';
%         case 2
%             plot_string = 'g.';
%         case 3
%             plot_string = 'r.';
%         case 4
%             plot_string = 'c.';
%         case 5
%             plot_string = 'm.';
%         case 6
%             plot_string = 'y.';
%         case 7
%             plot_string = 'k.';
%         end
%         plot(this_debug_data(:,1), this_debug_data(:,2), plot_string);
%         
% %         % now plot the best initial cluster candidate
% %         plot(orig_clusters(:,1), orig_clusters(:,2) , 'kd');
%         
%         % now plot the cluster means
%         plot(clusters(:,1), clusters(:,2) , 'kx');
%         
%         % increment the colour_index
%         colour_index = colour_index + 1;
%         if colour_index > 7
%             colour_index = 1;
%         end
%     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [result, included_clusters] = discard_small_clusters(temp_result, k, small_cluster_size)
%
% now discard the clusters that have less than 30 data points assigned to them
%
% included_clusters records which of the original k clusters were included

counter = 1;
included_clusters = [];
for i = 1 : k
    if temp_result.ClusterProbs(i) >= small_cluster_size % include this cluster
    	included_clusters = [included_clusters ; i];
        result.ClusterProbs(counter) = temp_result.ClusterProbs(i);
        result.Means(counter, :) = temp_result.Means(i,:);
        result.CovMats{counter} = temp_result.CovMats{i};
        counter  = counter + 1;
    end
end
result.N_S = sum(result.ClusterProbs);
result.NumClusters = size(result.Means,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rep_points] = get_representative_points(prop_rep_points, C, included_clusters, number_of_clusters, k)
%
% Gets a proportion prop_rep_points of the points assigned poins from each cluster

rep_points = cell(number_of_clusters,1);
actual_cluster_counter = 1;
for i = 1 : k
	if any(i == included_clusters(:)) % if we didn't dicard this cluster
		% work out how many representative points the fraction prop_rep_points corresponds to
		num_rep_points = floor(prop_rep_points * size(C{i,2},1));
		if num_rep_points > 0 % ignore clusters which are too small to return at least one point
			random_indices = randsample(1:size(C{i,2},1), num_rep_points);
			this_random_sample = C{i,2}(random_indices, :); % 
			rep_points{actual_cluster_counter} = this_random_sample;
			actual_cluster_counter = actual_cluster_counter + 1;
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dist] = euclidean_metric(x1, x2, varargin)
%
% Return the Euclidean distance between x1 and x2

square_diffs = (x1 - x2) .^ 2;
dist = sqrt(sum(square_diffs));

