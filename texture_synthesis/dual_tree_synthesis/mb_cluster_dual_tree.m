function mb_cluster_dual_tree(varargin)

% Clustering method based on the divide and conquer method implemented by
% CJR. This function implements the final gathering of the divided cluster
% results and performs the final clustering

% get a dir listing in the results dir of files matching the naming
% convention of our cluster results

args = u_packargs(varargin,... % the user's input
			 '0', ... % strict mode
			 {'ClusteringFunctionArgs',...
             'DualTreeDir',...
             'FinalFile'}, ... % the mandatory arguments
			 'ClusteringFunction', 'mb_k_means_clustering2',...
             'MaxFinalMemory', 512,...
             'NumLevels', 4, ...
             'Quiet', 0 ...
         );
     
clear varargin;

if ~args.Quiet
    display('version: 30/01/2008');
end

%get cluster indices
args.ClusteringFunctionArgs.Data = ...
    get_cluster_data(args.DualTreeDir, args.MaxFinalMemory, args.NumLevels);

%don't want to create large matrices so do standardisation column by column
stds = std(args.ClusteringFunctionArgs.Data);
%means = mean(args.ClusteringFunctionArgs.Data);

for dim = 1:length(stds);
    args.ClusteringFunctionArgs.Data(:,dim) =...
        args.ClusteringFunctionArgs.Data(:,dim) / stds(dim);
        %(args.ClusteringFunctionArgs.Data(:,dim) - means(dim)) ./ stds(dim);
end

if ~args.Quiet
    display('Cluster data sampled');
    display(['About to perform clustering on ',...
        num2str(size(args.ClusteringFunctionArgs.Data,1)), ' points']);
end

% Do the final clustering
cluster_model = feval(args.ClusteringFunction, args.ClusteringFunctionArgs); %#ok

% Save the clustering_result
save(args.FinalFile, 'cluster_model');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = get_cluster_data(dual_tree_dir, max_memory, num_levels)

% Set constant
size_of_data_element = 8; % assume each data element is a double, hence 8 bytes per element.

if dual_tree_dir(end) ~= filesep
	dual_tree_dir = [dual_tree_dir filesep];
end

% first get a directory listing and other info about the image directory
dual_tree_files = dir([dual_tree_dir,'*dual_tree*']);
num_trees = length(dual_tree_files);

% work out how much memory (in bytes) each window will occupy
mem_per_sample = 2* num_levels * size_of_data_element;

% work out how many points we can fit in max_memory
num_pts = floor((max_memory * 1024 * 1024) /  mem_per_sample);

% work out how many points per image this is
num_pts_per_tree = floor(num_pts / num_trees);

if num_pts_per_tree < 1
	error('This function assumes that we can sample at least one window from each image')
end

display(['Next Data Function: num_points_per_image = ' num2str(num_pts_per_tree)])

%pre-allocate for maximum amount of data
data = repmat(NaN, num_pts, 2*num_levels);

%Start counter
curr_idx = 1;

for ii = 1:num_trees
	
    % open the i-th DT file
	dual_tree = u_load([dual_tree_dir, dual_tree_files(ii).name]);
	
    %Work out how many pts to take from this tree
    [rows cols] = size(dual_tree{1}(:,:,1));
    pts_in_tree = rows*cols;
    curr_pts = min(num_pts_per_tree, pts_in_tree);
    
    pts = zeros(curr_pts, 2);
    
    %Randomly sample indices and convert [r c] pts
    [pts(:,1) pts(:,2)] = ind2sub([rows, cols], randsample(1:pts_in_tree, curr_pts));
    
    %Get DT feature vectors for these points
    data(curr_idx:curr_idx + curr_pts - 1, :) = ...
        mb_get_dt_vector(pts, dual_tree, num_levels);

    %Update counter to new position
    curr_idx = curr_idx + curr_pts;    
end

%clear up any space in data we haven't assigned (e.g. if less points in
%dual-trees than maximum size of data specified)
data(curr_idx:end,:) = [];


end



