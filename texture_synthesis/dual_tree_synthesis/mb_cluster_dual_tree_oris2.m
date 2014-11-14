function mb_cluster_dual_tree_oris2(varargin)

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
             'WindowSize', 11,...
             'Level', 2, ...
             'Quiet', 0,...
             'Overlap', 20,...
             'WindowSize2', 0,...
             'Real', 1, ...
             'Standardise', 1);
clear varargin;

if ~args.Quiet
    display('version: 30/01/2008');
end

%get cluster indices
cluster_idx = get_cluster_idx(...
    args.DualTreeDir, args.Level, ...
    args.WindowSize, args.MaxFinalMemory, args.Overlap, args.WindowSize2);

% save the indices as they may be useful for debugging
save([args.FinalFile, '_data_idx'], 'cluster_idx');

if ~args.Quiet
    display('Cluster indices generated');
end

%set up args struct for getting data from indices
data_args.LoadIdx = 0;
data_args.Idx = cluster_idx; clear cluster_idx;
data_args.DualTreeDir = args.DualTreeDir;
data_args.Level = args.Level;
data_args.WindowSize = args.WindowSize;
data_args.WindowSize2 = args.WindowSize2;
data_args.Standardise = args.Standardise;

% Get data for orientations 1-6, with either {1,2,3} or {4,5,6} real
data_args.Real = args.Real;
data_args.Orientation = 1; %
[data1] =...
    mb_get_dual_tree_subband_data(data_args);
data_args.Orientation = 2; %
[data2] =...
    mb_get_dual_tree_subband_data(data_args);
data_args.Orientation = 3; %
[data3] =...
    mb_get_dual_tree_subband_data(data_args);

data_args.Real = ~args.Real;
data_args.Orientation = 4; %
[data4] =...
    mb_get_dual_tree_subband_data(data_args);
data_args.Orientation = 5; %
[data5] =...
    mb_get_dual_tree_subband_data(data_args);
data_args.Orientation = 6; %
[data6] =...
    mb_get_dual_tree_subband_data(data_args);

args.ClusteringFunctionArgs.Data = [data1 data2 data3 data4 data5 data6];
clear data1 data2 data3

if ~args.Quiet
    display('Data sampled');
    display(['About to perform clustering on ',...
        num2str(size(args.ClusteringFunctionArgs.Data,1)), ' points']);
end

% Do the final clustering
cluster_model = feval(args.ClusteringFunction, args.ClusteringFunctionArgs); %#ok

%Save stats into results structure
cluster_model.Stats = [];
cluster_model.WindowSize = args.WindowSize;
cluster_model.WindowSize2 = args.WindowSize2;

% Save the clustering_result
save(args.FinalFile, 'cluster_model');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function all_indices = get_cluster_idx(...
    DualTreeDir, Level, WindowSize, MaxMemory, overlap, WindowSize2)

% Set constant
size_of_data_element = 8; % assume each data element is a double, hence 8 bytes per element.


if DualTreeDir(end) ~= filesep
	DualTreeDir = [DualTreeDir filesep];
end

% first get a directory listing and other info about the image directory
DualTreeFiles = dir([DualTreeDir,'*dual_tree*']);
number_of_trees = length(DualTreeFiles);

% work out how much memory (in bytes) each window will occupy
mem_per_sample = 6*(WindowSize^2 + WindowSize2^2) * size_of_data_element;

% work out how many points we can fit in MaxMemory
num_points = floor((MaxMemory * 1024 * 1024) /  mem_per_sample);

% work out how many points per image this is
num_points_per_tree = floor(num_points / number_of_trees);

if num_points_per_tree < 1
	error('This function assumes that we can sample at least one window from each image')
end

display(['Next Data Function: num_points_per_image = ' num2str(num_points_per_tree)])

all_indices = repmat(NaN, number_of_trees, num_points_per_tree);

for ii = 1:number_of_trees
	
    % open the i-th image file
	dual_tree = u_load([DualTreeDir, DualTreeFiles(ii).name]);
	
	% create a matrix of zeros of the same size as the image
	sampled = zeros(size(dual_tree{Level}(:,:,1)));
    
    % set the edges to ones so that we don't sample windows that lie off
    % the edge of the window
    sampled(1:overlap, :) = 1; % top portion
    sampled(end-overlap+1:end, :) = 1; % bottom portion
    sampled(:, 1:overlap) = 1; % left side
    sampled(:,end-overlap+1 : end) = 1;%#ok % right side    
	
    % sample from the image
    unsampled_indices = find(~sampled);
    points_in_this_tree = min(num_points_per_tree, length(unsampled_indices));

    all_indices(ii, 1:points_in_this_tree) = ...
        randsample(unsampled_indices, points_in_this_tree);

end


end



