function mb_cluster_gp(varargin)

% Clustering method based on the divide and conquer method implemented by
% CJR. This function implements the final gathering of the divided cluster
% results and performs the final clustering

% get a dir listing in the results dir of files matching the naming
% convention of our cluster results

args = u_packargs(varargin,... % the user's input
			 '0', ... % strict mode
			 {'ClusteringFunctionArgs',...
             'PyramidsDir',...
             'FinalFile'}, ... % the mandatory arguments
			 'ClusteringFunction', 'mb_k_means_clustering2',...
             'MaxFinalMemory', 512,...
             'WindowSize', 11,...
             'Level', 2, ...
             'Quiet', 0,...
             'Overlap', 20,...
             'WindowSize2', 0,...
             'Standardise', 1);
clear varargin;

if ~args.Quiet
    display(['mb_cluster_gp function started ', datestr(now)]);
end

%get cluster indices
cluster_idx = get_cluster_idx(...
    args.PyramidsDir, args.Level,...
    args.WindowSize, args.MaxFinalMemory, args.Overlap, args.WindowSize2);

% save the indices as they may be useful for debugging
save([args.FinalFile, '_data_idx'], 'cluster_idx');

if ~args.Quiet
    display('Cluster indices generated');
end

%set up args struct for getting data from indices
data_args.LoadIdx = 0;
data_args.Idx = cluster_idx; clear cluster_idx;
data_args.ImageDir = args.PyramidsDir;
data_args.Level = args.Level;
data_args.Orientation = 1;
data_args.WindowSize = args.WindowSize;
data_args.WindowSize2 = args.WindowSize2;
data_args.Standardise = args.Standardise;

%get data using cluster indices;
[args.ClusteringFunctionArgs.Data stats] =...
    get_band_data_for_clustering(data_args);

if ~args.Quiet
    display('Data sampled');
    display(['About to perform clustering on ',...
        num2str(size(args.ClusteringFunctionArgs.Data,1)), ' points']);
end

% Do the final clustering
final_clustering_result = feval(args.ClusteringFunction, args.ClusteringFunctionArgs); %#ok

%Save stats into results structure
final_clustering_result.Stats = stats;
final_clustering_result.WindowSize = args.WindowSize;
final_clustering_result.WindowSize2 = args.WindowSize2;

% Save the clustering_result
save(args.FinalFile, 'final_clustering_result');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function all_indices = get_cluster_idx(...
    pyr_dir, level, window_size, max_memory, overlap, window_size2)

% Set constant
size_of_data_element = 8; % assume each data element is a double, hence 8 bytes per element.


if pyr_dir(end) ~= filesep
	pyr_dir = [pyr_dir filesep];
end

% first get a directory listing and other info about the image directory
pyr_files = dir([pyr_dir,'*pyramid*']);
number_of_pyrs = length(pyr_files);

% work out how much memory (in bytes) each window will occupy
mem_per_sample = (window_size^2 + window_size2^2) * size_of_data_element;

% work out how many points we can fit in MaxMemory
num_points = floor((max_memory * 1024 * 1024) /  mem_per_sample);

% work out how many points per image this is
num_points_per_image = floor(num_points / number_of_pyrs);

if num_points_per_image < 1
	error('This function assumes that we can sample at least one window from each image')
end

display(['Next Data Function: num_points_per_image = ' num2str(num_points_per_image)])

all_indices = repmat(NaN, number_of_pyrs, num_points_per_image);

for ii = 1:number_of_pyrs
	
    % open the i-th image file
	pyramid = u_load([pyr_dir, pyr_files(ii).name]);
	
	% create a matrix of zeros of the same size as the image
	sampled = zeros(size(pyramid{level, 1}));
    
    % set the edges to ones so that we don't sample windows that lie off
    % the edge of the window
    sampled(1:overlap, :) = 1; % top portion
    sampled(end-overlap+1:end, :) = 1; % bottom portion
    sampled(:, 1:overlap) = 1; % left side
    sampled(:,end-overlap+1 : end) = 1;%#ok % right side    
	
    % sample from the image
    unsampled_indices = find(~sampled);
    points_in_this_image = min(num_points_per_image, length(unsampled_indices));

    all_indices(ii, 1:points_in_this_image) = ...
        randsample(unsampled_indices, points_in_this_image);

end


end