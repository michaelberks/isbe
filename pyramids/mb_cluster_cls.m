function mb_cluster_cls(varargin)

% Clustering method based on the divide and conquer method implemented by
% CJR. This function implements the final gathering of the divided cluster
% results and performs the final clustering

% get a dir listing in the results dir of files matching the naming
% convention of our cluster results

args = u_packargs(varargin,... % the user's input
			 '0', ... % strict mode
			 {'ClusteringFunctionArgs',...
             'PyramidsDir',...
             'CLSDir',...
             'FinalFile',...
             'CLSMode'}, ... % the mandatory arguments
			 'ClusteringFunction', 'mb_k_means_clustering2',...
             'MaxFinalMemory', 512,...
             'WindowSize', 11,...
             'Level', 2, ...
             'Orientation', 1,...
             'Quiet', 0,...
             'Overlap', 20,...
             'WindowSize2', 0,...
             'Standardise', 0);
clear varargin;

if ~args.Quiet
    display(['mb_cluster_cls function started ', datestr(now)]);
end

%get cluster indices
cluster_idx = get_cluster_idx(...
    args.PyramidsDir, args.CLSDir, args.Level, args.Orientation, args.WindowSize,...
    args.MaxFinalMemory, args.Overlap, args.WindowSize2, args.CLSMode);

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
data_args.Orientation = args.Orientation;
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
    pyr_dir, cls_dir, level, orientation, window_size, max_memory, overlap, window_size2, cls_mode)

% Set constant
size_of_data_element = 8; % assume each data element is a double, hence 8 bytes per element.


if pyr_dir(end) ~= filesep
	pyr_dir = [pyr_dir filesep];
end
if cls_dir(end) ~= filesep
	cls_dir = [cls_dir filesep];
end

% first get a directory listing and other info about the pyramid directory
pyr_files = dir([pyr_dir,'*pyramid*']);
number_of_pyrs = length(pyr_files);

% first get a directory listing and other info about the pyramid directory
cls_files = dir([cls_dir,'*cls*']);

if length(cls_files) ~= number_of_pyrs
    %error('There must be matching CLS files for every Gaussian pyramid');
end

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
	cls = u_load([cls_dir, cls_files(ii).name]);
    
    cls_map = mb_cls_map_orientations(cls);
	% create a matrix of zeros of the same size as the image
	sampled = cls_map{level-1} == orientation;
    
    if any(size(sampled) ~= size(pyramid{level,1}))
        error('CLS map doesn''t match pyramid level')
    end
    
    if cls_mode
        %If we're selecting CLS pixels we need these to be zeros in the
        %image map so take inverse of cls map
        sampled = ~sampled;
    end
    
    % set the edges to ones so that we don't sample windows that lie off
    % the edge of the window
    sampled(1:overlap, :) = 1; % top portion
    sampled(end-overlap+1:end, :) = 1; % bottom portion
    sampled(:, 1:overlap) = 1; % left side
    sampled(:,end-overlap+1 : end) = 1;%#ok % right side    
	
    % sample from the image
    unsampled_indices = find(~sampled);
    points_in_this_image = min(num_points_per_image, length(unsampled_indices));
    
    if points_in_this_image > 0
        all_indices(ii, 1:points_in_this_image) = ...
            randsample(unsampled_indices, points_in_this_image);
    end

end


end



