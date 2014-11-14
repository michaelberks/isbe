function [data ilp_data icp_data] = mb_get_dual_tree_data(dual_tree_dir, level, max_memory, get_ilp, get_icp)

%Set default arguments
if nargin < 4
    get_ilp = 0;
    get_icp = 0;
end

% Set constant
size_of_data_element = 8; % assume each data element is a double, hence 8 bytes per element.

if dual_tree_dir(end) ~= filesep
	dual_tree_dir = [dual_tree_dir filesep];
end

% first get a directory listing and other info about the image directory
DualTreeFiles = dir([dual_tree_dir,'*dual_tree*']);
number_of_trees = length(DualTreeFiles);

% work out how much memory (in bytes) each window will occupy
mem_per_sample = 6 * size_of_data_element;

% work out how many points we can fit in max_memory
num_pts = floor((max_memory * 1024 * 1024) /  mem_per_sample);

% work out how many points per image this is
num_pts_per_tree = floor(num_pts / number_of_trees);

if num_pts_per_tree < 1
	error('This function assumes that we can sample at least one window from each image')
end

display(['Next Data Function: num_points_per_image = ' num2str(num_pts_per_tree)])

%pre-allocate for maximum amount of data
data = repmat(NaN, num_pts, 6);

ilp_data = [];
if get_ilp
    ilp_data = repmat(NaN, num_pts, 6);
end
icp_data = [];
if get_ilp
    icp_data = repmat(NaN, num_pts, 6);
end

curr_idx = 1;

for ii = 1:number_of_trees
	
    % open the i-th image file
	dual_tree = u_load([dual_tree_dir, DualTreeFiles(ii).name]);
    if get_ilp || get_icp
        [ilp icp] = mb_dual_tree_transform(dual_tree);
    end
    
	[rows cols] = size(dual_tree{level}(:,:,1));
    curr_pts = min(num_pts_per_tree, rows*cols);
	
    dt_lev = reshape(dual_tree{level}, [], 6);
    
    idx = randsample(1:rows*cols, curr_pts);
    
    data(curr_idx:curr_idx + curr_pts - 1, :) = dt_lev(idx, :);
    if get_ilp
        ilp_lev = reshape(ilp{level}, [], 6);
        ilp_data(curr_idx:curr_idx + curr_pts - 1, :) = ilp_lev(idx, :);
    end
    if get_icp
        icp_lev = reshape(icp{level}, [], 6);
        icp_data(curr_idx:curr_idx + curr_pts - 1, :) = icp_lev(idx, :);
    end
    curr_idx = curr_idx + curr_pts;    
end

%clear up any space in data we haven't assigned (e.g. if less points in
%dual-trees than maximum size of data specified)
data(curr_idx:end,:) = [];
if get_ilp
    ilp_data(curr_idx:end,:) = [];
end
if get_ilp
    icp_data(curr_idx:end,:) = [];
end