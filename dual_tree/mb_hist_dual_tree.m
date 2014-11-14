function [data stds] = mb_hist_dual_tree(dual_tree_dir, max_memory, levels)

% Comment!!

%get cluster indices
data = ...
    get_cluster_data(dual_tree_dir, max_memory, levels);

stds = std(data);
%means = mean(data);

for dim = 1:2:length(stds)-1;
    data(:,dim) =...
        data(:,dim) / stds(dim);
        %(data(:,dim) - means(dim)) ./ stds(dim);
end
% figure;
% %Work through the dimensions of data
% for d = 1:6
%     subplot(2,3,d); hold on;
%     %figure; hold on;
%     
%     [bin_counts bin_centres] = hist3(data(:,[d d+1]), [100 100]);
%     imagesc(bin_counts); axis image;
%     set(gca, 'xtick', bin_centres{1}(10:10:end), 'xticklabel', bin_centres{1}(10:10:end));
%     set(gca, 'xtick', bin_centres{2}(10:10:end), 'xticklabel', bin_centres{2}(10:10:end));
%     title(['Dims ', num2str(d), ' and ', num2str(d+1)]);  
% end

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



