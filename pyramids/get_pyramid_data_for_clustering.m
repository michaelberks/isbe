function data = get_pyramid_data_for_clustering(varargin)

args = u_packargs(varargin,... % the user's input
			 '0', ... % strict mode
			 {'IdxFile',...
             'ImageDir'},...
             'ImageList', [],...
             'Dimensions', 27);
clear varargin;

if args.ImageDir(end) ~= '/'
	args.ImageDir = [args.ImageDir '/'];
end
if isempty(args.ImageList)
    args.ImageList = dir([args.ImageDir,'*.mat']);
end

%get number of images
num_images = length(args.ImageList);
    
%load cluster indices - each rows defines indices for one image
cluster_idx = u_load(args.IdxFile); %load variable cluster_idx

%pre-allocate memory for data
data = repmat(NaN, sum(~isnan(cluster_idx(:))), args.Dimensions); %#ok
curr_data_point = 1;

%looping over each image, sample at the given indices
for jj = 1:num_images

    % get cluster indices and removes any NaNs
    image_idx = cluster_idx(jj, :); %#ok
    image_idx(isnan(image_idx)) = [];
    n_pts = length(image_idx);
    
    %load image
    pyramid = u_load([args.ImageDir, args.ImageList(jj).name]);
    
    [p_rows, p_cols] = ind2sub(size(pyramid{1,1}), image_idx); 
    clear image_idx;
    
    % now sample the detailing_pyramid_coeffs at each of the shosen locations
    data(curr_data_point:curr_data_point+n_pts-1,:) = ...
        mb_get_pyramid_coefficients(pyramid, [p_rows', p_cols']);
    curr_data_point = curr_data_point+n_pts;
end

%remove unfilled rows from data (shouldn't happen now pre-allocating on
%sum(~isnan(cluster_idx)) as opposed to numel(cluster_idx)
data(curr_data_point:end, :) = [];

