function [data] = get_window_data_for_clustering(varargin)

args = u_packargs(varargin,... % the user's input
			 '0', ... % strict mode
			 {'IdxFile',...
             'ImageDir'},...
             'ImageList', [],...
             'Dimensions', 5);
clear varargin;

if args.ImageDir(end) ~= filesep
	args.ImageDir = [args.ImageDir filesep];
end
if isempty(args.ImageList)
    args.ImageList = dir([args.ImageDir,'*.bmp']);
end

%get number of images
num_images = length(args.ImageList);
    
%load cluster indices - each rows defines indices for one image
cluster_idx = u_load(args.IdxFile); %load variable cluster_idx

%pre-allocate memory for data
data = repmat(NaN, sum(~isnan(cluster_idx(:))), args.Dimensions.^2); %#ok
curr_data_point = 0;

%looping over each image, sample at the given indices
for jj = 1:num_images

    % get cluster indices and removes any NaNs
    image_idx = cluster_idx(jj, :); %#ok
    image_idx(isnan(image_idx)) = [];

    %load image
    image = imread([args.ImageDir, args.ImageList(jj).name]); 
    %convert indices to subscripts
    [rows, cols] = ind2sub(size(image), image_idx);
    
    for kk = 1:length(image_idx)

        %update data point counter
        curr_data_point = curr_data_point + 1;

        %sample from image
        sample = sample_window(image, args.Dimensions, rows(kk), cols(kk));
        data(curr_data_point,:) = sample(:)';

    end
    clear image image_idx rows cols sample
end
%remove unfilled rows from data (shouldn't happen now pre-allocating on
%sum(~isnan(cluster_idx)) as opposed to numel(cluster_idx)
data(curr_data_point+1:end, :) = [];