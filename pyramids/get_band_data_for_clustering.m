function [data stats] = get_band_data_for_clustering(varargin)

args = u_packargs(varargin,... % the user's input
			 '0', ... % strict mode
			 {'LoadIdx',...
             'Idx',...
             'ImageDir',...
             'Level',...
             'Orientation'},...
             'ImageList', [],...
             'WindowSize', 11,...
             'WindowSize2', 0,...
             'Standardise', 0);
clear varargin;

if args.ImageDir(end) ~= '/'
	args.ImageDir = [args.ImageDir '/'];
end
if isempty(args.ImageList)
    args.ImageList = dir([args.ImageDir,'*.mat']);
end

%get number of images
num_images = length(args.ImageList);

if args.LoadIdx
    %load cluster indices - each rows defines indices for one image
    cluster_idx = u_load(args.Idx); %load variable cluster_idx
else
    % Indices are already stored in args.Idx
    cluster_idx = args.Idx;
    args = rmfield(args, 'Idx'); %memory may be short so clear Idx from args
end

dimensions = args.WindowSize^2 + args.WindowSize2^2;

%pre-allocate memory for data
data = repmat(NaN, sum(~isnan(cluster_idx(:))), dimensions); %#ok
curr_data_point = 0;

%looping over each image, sample at the given indices
for jj = 1:num_images

    % get cluster indices and removes any NaNs
    image_idx = cluster_idx(jj, :); %#ok
    image_idx(isnan(image_idx)) = [];

    %load image
    pyramid = u_load([args.ImageDir, args.ImageList(jj).name]); 
    sample_image = pyramid{args.Level, args.Orientation};
    
    if args.WindowSize2
        if args.Level + 1 == size(pyramid, 1);
            sample_image2 = pyramid{args.Level+1, 1};
        else
            sample_image2 = pyramid{args.Level+1, args.Orientation};
        end
    end
    
    clear pyramid;
    
    %convert indices to subscripts
    [rows, cols] = ind2sub(size(sample_image), image_idx);
    
    for kk = 1:length(image_idx)

        %update data point counter
        curr_data_point = curr_data_point + 1;

        %sample from image at current level
        sample = sample_window(sample_image, args.WindowSize, rows(kk), cols(kk));
        sample2 = [];
        %sample from lower level if required
        if args.WindowSize2
            sample2 = ...
                sample_window(sample_image2, args.WindowSize2,...
                ceil(rows(kk)/2), ceil(cols(kk)/2));
        end
            
        data(curr_data_point,:) = [sample(:)' sample2(:)'];

    end
    clear image image_idx rows cols sample
end
%remove unfilled rows from data (shouldn't happen now pre-allocating on
%sum(~isnan(cluster_idx)) as opposed to numel(cluster_idx)
data(curr_data_point+1:end, :) = [];
stats = [];

if args.Standardise
    stats.Mean = mean(reshape(data(:, 1:args.WindowSize^2), 1, []));
    stats.SD = std(reshape(data(:, 1:args.WindowSize^2), 1, []));
    data(:,1:args.WindowSize^2) = ...
        (data(:,1:args.WindowSize^2) - stats.Mean) / stats.SD;
    
    if args.WindowSize2
        stats.Mean2 = mean(reshape(data(:, args.WindowSize^2+1:end), 1, []));
        stats.SD2 = std(reshape(data(:, args.WindowSize^2+1:end), 1, []));
        data(:, args.WindowSize^2+1:end) = ...
            (data(:, args.WindowSize^2+1:end) - stats.Mean2) / stats.SD2;
    end
end
