function [data stats] = mb_get_dual_tree_subband_data(varargin)

args = u_packargs(varargin,... % the user's input
			 '0', ... % strict mode
			 {'LoadIdx',...
             'Idx',...
             'DualTreeDir',...
             'Level',...
             'Orientation'},...
             'DualTreeList', [],...
             'WindowSize', 11,...
             'WindowSize2', 0,...
             'Real', 1, ....
             'Standardise', 0);
clear varargin;

if args.DualTreeDir(end) ~= '/'
	args.DualTreeDir = [args.DualTreeDir '/'];
end
if isempty(args.DualTreeList)
    args.DualTreeList = dir([args.DualTreeDir,'*dual_tree*']);
end

%get number of images
num_trees = length(args.DualTreeList);

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
for jj = 1:num_trees

    % get cluster indices and removes any NaNs
    image_idx = cluster_idx(jj, :); %#ok
    image_idx(isnan(image_idx)) = [];

    %load image
    dual_tree = u_load([args.DualTreeDir, args.DualTreeList(jj).name]); 
    
    if args.Real
        sample_image = real(dual_tree{args.Level}(:,:, args.Orientation));
    else
        sample_image = imag(dual_tree{args.Level}(:,:, args.Orientation));
    end 
    
    if args.WindowSize2
        if args.Level + 1 == size(dual_tree, 1);
            if args.Real
                sample_image2 = real(dual_tree{args.Level+1, 1});
            else
                sample_image2 = imag(dual_tree{args.Level+1, 1});
            end                
        else
            if args.Real
                sample_image2 = real(dual_tree{args.Level+1}(:,:,args.Orientation));
            else
                sample_image2 = imag(dual_tree{args.Level+1}(:,:,args.Orientation));
            end
        end
    end
    
    clear dual_tree;
    
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
    clear sample_image image_idx rows cols sample
end
%remove unfilled rows from data (shouldn't happen now pre-allocating on
%sum(~isnan(cluster_idx)) as opposed to numel(cluster_idx)
data(curr_data_point+1:end, :) = [];

%Compute mean and SD for data and save such that original_data = ....
% (data*SD) + mean
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
else
    stats.Mean = 0;
    stats.SD = 1;
    stats.Mean2 = 0;
    stats.SD2 = 1;
end
