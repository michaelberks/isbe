function [data stats] = mb_get_dual_tree_data_all_subbands(varargin)

args = u_packargs(varargin,... % the user's input
			 '0', ... % strict mode
			 {'LoadIdx',...
             'Idx',...
             'DualTreeDir',...
             'Level'},...
             'DualTreeList', [],...
             'WindowSize', 5,...
             'WindowSize2', 3,...
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

%Work out dimensions for data ((re + im) * (number of points
%in this level window + number of points in next level window) in a single
%orientation band
dimensions = 2*(args.WindowSize^2 + args.WindowSize2^2);

%pre-allocate memory for data
data = repmat(NaN, sum(~isnan(cluster_idx(:))), 6*dimensions); %#ok

% Set up the expected phase shifts for interpolating each subband:
w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.15;

tree_data_point = 0;
%looping over each image, sample at the given indices
for jj = 1:num_trees
    
    %load image
    dual_tree = u_load([args.DualTreeDir, args.DualTreeList(jj).name]);
    
    %get dimensions of sub-bands in this level
    [R C] = size(dual_tree{args.Level}(:,:,1));
    
    % get cluster indices and removes any NaNs
    image_idx = cluster_idx(jj, :); %#ok
    image_idx(isnan(image_idx)) = [];
    
    %convert indices to subscripts
    [rows, cols] = ind2sub([R C], image_idx);
  
    %for each oriented sub-band
    for band = 1:6,
            
        %Get real and imaginary parts of sub-band
        sample_image_re1 = real(dual_tree{args.Level}(:,:,band));
        sample_image_im1 = imag(dual_tree{args.Level}(:,:,band)); 
    
        %Check if we're sampling from the next coarsest sub-band
        if args.WindowSize2
            %if so, interpolate coarser level
            temp = cpxinterp2(dual_tree{args.Level+1}(:,:,band), [-.25 .25], w(band,:),'spline');        
            dt_int = temp(1:R, 1:C);
            clear temp;
            
            %get real and imaginary parts of each sub-band
            sample_image_re2 = real(dt_int);
            sample_image_im2 = imag(dt_int);
        end        

        %Now work through each index point and sample the data
        curr_data_point = tree_data_point;
        for kk = 1:length(image_idx)

            %update data point counter
            curr_data_point = curr_data_point + 1;

            %sample from image at current level
            sample_win_re1 = sample_window(sample_image_re1, args.WindowSize, rows(kk), cols(kk));
            sample_win_im1 = sample_window(sample_image_im1, args.WindowSize, rows(kk), cols(kk));
            sample_win_re2 = [];
            sample_win_im2 = [];

            %sample from lower level if required
            if args.WindowSize2
                sample_win_re2 = ...
                    sample_window(sample_image_re2, args.WindowSize2, rows(kk), cols(kk));
                sample_win_im2 = ...
                    sample_window(sample_image_im2, args.WindowSize2, rows(kk), cols(kk));
            end
            
            %save samples in data structure
            data(curr_data_point,(band-1)*dimensions + 1:band*dimensions) =...
                [sample_win_re1(:)' sample_win_im1(:)' sample_win_re2(:)' sample_win_im2(:)'];

        end   
    end
    tree_data_point = curr_data_point;    
    clear dual_tree;
end
%remove unfilled rows from data (shouldn't happen now pre-allocating on
%sum(~isnan(cluster_idx)) as opposed to numel(cluster_idx)
data(curr_data_point+1:end, :) = [];

%Compute mean and SD for data and save such that original_data = ....
% (data*SD) + mean
if args.Standardise
    stats.Mean = mean(data);
    stats.SD = std(data);
    data = (data - repmat(stats.Mean, size(data,1), 1)) ./ ...
        repmat(stats.SD, size(data,1), 1);
else
    stats.Mean = 0;
    stats.SD = 1;
end
