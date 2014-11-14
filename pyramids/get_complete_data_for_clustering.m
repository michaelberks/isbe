function data = get_complete_data_for_clustering(varargin)

args = u_packargs(varargin,... % the user's input
			 '0', ... % strict mode
			 {'LoadIdx',...
             'Idx',...
             'ImageDir',...
             'Level'},...
             'ImageList', [],...
             'WindowSize', 15,...
             'NumLevels', 5,...
             'NumOrientations', 5,...
             'MaxMemory', 128);
clear varargin;

if args.ImageDir(end) ~= '/'
	args.ImageDir = [args.ImageDir '/'];
end
if isempty(args.ImageList)
    args.ImageList = dir([args.ImageDir,'*.mat']);
end

%Calculate the overlap needed to ensure sample window won't go out of
%bounds
overlap = (args.WindowSize + 1) / 2;

%calculate dimension of complete data
pyr_dimension = (args.NumLevels - args.Level + 1)*args.NumOrientations + 1;

if args.Level == 1 || args.Level == args.NumLevels
    dimension = pyr_dimension + args.WindowSize^2;
else
	dimension = pyr_dimension + (args.WindowSize^2) * args.NumOrientations;
end

%get number of images
num_images = length(args.ImageList);

% work out how much memory (in MB) each window will occupy
mem_per_sample = dimension * 2^-17;

% work out how many points we can fit in args.MaxMemory
num_points = floor(args.MaxMemory /  mem_per_sample);

if num_points < 1
	error('This function assumes that we can sample at least one window from each image')
end

% work out how many points per image this is
num_points_per_image = floor(num_points / num_images);

if args.LoadIdx
    %load cluster indices - each rows defines indices for one image
    cluster_idx = u_load(args.Idx); %load variable cluster_idx
else
    % Indices are already stored in args.Idx
    cluster_idx = args.Idx;
    args = rmfield(args, 'Idx'); %memory may be short so clear Idx from args
end

%pre-allocate memory for data
data = repmat(NaN, num_points_per_image*num_images, dimension);

start_point = 1;

%looping over each image, sample at the given indices
for jj = 1:num_images
    %load pyramid
    pyramid = u_load([args.ImageDir, args.ImageList(jj).name]);
    
    % get cluster indices and removes any NaNs
    image_idx = cluster_idx(jj, :); %#ok
    image_idx(isnan(image_idx)) = [];
    
    %convert indices to rows/cols in finest pyramid level
    [rows, cols] = ind2sub(size(pyramid{1, 1}), image_idx');
    
    %only need to downsample for levels 3 and above
    if args.Level > 2
        %downsample rows/cols to args.Level
        rows = ceil(rows/2^(args.Level-2));
        cols = ceil(cols/2^(args.Level-2));
        
        %having down-sampled rows/cols move closer to edge and can force
        %sample window out of bounds
        scrap = rows < overlap | cols < overlap;        
        rows(scrap) = [];
        cols(scrap) = [];
        
        %make rows and cols unique
        rc = unique([rows cols], 'rows');
        rows = rc(:,1); cols = rc(:,2); clear rc;
    end
    sub_length = length(rows);
    
    %sub-sample the rows/cols if we have to
    if  sub_length > num_points_per_image
        rows = rows(randsample(sub_length, num_points_per_image));
        cols = cols(randsample(sub_length, num_points_per_image));
        num_points_this_image = num_points_per_image;
    else
        num_points_this_image = sub_length;
    end
    
    %get pyramid coeffs from all levels coarser than args.Level
    if args.Level > 1
        pyr_coeffs = mb_get_pyramid_coefficients(pyramid, 2^(args.Level-2)*[rows, cols]);
    else
        pyr_coeffs = mb_get_pyramid_coefficients(pyramid, [rows, cols]);
    end
    
    data(start_point:start_point+num_points_this_image-1, 1:pyr_dimension) =...
        pyr_coeffs(:, end-pyr_dimension+1:end);
    
    %get co-effs from each orientation of args.Level
    for ori = 1:args.NumOrientations
        
        image = pyramid{args.Level, ori};

        ori_start_dim = pyr_dimension + (args.WindowSize^2) * (ori - 1) + 1; 
        ori_end_dim = ori_start_dim + (args.WindowSize^2) - 1;
        
        %start count again for each orientation
        curr_data_point = start_point;
        for kk = 1:num_points_this_image

            %sample from image
            sample = sample_window(image, args.WindowSize, rows(kk), cols(kk));
            data(curr_data_point, ori_start_dim:ori_end_dim) = sample(:)';
            
            %update data point counter
            curr_data_point = curr_data_point + 1;
        end
        if args.Level == 1 || args.Level == args.NumLevels
            %Only one sub-band in first and last levels
            break
        end
    end
    start_point = start_point + num_points_this_image;
end
%remove unfilled rows from data
data(curr_data_point+1:end, :) = [];

