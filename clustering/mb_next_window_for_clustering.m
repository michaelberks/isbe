function [all_data, finished_sampling_flag] = mb_next_window_for_clustering(varargin)
%
% A function that uniformly samples windows from images that are stored
% in a particular directory for use with the function MB_CLUSTER_LARGE_DATA_SET,
% and indirectly by the texture modelling by Gaussian mixture model functions.
%
% MB_NEXT_WINDOW_FOR_CLUSTERING uses the u_packargs interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% 'ImageDir'
%			- the full path to the directory with the image(s) in. NOTE: This must contain
%			ONLY images, and only filetypes supported by IMREAD.
% 'TempStorageDir'
%			- an empty directory that can be used to record which windows have been
%			sampled already. 
% 'WindowSize'
%			- a scalar specifing the window width (and since the windows
%			are square, the height also) -- this must be odd (so the windows
%			have a centre pixel).
% 'MaxMemory'
%			- the maximum amount of memory to use for the returned data matrix,
%			data (in MBytes)
% 'CleanUp'
%			- if set to true, the function will not return any more data, but will clean up its
%			temporary files. Defaults to false. 
%
% MB_NEXT_WINDOW_FOR_CLUSTERING returns the windows (as vectors) in data
% and sets finished_sampling_flag to (1==1) if there is no more data left, or (1==0) if there is.

% Set constant
size_of_data_element = 8; % assume each data element is a double, hence 8 bytes per element.

% pack the args
args = u_packargs(varargin, ... % what the user sent us
			'0',...
			{...  % The mandatory arguments
            'ImageFiles'},... % The optional arguments
            'ImageDir', [mberksroot, 'background\images\'],...
            'TempStorageDir', [mberksroot, 'background\temp\'],...
			'WindowSize', [], ...
			'MaxMemory', 128, ...
			'CleanUp', (1==0));

% ensure that the paths end in the slash appropriate for this platform
if args.TempStorageDir(end) ~= filesep
	args.TempStorageDir = [args.TempStorageDir filesep];
end
if args.ImageDir(end) ~= filesep
	args.ImageDir = [args.ImageDir filesep];
end

% check to see if we've been asked to clean up after ourselves.
if args.CleanUp
	delete([args.TempStorageDir 'temp_sampled_record_*']); % delete all the files that match the naming convention for the files that store where we've been
	disp('Deleting all the temporary files used by MB_NEXT_WINDOW_FOR_CLUSTERING');
	return; % now just return
end

% first get a directory listing and other info about the image directory
number_of_images = length(args.ImageFiles);

% if we haven't created the temporary files yet, do so
tempdir_list = dir([args.TempStorageDir 'temp_sampled_record_*']); % get a list of the files that match the naming convention for the files that store where we've been
if length(tempdir_list) ~= number_of_images
    disp('mb_next_window_for_clustering: Creating the temp files to store where we''ve sampled');
	create_temp_files(number_of_images, args);
end

% work out how much memory (in bytes) each window will occupy
mem_per_sample = (args.WindowSize * args.WindowSize) * size_of_data_element;

% work out how many points we can fit in args.MaxMemory
num_points = floor((args.MaxMemory * 1024 * 1024) /  mem_per_sample);
if num_points < 1
	error('This function assumes that we can sample at least one window from each image')
end

% work out how many points per image this is
num_points_per_image = floor(num_points / number_of_images);
disp(datestr(now))
disp(['Next Data Function: num_points_per_image = ' num2str(num_points_per_image)])

% now loop through the images sampling the windows
finished_sampling_flag = (1==0);
full_count = 0;
all_data = []; % <SOME PRE-ALLOCATION MAY BE USEFUL HERE>

%h = timebar(0,'Getting samples from images...');
for i = 1:number_of_images
	% load the image
	image = imread([args.ImageDir args.ImageFiles(i).name]);
    disp(datestr(now))
    disp(['Next Data Function: Loaded the image: ' args.ImageFiles(i).name])
	
	% load the corresponding temp file - avariable called 'sampled'
    load([args.TempStorageDir  'temp_sampled_record_' zerostr(i,3)]);
    
	% see if we can sample
	if all(sampled(:))
		full_count = full_count + 1;
		if full_count == number_of_images
			finished_sampling_flag = (1==1);
            disp(datestr(now))
            disp('Next Data Function: We have finished sampling from this image')
			break; % we have finished sampling all the data
		end
	else
		% sample from the image
		unsampled_indices = find(~sampled);
		if length(unsampled_indices) < num_points_per_image
			num_points_per_image = length(unsampled_indices);
		end

		% The next two lines of code achieve random sampling without replacement
		% get random indices into the list of unsampled pixels
		[rows, cols] = ind2sub(size(image), randsample(unsampled_indices, num_points_per_image));
		
		% now actually do the sampling
		% pre-allocation
		data = zeros(length(rows), (args.WindowSize * args.WindowSize));
		for j = 1 : length(rows)
            
			this_sample = sample_window(image, args.WindowSize, rows(j), cols(j));
			data(j,:) = this_sample(:)';
            % now record where we've just sampled
            sampled(rows(j), cols(j)) = 1; %#ok
		end
        all_data = [all_data; data];
        disp(['Next Data Function: There are ' num2str(sum(sampled(:)==0)) ' unsampled points in this image'])
		% save where we've been
		save([args.TempStorageDir  'temp_sampled_record_' zerostr(i, 3)], 'sampled');
        disp(datestr(now))
        disp(['Next Data Function: Sampled ' num2str(size(data,1)) ' points, and have recorded where these came from'])
	end
% 	timebar(i/number_of_images,h);
end
% close(h);

all_data = double(all_data);

% if we have finished sampling, we can delete the temporary files
if finished_sampling_flag
	delete([args.TempStorageDir 'temp_sampled_record_*']); % delete all the files that match the naming convention for the files that store where we've been
	disp('Deleting all the temporary files used by MB_NEXT_WINDOW_FOR_CLUSTERING');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = create_temp_files(number_of_images, args)
%
% Create the temporary files
for i = 1:number_of_images
	% open the i-th image file
	image = imread([args.ImageDir args.ImageFiles(i).name]);
	
	% create a matrix of zeros of the same size as the image
	sampled = zeros(size(image));
    
    % set the edges to ones so that we don't sample windows that lie off the edge of the window
    overlap = (args.WindowSize-1)/2;
    sampled(1:overlap, :) = 1; % top portion
    sampled(end-overlap+1:end, :) = 1; % bottom portion
    sampled(:, 1:overlap) = 1; % left side
    sampled(:,end-overlap+1 : end) = 1; %#ok right side %
	
	% save this as a matlab file in the temp dir
	save([args.TempStorageDir  'temp_sampled_record_' zerostr(i,3)], 'sampled');
end
