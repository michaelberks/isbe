function [responses orientations] = mb_get_gabor_data(image_dir, max_memory)

% Set constant
size_of_data_element = 8; % assume each data element is a double, hence 8 bytes per element.

if image_dir(end) ~= filesep
	image_dir = [image_dir filesep];
end

% first get a directory listing and other info about the image directory
image_files = dir([image_dir,'*.bmp*']);
number_of_images = length(image_files);

% work out how much memory (in bytes) each window will occupy
mem_per_sample = size_of_data_element;

% work out how many points we can fit in max_memory
num_pts = floor((max_memory * 1024 * 1024) /  mem_per_sample);

% work out how many points per image this is
num_pts_per_image = floor(num_pts / number_of_images);

if num_pts_per_image < 1
	error('This function assumes that we can sample at least one window from each image')
end

display(['Next Data Function: num_points_per_image = ' num2str(num_pts_per_image)])

%pre-allocate for maximum amount of data
responses = repmat(NaN, num_pts, 1);
orientations = repmat(NaN, num_pts, 1);

curr_idx = 1;

for ii = 1:number_of_images
	
    % open the i-th image file
	image_in = imread([image_dir, image_files(ii).name]);
    [response orientation] = gabor_filter('ImageIn', image_in);
    
	[rows cols] = size(image_in);
    curr_pts = min(num_pts_per_image, rows*cols);
    
    idx = randsample(1:rows*cols, curr_pts);
    
    responses(curr_idx:curr_idx + curr_pts - 1, :) = response(idx);
    orientations(curr_idx:curr_idx + curr_pts - 1, :) = orientation(idx);
    
    curr_idx = curr_idx + curr_pts;    
end

%clear up any space in data we haven't assigned (e.g. if less points in
%dual-trees than maximum size of data specified)
responses(curr_idx:end,:) = [];
orientations(curr_idx:end,:) = [];