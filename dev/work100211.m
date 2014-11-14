clear
%profile on;
win_size = 3;
num_levels = 6;
pad_w = floor(win_size/2);
win_idx = -pad_w:pad_w;
do_max = 1;
bg_ratio = 0.5;

if do_max
    size_sample_vector = 2*win_size*win_size*num_levels;
else
    size_sample_vector = 2*win_size*win_size*num_levels*6;
end

%create a test image
bg = u_load('M:\chen\data\normal_smooth512\bg001.mat');
[row col] = size(bg);
[bar_image, dummy, label] = create_ellipse_bar(10, 10, 45, row, col);

%Add bar to background
image_in = bar_image + bg;

%Make an edge mask to avoid edge of image
edge_mask = true(row, col); 
edge_mask([1:pad_w end-pad_w+1], :) = false; 
edge_mask(:, [1:pad_w end-pad_w+1]) = false;

%Use label to get indices for this image
bar_idx = find(label & edge_mask); %use all bar pixels
bg_idx = find(~label & edge_mask);

%Now select a random selection of the background pixels
rand_idx = randperm(length(bg_idx));
bg_idx = bg_idx(rand_idx(1:round(length(bar_idx)*bg_ratio)));

image_idx = [bar_idx; bg_idx];
num_samples_image = length(image_idx);
[rows cols] = ind2sub([row col], image_idx);

%Compute dual-tree of image
dt = dtwavexfm2(image_in, num_levels);

%create storage for data for both methods
training_data1 = zeros(num_samples_image, size_sample_vector);
training_data2 = zeros(num_samples_image, size_sample_vector);

%Now compare the two methods
%1) Using new method only interpolating at sample points
tic;

%Make copies of sample rows and cols at positions of local window patch
rr = repmat(rows*ones(1,win_size) + ones(num_samples_image,1)*win_idx, 1, win_size);
cc = kron(cols*ones(1,win_size) + ones(num_samples_image,1)*win_idx, ones(1,win_size));

%Get interpolated dual-tree coefficients
dt_samples = dt_to_pixel_subset(dt, rr, cc);

%Check whether to calculate maximum
if do_max
    dt_samples = squeeze(max(dt_samples,[],3));
end

%Save the samples in the training data
curr_sample = 1;
training_data1(curr_sample:num_samples_image+curr_sample-1,1:end/2)...
    = reshape(dt_samples, num_samples_image, []);

%Increment the image counter
curr_sample = curr_sample + num_samples_image;
toc;

%2) Using old method interpolating at all points
tic;
dt_full = dt_to_full_image(dt); clear bg_dt;

if do_max
    dt_full = squeeze(max(dt_full, [], 3));
end
pad_dt = padarray(dt_full, [pad_w pad_w], 'replicate');

curr_sample = 1;
%loop through image pixels selecting points
for ii = 1:num_samples_image

    rr = rows(ii) + pad_w;
    cc = cols(ii) + pad_w;

    %Extract label and sample vector
    sample_vector = pad_dt(rr+win_idx, cc+win_idx, :);

    %Convert sample vector to mag/phase form and save in main data
    training_data2(curr_sample,1:end/2) = sample_vector(:).';
    curr_sample = curr_sample + 1;

end

toc;

%Finally check we've sampled the same data
max(abs(training_data1(:) - training_data2(:))) %should be zero!

% %Convert into magnitude and phase
% training_data(:,1+end/2:end) = angle(training_data(:,1:end/2));
% training_data(:,1:end/2) = abs(training_data(:,1:end/2));

%profile viewer;
%%

