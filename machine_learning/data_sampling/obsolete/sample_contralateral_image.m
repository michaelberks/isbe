function [training_data training_labels test_data test_labels row col image_idx_abnormal image_idx_normal] = sample_contralateral_image(varargin)
%
% SAMPLE_TRAINING_DATA_ALLSUBBANDS create training data which use all
% subbands coefficients of DTCWT
% USAGE:
% [training_data training_labels parameters] =...
%   sample_training_data_allsubbands('num_samples', 1000, 'bg_dir',...
%   'E:\DTCWTmini\data\normal_smooth512\', 'save_path',
%   'E:\DTCWTmini\data\synthetic_512\', 'plot', 1);
%
% Inputs:
%      num_samples - the total number of samples
%      bg_dir -  the patch of background patch
%      halfwidth - halfwidth of Gaussian profile at half its maximum height
%
%      contrast - maximum height of (scaled) Gaussian profile
%
%      orientation - orientation of bar in image in degrees
%
%      row - number of rows in image
%
%      col - number of columns in image
%
%
% Outputs:
%      image_out - image containing spicule
%
%      label - label of spicule (1) vs background (0)
%
%      label_centre - the label of centre line (1) vs background (0)
%
% Example:
%
% Notes:
% See also: sample_bar_training_data
%
% Created: 08-February-2010
% Author: Zezhi Chen and Michael Berks
% Email : zezhi.chen@manchester.ac.uk
% Phone : +44 (0)161 275 7669
% Copyright: (C) University of Manchester
%


% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    {'pair_name'},... % the mandatory arguments
    'win_size', 3,...
    'num_levels', 6,...
    'do_max', 0,...
    'num_train', 1e3,...
    'num_test', 1e3,...
    'feature_type', 'all',...
    'save_path', [], ...
    'plot', 0);

win_size = args.win_size; % the size of scanning window
num_levels = args.num_levels; %number of levels of DT-CWT from which to extract coefficients
num_train = args.num_train; %number of feature vectors to sample for training from each region
num_test = args.num_test; %number of feature vectors to sample for training from each region

%work out size of feature vectors
if strcmp(args.feature_type, 'all')
    feature_length = 2;
else
    feature_length = 1;
end

if args.do_max
    size_sample_vector = feature_length*win_size*win_size*num_levels;
else
    size_sample_vector = feature_length*win_size*win_size*num_levels*6;
end

contralateral_pair = u_load(args.pair_name);
num_pixels = numel(contralateral_pair.abnormal_roi);

if num_train > num_pixels
    display(['Warning: Number of training samples (' num2str(num_train) ') is greater than' ...
        ' the number of pixels in the image (' num2str(num_train) '). No test data will be selected']);
    
    num_train = num_pixels;
end
if num_train + num_test > num_pixels
    num_test = num_pixels - num_train;
end

%Create memory storage for data
training_data = zeros(2*num_train, size_sample_vector);
training_labels = [true(num_train, 1); false(num_train, 1)];

test_data = zeros(2*num_test, size_sample_vector);
test_labels = [true(num_test, 1); false(num_test, 1)];

%Extract training and test data from abnormal patch
[training_data(1:num_train,:) test_data(1:num_test,:) row col image_idx_abnormal] = ...
    sample_data(contralateral_pair.abnormal_roi, num_train, num_test, win_size, num_levels, args.feature_type, args.do_max);

%Extract training and test data from normal patch
[training_data(num_train+1:end,:) test_data(num_test+1:end,:) row col image_idx_normal] = ...
    sample_data(contralateral_pair.normal_roi, num_train, num_test, win_size, num_levels, args.feature_type, args.do_max);

%--------------------------------------------------------------------------
function [training_data test_data row col image_idx] =...
        sample_data(image_in, num_train, num_test, win_size, num_levels, feature_type, do_max)

%Constants computed from arguments
pad_w = floor(win_size/2); %half size of window size
win_idx = -pad_w:pad_w;

%get training and testing idx from abnormal image
[row col] = size(image_in);

image_idx = randperm(row*col)';
image_idx(num_train+num_test+1:end) = [];

[rows cols] = ind2sub([row col], image_idx);
    
% Create DT-CWT of image
dt = dtwavexfm2b(image_in, num_levels); clear image_in;

%Make copies of sample rows and cols at positions of local window patch
rr = repmat(rows*ones(1,win_size) + ones(num_train+num_test,1)*win_idx, 1, win_size);
cc = kron(cols*ones(1,win_size) + ones(num_train+num_test,1)*win_idx, ones(1,win_size));
    
%Get interpolated dual-tree coefficients
dt_samples = dt_to_pixel_subset(dt, rr, cc); clear dt rr cc;

if do_max
    %get the maximum response across orientations
    dt_samples = squeeze(max(dt_samples, [], 3)); clear dt;
end
    
temp_samples=reshape(dt_samples, num_train+num_test, []);
    
switch feature_type
    case 'all'
        training_data = [abs(temp_samples(1:num_train,:)) angle(temp_samples(1:num_train,:))];
        test_data = [abs(temp_samples(num_train+1:end,:)) angle(temp_samples(num_train+1:end,:))];
        
    case 'real'
        training_data = real(temp_samples(1:num_train,:));
        test_data = real(temp_samples(num_train+1:end,:));

    case 'mag'
        training_data = abs(temp_samples(1:num_train,:));
        test_data = abs(temp_samples(num_train+1:end,:));

    case 'phase'
        training_data = angle(temp_samples(1:num_train,:));
        test_data = angle(temp_samples(num_train+1:end,:));

    otherwise
        warning(['Feature type: ', args.feature_type, ' not recognised']); %#ok
        training_data = [abs(temp_samples(1:num_train,:)) angle(temp_samples(1:num_train,:))];
        test_data = [abs(temp_samples(num_train+1:end,:)) angle(temp_samples(num_train+1:end,:))];
end
