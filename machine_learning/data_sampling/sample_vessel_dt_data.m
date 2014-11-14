function [training_data training_labels] = sample_vessel_dt_data(varargin)
%SAMPLE_VESSEL_DT_DATA *Insert a one line summary here*
%   [] = sample_mammo_training_data(varargin)
%
% SAMPLE_VESSEL_DT_DATA uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 17-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    {'num_samples',... % the mandatory arguments
    'saved_data_dir'}, ...
    'rgb_channel', 'all',...
    'win_size', 3,...
    'num_levels', 4,...
    'feature_type', 'conj',...
    'feature_shape', 'rect',...
    'rotate', 0,...
    'do_max', 0,...
    'data_list', [], ...
    'label_list', [], ...
    'data_name', 'X',...
    'label_name', 'y',...
    'save_path', [],...
	'pca', []);

clear varargin;

num_samples = args.num_samples; %total number of feature vectors to sample
data_list = args.data_list;
label_list = args.data_list;


if length(args.num_levels) == 1;
    levels = 1:args.num_levels;
else
    levels = args.num_levels;
end

%Convert DT representation and allocate to main data
rep_args = get_substructure(args,...
    {'feature_type', 'do_max', 'win_size', 'pca'});
    
sampling_args.feature_shape = args.feature_shape;
sampling_args.feature_type = args.feature_type;
sampling_args.do_max = args.do_max;
sampling_args.rotate = args.rotate;
sampling_args.win_size = args.win_size;

switch args.rgb_channel %could use if but keep switch in case of further options
    case 'all'
        num_channels = 3;
        
    case 'rgb'
        num_channels = 1;
        
    otherwise %could be 'r', 'g', 'b'
        num_channels = 1;
end
samples_per_channel = args.num_levels*compute_dt_feature_size(sampling_args);
size_sample_vector = num_channels*samples_per_channel;

if isempty(data_list)
    if num_channels == 3
        data_list = [...
            dir([args.saved_data_dir, '/' args.data_name '*_r_*']) ...
            dir([args.saved_data_dir, '/' args.data_name '*_g_*']) ...
            dir([args.saved_data_dir, '/' args.data_name '*_b_*'])];
    else
        data_list = dir([args.saved_data_dir, '/' args.data_name '*_' args.rgb_channel '_*']);
    end
end
if isempty(label_list)
    label_list = dir([args.saved_data_dir, '/' args.label_name '*']);
end

%Get list of images in image dir
%Create data
training_data = zeros(num_samples, size_sample_vector);
training_labels = zeros(num_samples, 1);

curr_sample = 1;
num_data = length(data_list);
data_order = randperm(num_data);

%loop through each image sampling data
for kk = 1:num_data
    
    %num_samples_image = binornd((num_samples + 1 - curr_sample), 1/(num_images+1-kk), 1, 1);
    num_samples_data = ...
        sample_from_binomial((num_samples + 1 - curr_sample), 1/(num_data+1-kk), 1);
    
    %display(['Samples in image ' num2str(kk) ': ', num2str(num_samples_image)]);
    
    if ~num_samples_data
        continue;
    end
    
    %load in orientations
    y = u_load([args.saved_data_dir label_list(data_order(kk)).name]);
    
    %Check we have enough samples in data
    num_samples_data = min(num_samples_data, size(y,1)); %just in case...
    
    %Select random samples to keep
    rand_idx = randperm(size(y,1));
    rand_idx = rand_idx(1:num_samples_data);
    
    %Copy training labels into main structure
    training_labels(curr_sample:num_samples_data+curr_sample-1,:) = ...
        y(rand_idx,:);
    
    %load in training data from X files
    for jj = 1:num_channels
        X = u_load([args.saved_data_dir data_list(data_order(kk),jj).name]);
        
        %Discard unwanted window dimensions if using 1x1 windows
        if args.win_size == 1
            X = X(:,5,:,:);
        end
        
        %Discard unwanted levels
        X = X(:,:,:,levels);
        
        %Workout columns to store data
        cols = (1:samples_per_channel) + (jj-1)*samples_per_channel;
        
        %Discard unwanted rows, convert DT form and save in main output structure
        training_data(curr_sample:num_samples_data+curr_sample-1,cols) = ...
            convert_dt_representation(X(rand_idx,:,:,:), rep_args);
        
    end
   
    %Update the current sample count
    curr_sample = curr_sample + num_samples_data;
end