function [sample_vector_size] = ...
    get_sample_vector_size(decomposition_args, sampling_args)

f_debug = (nargin == 0 && nargout == 0);
if f_debug, test_script(); return; end

[sample_vector_size] = func(decomposition_args, sampling_args);


%% The function
function [sample_vector_size] = func(decomposition_args, sampling_args)

samples_per_pixel = get_samples_per_pixel(decomposition_args);
n_pixels = sampling_args.win_size^2;
[n_channels, channels] = get_channels(sampling_args);
sample_vector_size = samples_per_pixel * n_pixels * n_channels;


%% Test script
function test_script()
clc;
