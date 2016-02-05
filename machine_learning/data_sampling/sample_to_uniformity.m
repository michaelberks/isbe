function [rand_idx] = sample_to_uniformity(values, num_samples, start_value, end_value)
%SAMPLE_TO_UNIFORMITY *Insert a one line summary here*
%   [rand_idx] = sample_to_uniformity(values, num_samples, start_value, end_value)
%
% Inputs:
%      values - vector of values for existing sample
%
%      num_samples - size of output indices
%
%      start_value - of range to sample uniformly
%
%      end_value - of range to sample uniformly
%
%
% Outputs:
%      rand_idx - output indices, size num_samples x 1, selected so that
%       hist(values(rand_idx)) is approximately uniform
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 02-Dec-2015
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

% Sort values and get associated sorted indexes
[values_sorted sorted_idx] = sort(values);

%Generate a sorted list of random doubles in the desired range
rand_dbls = start_value + (end_value-start_value)*sort(rand(num_samples,1));

%Pre-allocate output
rand_idx = zeros(num_samples,1);

%Now loop and fill the rand samples
curr_idx = 1;
for i_sample = 1:num_samples
    while (values_sorted(curr_idx) < rand_dbls(i_sample))
        curr_idx = curr_idx + 1;
    end
    
    rand_idx(i_sample) = sorted_idx(curr_idx);
end

%Finally sort the rand_idx
rand_idx = sort(rand_idx);
