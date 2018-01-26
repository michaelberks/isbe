function [bin_counts, bin_centres] = compute_weighted_histogram(input_data, input_weights, bin_centres)
%COMPUTE_WEIGHTED_HISTOGRAM compute a weighted histogram using the sparse
%matrix trick
%   [bin_counts, bin_centres] = compute_weighted_histogram(input_data, bin_centres)
%
% Inputs:
%      input_data - Used to determine each bin
%      input_weight - Weights used to increment the bin counts. Must be
%      same size as input_data
%
%      bin_centres - If a single integer value N is given, will select N
%      evenly spaced bins over the range of input_data. Otherwise must be
%      evenly spaced and monotonically increasing
%
%
% Outputs:
%      bin_counts - *Insert description of input variable here*
%
%      bin_centres - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 18-Sep-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
input_data = input_data(:);
input_weights = input_weights(:);

%Set up the bins
if numel(bin_centres) == 1
    bin_centres = linspace(min(input_data), max(input_data), bin_centres);
end
num_bins = length(bin_centres);

%Compute the bin index for each input data
bin_spacing = bin_centres(2) - bin_centres(1);
bin_idx = round((input_data - bin_centres(1)) / bin_spacing) + 1;
bin_idx(bin_idx < 1) = 1;
bin_idx(bin_idx > num_bins) = num_bins;

%Compute weighted hist using a sparse matrix
bin_counts = sparse(bin_idx, 1, input_weights, num_bins, 1);
