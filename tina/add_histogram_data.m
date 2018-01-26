function [hist_struc] = add_histogram_data(hist_struc, data, weights)
%ADD_HISTOGRAM_DATA *Insert a one line summary here*
%   [hist_struc] = add_histogram_data(hist_struc, data, weights)
%
% Inputs:
%      hist_struc - *Insert description of input variable here*
%
%      data - *Insert description of input variable here*
%
%      weights - *Insert description of input variable here*
%
%
% Outputs:
%      hist_struc - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 29-Jun-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

%Make sure data is col vector and stored as double
if size(data,2) > 1
    data = data(:);
end
if ~isa(data, 'double')
    data = double(data);
end
n = numel(data);

%Make sure we have valid weights
if exist('weights', 'var') 
    if numel(weights) ~= n
        error('If weights are supplied they must match the size of the input data');
    elseif size(weights,2) > 1
        weights = weights(:);
    end
else
    weights = ones(n,1);
end
    
%First compute data under/over or in range
idx_under = data < hist_struc.xmin;
idx_over = data > hist_struc.xmax;
idx_valid = ~idx_under & ~idx_over;

%Update the over under totals (but don't include these in hist summary
%stats)
hist_struc.under = hist_struc.under + sum(weights(idx_under));
hist_struc.over = hist_struc.over + sum(weights(idx_over));

%First compute nearest bin index for each data point
lb_idx = floor((data(idx_valid) - hist_struc.xmin)/ hist_struc.xincr) + 1;
ub_idx = lb_idx + 1;

ub_weights = (data(idx_valid) - hist_struc.xmin - (lb_idx-1)*hist_struc.xincr) / hist_struc.xincr;
lb_weights = 1 - ub_weights;

lb_weights = lb_weights .* weights(idx_valid);
ub_weights = ub_weights .* weights(idx_valid);

idx_nonzero_l = abs(lb_weights) > 0 & lb_idx >= 1;
idx_nonzero_u = abs(ub_weights) > 0 & ub_idx <= hist_struc.nbins;

%Finally use sparse trick to compute histogram
counts = full(sparse(lb_idx(idx_nonzero_l), 1, lb_weights(idx_nonzero_l), hist_struc.nbins, 1)) +...
    full(sparse(ub_idx(idx_nonzero_u), 1, ub_weights(idx_nonzero_u), hist_struc.nbins, 1));
hist_struc.xcounts = hist_struc.xcounts + counts;

%Update summary statistics
curr_contents = hist_struc.contents;
hist_struc.contents = hist_struc.contents + sum(weights(idx_valid));
hist_struc.entries = hist_struc.entries + sum(idx_valid);
hist_struc.mean = (hist_struc.mean * curr_contents +...
    sum(weights(idx_valid).*(data(idx_valid)))) / hist_struc.contents;

hist_struc.mean2 = (hist_struc.mean2 * curr_contents + ...
    sum(weights(idx_valid).*data(idx_valid).*data(idx_valid))) / hist_struc.contents;


