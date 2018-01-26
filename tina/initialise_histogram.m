function [hist_struc] = initialise_histogram(data, bins, bin_lims, num_bins, bin_width, weights)
%INITIALISE_HISTOGRAM Set up a histogram structure - optionally adding data
%   [hist_struc] = initialise_histogram(data, bins, bin_lims, bin_width, num_bins)
%
% Inputs:
%      data - *Insert description of input variable here*
%
%      bins - *Insert description of input variable here*
%
%      bin_lims - *Insert description of input variable here*
%
%      bin_width - *Insert description of input variable here*
%
%      num_bins - *Insert description of input variable here*
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

%Check if we've got data
if exist('data', 'var') && ~isempty(data) && isnumeric(data)
    valid_data = true;
    data = double(data(:));
else
    valid_data = false;
end

%Set up bins
if ~exist('bins', 'var') || isempty(bins)

    %If bin limits not supplied, we must have data
    if ~exist('bin_lims', 'var') || length(bin_lims) ~= 2
        if ~valid_data 
            error('If histogram bins are not explicitly supplied, bin_lims must be provided as a 2-element vector');
        else
            bin_lims = [min(data) max(data)];
        end
    end
        
    %Now compute bins either from the number of bins, the bin widths or
    %automatically
    if exist('num_bins', 'var') && length(num_bins) == 1
        bin_width = (bin_lims(2) - bin_lims(1)) / (num_bins - 1);
        
    elseif exist('bin_width', 'var') && length(bin_width) == 1
        num_bins = (bin_lims(2) - bin_lims(1)) / bin_width + 1;
        
    elseif valid_data
        %Auto estimate number of bins
        preferred_bin_width = 2*iqr(data(:))*numel(data)^(-1/3);
        num_bins = round((bin_lims(2) - bin_lims(1)) / preferred_bin_width + 1);
        bin_width = (bin_lims(2) - bin_lims(1)) / (num_bins - 1);
    else
        error('If num_bins or bin_width not supplied, valid data is required to automatically compute number of bins')
    end
    bins = linspace(bin_lims(1), bin_lims(2), num_bins);
else
    bin_lims = [min(bins) max(bins)];
    bin_width = bins(2) - bins(1);
end

hist_struc.xmin = bin_lims(1);
hist_struc.xmax = bin_lims(2);
hist_struc.ymin = NaN;
hist_struc.ymax = NaN;

hist_struc.nbins = num_bins;
hist_struc.xbins = bins;
hist_struc.ybins = [];

hist_struc.xincr = bin_width;
hist_struc.yincr = NaN;

%Now set up the data terms
hist_struc.mean = 0;
hist_struc.mean2 = 0;

hist_struc.entries = 0;
hist_struc.contents = 0;
hist_struc.under = 0;
hist_struc.over = 0;
hist_struc.above = 0;
hist_struc.below = 0;
hist_struc.xcounts = 0;
hist_struc.ycounts = 0;

%If we have data add it to the sample
if valid_data
    n = numel(data);
    if exist('weights', 'var') 
        if numel(weights) ~= n
            error('If weights are supplied they must match the size of the input data');
        end
    else
        weights = ones(n,1);
    end
    [hist_struc] = add_histogram_data(hist_struc, data, weights);
end
    
