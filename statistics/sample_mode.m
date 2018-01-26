function [x_mode] = sample_mode(x, n_hist_bins, n_kdf_bins, debug_mode)
%SAMPLE_MODE *Insert a one line summary here*
%   [n_kdf_bins] = sample_mode(x, n_his_bins)
%
% Inputs:
%      x - *Insert description of input variable here*
%
%      n_his_bins - *Insert description of input variable here*
%      n_kdf_bins - *Insert description of input variable here*
%
%
% Outputs:
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 11-Dec-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('n_hist_bins', 'var') || isempty(n_hist_bins)
    n_hist_bins = 20;
end
if ~exist('n_kdf_bins', 'var') || isempty(n_kdf_bins)
    n_kdf_bins = 20;
end
if ~exist('debug_mode', 'var') || isempty(debug_mode)
    debug_mode = false;
end

[hist_counts, hist_bins] = hist(x(:), n_hist_bins);

[~,max_bin] = max(hist_counts);

if max_bin == 1
    xs = 2*hist_bins(1) - hist_bins(2);
    xe = hist_bins(2);
elseif max_bin == n_hist_bins
    xs = hist_bins(end);
    xe = 2*xs - hist_bins(end-1);
else
    xs = hist_bins(max_bin-1);
    xe = hist_bins(max_bin+1);
end

hist_bins_fine = linspace(xs, xe, n_kdf_bins);
include_pts = x(:) > xs & x(:) < xe;
% [kdf] = build_1d_kernel_distribution(x(include_pts), grid_x, 0);
% 
% [~,max_xi] = max(kdf.D_f);
% x_mode = kdf.x(max_xi);
bin_width2 = (hist_bins_fine(2) - hist_bins_fine(1))/2;
[hist_counts_fine_pos] = hist(x(include_pts)+bin_width2, hist_bins_fine);
[hist_counts_fine_neg] = hist(x(include_pts)-bin_width2, hist_bins_fine);
hist_counts_fine = (hist_counts_fine_pos + hist_counts_fine_neg) / 2;
[~,max_bin] = max(hist_counts_fine);
x_mode = hist_bins_fine(max_bin);

if debug_mode
   figure;
   subplot(1,2,1); bar(hist_bins, hist_counts);
   subplot(1,2,2); bar(hist_bins_fine, hist_counts_fine);
   %subplot(1,2,2); plot(kdf.x, kdf.D_f); hold all;
end
    
