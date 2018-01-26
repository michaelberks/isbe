function [sigma2] = rician_noise_estimates(local_stats, estimate_type, n_hist_bins, n_kdf_bins, debug_mode)
%RICIAN_NOISE_ESTIMATES *Insert a one line summary here*
%   [sigma2] = rician_noise_estimates(local_stats, estimate_type)
%
% Inputs:
%      local_stats - *Insert description of input variable here*
%
%      estimate_type - *Insert description of input variable here*
%
%
% Outputs:
%      sigma2 - *Insert description of input variable here*
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

%Default to local mean estimator
if ~exist('estimate_type', 'var') || isempty(estimate_type)
    estimate_type = 'local_mean';
end
if ~exist('n_hist_bins', 'var') || isempty(n_hist_bins)
    n_hist_bins = 20;
end
if ~exist('n_kdf_bins', 'var') || isempty(n_kdf_bins)
    n_kdf_bins = 20;
end
if ~exist('debug_mode', 'var') || isempty(debug_mode)
    debug_mode = false;
end

%Estimate mode of local stats
local_stats_mode = sample_mode(local_stats, n_hist_bins, n_kdf_bins, debug_mode);

switch estimate_type
    
    case 'local_mean'
        sigma2 = (2/pi)*local_stats_mode^2;
    case 'local_x2'
        sigma2 = local_stats_mode / 2;
    case 'local_var_bg'
        sigma2 = 2*local_stats_mode / (4 - pi);
    case 'local_var_im'
        sigma2 = local_stats_mode;
end
