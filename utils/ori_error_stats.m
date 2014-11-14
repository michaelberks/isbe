function [err_stats] = ori_error_stats(ori_errors)
%ORI_ERROR_STATS *Insert a one line summary here*
%   [err_stats] = ori_error_stats(ori_errors)
%
% Inputs:
%      ori_errors - *Insert description of input variable here*
%
%
% Outputs:
%      err_stats - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 11-Sep-2012
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

% remove any NaNs
nan_inds = isnan(ori_errors);
if any(nan_inds)
    warning('ori_error:nans', ...
            'Removing %i NaNs from orientations',sum(nan_inds(:)) );
end

% convert error to degrees first
orierr_deg = ori_errors(~nan_inds)*(180/pi);

% number of samples used
err_stats.n_samples = length(orierr_deg);

% mean and standard deviation of error in range (-180..180]
err_stats.mean = mean(orierr_deg);
err_stats.sd = std(orierr_deg);

% stats for absolute error in range [0..180]
orierr_deg = sort(abs(orierr_deg));
err_stats.abs_mean = mean(orierr_deg);
err_stats.abs_median = median(orierr_deg);
err_stats.abs_range = orierr_deg([1,end])';
err_stats.abs_percentiles = ...
    orierr_deg( round([0.01:0.01:1]*length(orierr_deg)) )';
