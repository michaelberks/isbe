function [mle_dist_params] = dist_params_from_hist(hist_bins, hist_counts, dist_name, n_sample_pts, do_opt)
%DIST_PARAMS_FROM_HIST *Insert a one line summary here*
%   [mle_dist] = dist_params_from_hist(hist_bins, hist_counts, dist_name)
%
% Inputs:
%      hist_bins - locations of histogram bine (eg x), vector length n
%
%      hist_counts - count for each histogram bin (eg estimate of p(x)),
%      vector length n.
%
%      dist_name - name of a distribution that can fitted using matlab's
%      stats toolbox function mle see MLE for more details
%
%      n_samples - number of samples used to generate random sample to
%      which mle estimate is fitted
%
%
% Outputs:
%      mle_dist - Maximum likelihood of parameters for chosen distribution
%
%
% Example:
%
% Notes: MLE
%
% See also:
%
% Created: 19-Dec-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 


%If number of samples not specified use square of histogram length
if ~exist('n_samples', 'var') || isempty(n_sample_pts)
    n_sample_pts = length(hist_bins)^2;
end

%If number of samples not specified use square of histogram length
if ~exist('do_opt', 'var') || isempty(do_opt)
    do_opt = 0;
end

%Generate random sample using empirical distribution fucntion (ie the
%cumulative sum of the hisogram histogram counts)
edf = cumsum(hist_counts);

%Draw a uniformly distributed random sample then use 1d interpolation
%against cdf to invert
sample_uni = rand(n_sample_pts, 1);
s_idx = max(1, find(edf >= min(sample_uni), 1)-1);
e_idx = find(edf > max(sample_uni), 1);
if isempty(e_idx)
    e_idx = length(hist_bins);
end
sample_edf = interp1(edf(s_idx:e_idx), hist_bins(s_idx:e_idx), sample_uni);

%Now use Matlab's MLE function to estimate paramaters distribution that
%best fit this data
mle_dist_params = mle(sample_edf, 'distribution', dist_name);

if do_opt
    %Set options for optimisation    
    obj_fun = @(x)pdf_to_hist_counts_sse(x,...
            dist_name, hist_bins(:), hist_counts(:)); %auxilliary variables required to compute objective

    mle_dist_params = fminunc(obj_fun, mle_dist_params); 
end

function sse = pdf_to_hist_counts_sse(x, dist_name, hist_bins, hist_counts)

switch length(x)
    case 1
        p_x = pdf(dist_name, hist_bins, x(1));
    
    case 2
        p_x = pdf(dist_name, hist_bins, x(1), x(2));
        
    case 3
        p_x = pdf(dist_name, hist_bins, x(1), x(2), x(3));    
end

sse =  sum((p_x - hist_counts).^2); 