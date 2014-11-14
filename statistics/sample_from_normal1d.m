function [rand_samp] = sample_from_normal1d(mu, sigma, no_samples)
% Generate a sample from normal (Gaussian) distribution defined by {mu,
% covar}
if nargin < 3
    no_samples = 1;
end

%take uniform sample
rand_samp = mu + sigma * randn(1, no_samples);
