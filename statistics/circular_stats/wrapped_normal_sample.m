function [sample] = wrapped_normal_sample(mu, rho, n_pts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OR: Sample from standard normal distribution and wrap it

sigma2 = -2*log(rho);
sample_norm = sample_from_normal1d(mu, sigma2, n_pts);

sample = mod(sample_norm, 2*pi);