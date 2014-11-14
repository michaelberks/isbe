function [sample, cdf] = cardioid_sample(mu, rho, n_pts, cdf)

% Invert the CDF of the von mises by integrating the von Mises as values
% across its PDF, then interpolate;

% Can choose to load in the integrated function as look-up table rather than integrating again
% each time

if nargin < 4
    %Create look-up table
    lookup_pts = 1000;
    theta = linspace(0, 2*pi, lookup_pts);
    cdf = cardioid_cdf(theta, mu, rho);
else
    theta = linspace(0, 2*pi, length(cdf));
end
        
        
sample_uni = rand(n_pts, 1);
sample = interp1(cdf, theta, sample_uni);  