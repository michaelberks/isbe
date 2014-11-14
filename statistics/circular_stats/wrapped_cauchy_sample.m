function [sample, cdf] = wrapped_cauchy_sample(mu, rho, n_pts, cdf)

% Invert the CDF of the wrapped_cauchy by integrating the wrapped_cauchy as values
% across its PDF, then interpolate;

% Can choose to load in the integrated function as look-up table rather than integrating again
% each time

if nargin < 4
    %Create look-up table
    lookup_pts = 1000;
    theta = linspace(0, 2*pi, lookup_pts);
    cdf = zeros(lookup_pts, 1);
    F = @(theta)wrapped_cauchy_pdf(theta, mu , rho);
    
    for ii = 1:lookup_pts
        cdf(ii) = quadl(F,0,theta(ii));
    end
else
    theta = linspace(0, 2*pi, length(cdf));
end
        
        
sample_uni = rand(n_pts, 1);
sample = interp1(cdf, theta, sample_uni);  