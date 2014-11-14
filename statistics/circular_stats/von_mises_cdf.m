function [P] = von_mises_cdf(theta, mu, kappa, n)
% Defines the probability density function for a random variable theta that
% has a Von Mises distribution M(mu, k);

if nargin < 4
    n = 1e2;
end

%Compute the modified bessel function of the first kind and order 0
I_0 = besseli(0, kappa);

thetas = linspace(theta(1), theta(2), n);

p = (1 / (2*pi*I_0))*exp(kappa*cos(thetas - mu));

P = sum(p) * (theta(2)-theta(1)) / n;
