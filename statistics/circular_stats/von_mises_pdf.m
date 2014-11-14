function [p] = von_mises_pdf(theta, mu, kappa)
% Defines the probability density function for a random variable theta that
% has a Von Mises distribution M(mu, k);

%Compute the modified bessel function of the first kind and order 0
I_0 = besseli(0, kappa);

p = (1 / (2*pi*I_0))*exp(kappa*cos(theta - mu));
