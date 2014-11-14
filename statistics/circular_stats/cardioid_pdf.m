function [p] = cardioid_pdf(theta, mu, rho)
% Defines the probability density function for a random variable theta that
% has a Cardioid distribution C(mu, rho);

%Compute the modified bessel function of the first kind and order 0
p = (1 / (2*pi))*( 1 + 2*rho*cos(theta - mu) );