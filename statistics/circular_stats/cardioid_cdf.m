function [C] = cardioid_cdf(theta, mu, rho)
% Defines the cumulative distribution function for a random variable theta that
% has a Cardioid distribution C(mu, rho);

C = theta/(2*pi) + (rho/pi)*(sin(theta - mu) + sin(mu));