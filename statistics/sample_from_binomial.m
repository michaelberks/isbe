function [binom_samp] = sample_from_binomial(N, p, num_samples, normal_approx)
% Generate a sample from binomial (Gaussian) distribution defined by N and
% p. If normal_approx flag is set, will apply the normal approximation when
% suitable criteria are met (Np > 10 and N(1-p) > 10)

if nargin < 4
    normal_approx = true;
end

%Check for normal approximation
if normal_approx && (N*p > 10) && (N*(1-p) > 10)
    binom_samp = round(sample_from_normal(N*p, N*p*(1-p), num_samples));
else
    %perform bernoullli tests to get proper binomail distribution
    binom_samp = sum(rand(num_samples,N) < p, 2);
end
