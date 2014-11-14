function [rand_samp] = sample_from_probs(probs, num_samples)
% Generate a sample from binomial (Gaussian) distribution defined by N and
% p. If normal_approx flag is set, will apply the normal approximation when
% suitable criteria are met (Np > 10 and N(1-p) > 10)
probs(isnan(probs)) = 0;

if any(probs(:))
    cdf = cumsum(probs(:));
    cdf = cdf / cdf(end);
    u = rand(num_samples,1);
    rand_samp = zeros(num_samples,1);
    for ii = 1:num_samples
        rand_samp(ii) = sum(cdf < u(ii)) + 1;
    end
else
    warning('ASYM:sample_from_probs', 'All probabilities are zero or NaN, no samples returned (rand_samp=[])');
    rand_samp = [];
end