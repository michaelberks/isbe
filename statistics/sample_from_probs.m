function [rand_samp] = sample_from_probs(probs, num_samples)
% Generate a sample from estimate of pdf
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