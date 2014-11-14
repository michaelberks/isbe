function value = sample_rand_exp(range, decay_rate, dims)

if nargin < 3
    dims = 1;
end

mu = (range(2) - range(1)) / (2*log(decay_rate));
value = range(1) + exp_rand_sample(mu, dims);