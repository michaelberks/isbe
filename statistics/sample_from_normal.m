function [rand_samp] = sample_from_normal(mu, covar, no_samples)
% Generate a sample from normal (Gaussian) distribution defined by {mu,
% covar}
if nargin < 3
    no_samples = 1;
end

dim = size(covar, 1);
% If m is a column vector, change to a row vector
if size(mu, 2) ~= dim
    mu = mu';
end

if dim > 1
    %e_vecs and e_vals are both dim x dim matrices
    [e_vecs, e_vals] = eig(covar);
else
    e_vecs = 1;
    e_vals = covar;
end
e_vecs = e_vecs*sign(e_vals);
e_vals = e_vals*sign(e_vals);

%take uniform sample
uni_samp = randn(dim, no_samples);

%multiply by eigen vectors and values and translate by mean
%rand_samp = (e_vecs*(sqrtm(e_vals) * uni_samp))' + repmat(mu, no_samples, 1);
rand_samp = bsxfun(@plus, (e_vecs*(sqrtm(e_vals) * uni_samp))', mu); %Use bsxfun as slightly faster than repmat
