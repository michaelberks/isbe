function prob_vec = ncm_parzen(centres, width, rng, weights)

if (nargin==0 && nargout==0), test_script(); return; end

centres = centres(:);
n_centres = size(centres, 1);

if ~exist('weights','var'), weights = ones(n_centres, 1); end

weights = weights / sum(weights(:));

prob_vec = zeros(1, length(rng));
for i = 1:n_centres
    prob_vec = prob_vec + ...
               weights(i) * one_hump(centres(i,1), width, rng);
end


%% Get a single parzen window hump
function y = one_hump(mean, sigma, x)

sigma_sqr = sigma*sigma;

y = exp(-0.5 * (x-mean).*(x-mean) / sigma_sqr);
y = y / sqrt(2 * pi * sigma_sqr);


%% Test script
function test_script()

clc;

centres = rand(1,30);
xrng = -1:0.01:2;
pv = ncm_parzen(centres, 0.1, xrng);

figure(1); clf;
    plot(xrng, pv, 'b-');
    