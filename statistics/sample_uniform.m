function value = sample_uniform(range, dims)
if nargin < 2
    dims = 1;
end
value = range(1) + (range(2) - range(1))*rand(dims);