function indices = sample_from_pvec(p_vec, N)
% Given a discrete histogram (a vector of probabilities), draw samples from
% the distribution and return a vector of indices

if nargin<2, N = 1; end

f_debug = (nargin==0 && nargout==0);
if f_debug
	p_vec = rand(1,3);
	p_vec = p_vec/sum(p_vec);
	N = 100000;
end

% normalize score vector and compute cumulative probability
k = length(p_vec);
p_vec = p_vec / sum(p_vec);
cp_vec = cumsum(p_vec);

% sample from uniform distribution
samples = rand(N,1);

% find the bin into which each sample falls
indices = sum( samples(:,ones(1,k))>=cp_vec(ones(N,1),:), 2 ) + 1;
 
% figure(10); clf; hold on;
% 	plot(p_vec);
% 	plot(indices(1),p_vec(indices(1)),'r+');

if f_debug
	clc;
	disp(p_vec);
	disp(hist(indices,1:k)/N);
	clear;
end
