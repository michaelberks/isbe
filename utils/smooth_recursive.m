function out = smooth_recursive(in, sigma)
% Apply recursive gaussian smoothing as described in:
%   "Recursive implementation of the Gaussian filter"
%   Ian T Young and Lucas J van Vliet
%   Signal Processing 44, pp 139-151, 1995

if (nargin==0 && nargout==0), test_script(); return; end
out = func(in, sigma);


%% The function
function out = func(in, sigma)
output_sz = size(in);

if (sigma > 2.5)
    q = 0.98711*sigma - 0.96330;
else
    q = 3.97156 - 4.14554*sqrt(1 - 0.26891*sigma);
end

b0 = 1.57825 + 2.44413*q + 1.4281*q*q + 0.422205*q*q*q;
b1 = (2.44413*q + 2.85619*q*q + 1.26661*q*q*q) / b0;
b2 = -(1.4281*q*q + 1.26661*q*q*q) / b0;
b3 = (0.422205*q*q*q) / b0;
B = 1 - (b1 + b2 + b3);

% Forward pass
in = [zeros(3,1); in(:); zeros(3,1)];
w = zeros(size(in));
for i = 4:length(in)-3
    w(i) = B*in(i) + (b1*w(i-1) + b2*w(i-2) + b3*w(i-3));
end

% Backward pass
out = zeros(size(in));
for i = length(in)-3:-1:4
    out(i) = B*w(i) + (b1*out(i+1) + b2*out(i+2) + b3*out(i+3));
end

out = reshape(out(4:end-3), output_sz);


%% Test script
function test_script()
clc;

hw = 50; n = 0;
f = zeros(1,2*hw+1);
f(hw+1-n:hw+1+n) = 255;

g = smooth_recursive(f, 0.5);
disp([g(1:hw+1); g(end:-1:hw+1)]);

figure(1); clf; hold on;
    plot(-hw:hw, f, 'b-');
    plot(-hw:hw, g, 'r-');
    