function [filters] = gaussian_filters_2d_hilbert(varargin)
% Auxiliary function to make the 2d gaussian derivative filters

f_debug = (nargin == 0 && nargout == 0);
if f_debug, test_script(); return; end

[filters] = func(varargin{:});


%% The function
function [filters] = func(n_angles, sigma, width)

if ~exist('n_angles','var'), n_angles = 4; end
if ~exist('sigma','var'), sigma = 1; end;
if ~exist('width','var'), width = 5*sigma; end

filters = zeros(2*width+1, 2*width+1, n_angles);
for i = 1:n_angles
	theta = pi * (i-1)/n_angles;
    filters(:,:,i) = flt(width, theta, sigma);
end

function f = flt(scale, theta, sigma)

% flip y coordinate as we're dealing with images
[mx,my] = meshgrid(-scale:scale, -scale:scale);
xx = cos(theta)*mx - sin(theta)*my;
yy = -sin(theta)*mx - cos(theta)*my;

% Formula for the kernel taken from Appendix G of:
%   "The design and use of steerable filters",
%   W. T. Freeman and E. H. Adelson,
%   PAMI 13(9), September 1991
%
% To get to canonical form, x -> x / sqrt(2*sigma*sigma)
%                           y -> y / sqrt(2*sigma*sigma)

sxx = xx / sqrt(2*sigma*sigma);
syy = yy / sqrt(2*sigma*sigma);
kernel = -2.205*sxx + 0.9780*sxx.^3;
gaussian_envelope = exp(-(sxx.^2 + syy.^2));
f = kernel .* gaussian_envelope;

% % This is the form given in Ayres and Rangayyan, JEI 2007 - it's wrong
% sxx = xx / sigma;
% kernel = -2.205*sxx + 0.9780*sxx.^3;
% gaussian_envelope = exp(-(xx.^2 + yy.^2) / (2*sigma*sigma));
% f = (kernel / sigma) .* gaussian_envelope;


%% Test script
function test_script()
clc;

n_angles = 4;
sigma = 10;

filters = gaussian_filters_2d_hilbert(n_angles, sigma);

figure(1); colormap(gray(256));
for i = 1:n_angles
    subplot(1,n_angles,i);
    imagesc(filters(:,:,i)); axis('image','ij');
    
    display(rank(filters(:,:,i)));
end

r = rank(filters(:,:,2));
[u,s,v] = svd(filters(:,:,2));

figure(); clf; colormap(gray(256));
for i = 1:r
    img = u(:,i)*v(:,i)';
    subplot(2,2,i);
    imagesc(img); axis('image','ij');
end

