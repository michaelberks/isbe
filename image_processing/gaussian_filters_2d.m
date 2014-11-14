function [filters] = gaussian_filters_2d(varargin)
% aux function to make the 2d gaussian derivative filters

f_debug = (nargin==0 && nargout==0);
if f_debug, test_script(); return; end

[filters] = func(varargin{:});


%% The function
function [filters] = func(n_angles, sigma, width)

if ~exist('n_angles','var'), n_angles = 3; end

if ~exist('sigma','var')
    % use default sigma and width
	[g,dg,ddg] = gaussian_filters_1d(); 
elseif ~exist('width','var')
    % use supplied sigma and default width
	[g,dg,ddg] = gaussian_filters_1d(sigma);
else
    % use supplied sigma and supplied width
	[g,dg,ddg] = gaussian_filters_1d(sigma,width);
end

% define 2D filters in terms of 1D filters
Gxx = g'*ddg;
Gxy = dg'*dg;
Gyy = ddg'*g;

filters = zeros(size(Gxx,1),size(Gxx,2),n_angles);
for i = 1:n_angles
	theta = pi * (i-1)/n_angles;
	filters(:,:,i) = cos(theta)^2 * Gxx + ...
					 sin(theta)^2 * Gyy + ...
					 sin(2*theta) * Gxy;
end


%% Test script
function test_script()
clc;

n_angles = 3;
sigma = 3;

filters = func(n_angles, sigma);

figure(); colormap(gray(256));
for i = 1:n_angles
    subplot(1,n_angles,i);
    imagesc(filters(:,:,i));
    axis('image');
end
