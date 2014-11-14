function [filters] = gabor_filters(varargin)
% Auxiliary function to make the 2d gabor filters

f_debug = (nargin == 0 && nargout == 0);
if f_debug, test_script(); return; end

[filters] = func(varargin{:});


%% The function
function [filters] = func(n_angles, sigma_x, sigma_y, ...
                          tau, halfwidth, remove_bias)

if ~exist('n_angles','var'), n_angles = 4; end
if ~exist('sigma_x','var'), sigma_x = 1; end;

if ~exist('sigma_y','var'), sigma_y = sigma_x; end;
if ~exist('tau','var'), tau = sigma_x*2.35; end
if ~exist('halfwidth','var'), halfwidth = 5 * max(sigma_x, sigma_y); end
if ~exist('remove_bias','var'), remove_bias = true; end

filters = zeros(2*halfwidth+1, 2*halfwidth+1, n_angles);
for i = 1:n_angles
	theta = (i-1) * pi/n_angles;
    filters(:,:,i) = flt(halfwidth, theta, sigma_x, sigma_y, ...
                         tau, remove_bias);
end


function [gabor] = flt(halfwidth, theta, sx, sy, tau, remove_bias)

[mx,my] = meshgrid(-halfwidth:halfwidth, -halfwidth:halfwidth);

% flip y coordinate as we're dealing with images
xx = cos(theta)*mx - sin(theta)*my;
yy = -sin(theta)*mx - cos(theta)*my;

k	= 1/(sqrt(2*pi)*sx*sy); %Normalising factor with extra sigma
%k	= 1/sqrt(2*pi*sx*sy); %'Correct' normalising factor

gaussian_envelope = k * ...
                    exp(-0.5*((xx.^2 / sx^2) + (yy.^2 / sy^2)));

real_filter = gaussian_envelope .* cos(2*pi*xx/tau);
imag_filter = gaussian_envelope .* sin(2*pi*xx/tau);

if remove_bias
    real_bias = mean(real_filter(:));
    real_filter = real_filter - real_bias;
end

gabor = complex(real_filter, imag_filter);


%% Test script
function test_script()
clc; 

n_angles = 6;
sigma = 8;
% tau = 20;

filters = gabor_filters(n_angles, sigma);

figure; colormap(gray(256));
for i = 1:n_angles
    subplot(floor(sqrt(n_angles)),ceil(n_angles/floor(sqrt(n_angles))),i);
    imagesc(imag(filters(:,:,i))); axis('image','ij','off');
end
figure; colormap(gray(256));
for i = 1:n_angles
    subplot(floor(sqrt(n_angles)),ceil(n_angles/floor(sqrt(n_angles))),i);
    imagesc(real(filters(:,:,i))); axis('image','ij','off');
end
