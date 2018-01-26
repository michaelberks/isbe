function [kernel_dist] = build_1d_kernel_distribution(sample_x, grid_x, do_adaptive, debug_mode)
%BUILD_MASS_DISTRIBUTION_MODEL *Insert a one line summary here*
%   [mass_centres] = build_mass_distribution_model(mean_shape,mass_list)
%
% Inputs:
%      mean_shape- *Insert description of input variable here*
%
%      mass_centres- *Insert description of input variable here*
%
%
% Outputs:
%      mass_centres- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 17-Nov-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

%Set default parameters
if ~exist('do_adaptive', 'var') || isempty(do_adaptive);
    do_adaptive = 0;
end
if ~exist('grid_x', 'var') || isempty(grid_x);
    grid_x = linspace(min(sample_x * 1.1), max(sample_x * 1.1), 100);
end
if ~exist('debug_mode', 'var')
    debug_mode = 0;
end
if debug_mode
    do_adaptive = 1;
end

%Work out the number of points in the breast shape and the number of mass
%centres we have
n_pts = length(grid_x);
N = length(sample_x);

%--------------------------------------------------------------------------
%Compute the fixed kernel desnity disribution D

%First compute the optimal sigmas for the fixed Gaussian kernels in the x
%and y directions
xi = sample_x(:,1);
xm = mean(xi);
sigma_x = 1.06*(N^-.2)*sqrt(sum((xi - xm).^2) / (N - 1));

%Compute the density at each point in the breast shape
D = zeros(n_pts,1);
for ii = 1:n_pts
    u = grid_x(ii);

    %D is the sum of densities at each point attributed to the kernel
    %fitted at each mass centre
    D(ii) = sum(exp(-(u-xi).^2 / (2*sigma_x^2)));
end

%Compute the normalisation constant
A = sum(D);
D = D / A;

kernel_dist.D_f = D;
kernel_dist.x = grid_x;

%--------------------------------------------------------------------------
%Now compute adaptive kernels
if do_adaptive
    %First compute the density at each mass centre from the fixed kernel
    %distribution
    D_mass = zeros(N,1);
    for ii = 1:N
        D_mass(ii) = sum(exp(-(xi(ii)-xi).^2 / (2*sigma_x^2)));
    end
    D_mass = D_mass / A;
    %D_mass = interp2(uu, vv, D, xi, yi, 'linear');

    %Compute the geometric mean of the initial density estimate
    g = exp( sum(log(D_mass)) / N );

    %Compute the scaling factor to adapt the kernel fitted at each mass centre
    alpha = .5; %make alpha an optional input argument
    Li = (D_mass / g) .^ -alpha;

    %Compute the adaptive kernel density and normalise
    D_a = zeros(n_pts,1);
    for ii = 1:n_pts
        u = grid_x(ii);
        D_a(ii) = sum(exp(-(u-xi).^2 ./ (2*(Li.*sigma_x).^2)));
    end
    A_a = sum(D_a);
    D_a = D_a / A_a;

    %Save the location model information
    kernel_dist.D_a = D_a;
end

%If debugging, display the density distributions in the breast shape
if debug_mode
    figure; hold all; 
    plot(grid_x, D, 'linewidth', 2);
    plot(grid_x, D_a, 'linewidth', 2);
    title('Kernel smoothed distributions');
    xlabel('X');
    ylabel('k(X)');
    legend({'Fixed kernel size', 'Adaptive kernel size'});
end
