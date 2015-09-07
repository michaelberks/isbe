function [kernel_dist] = build_2d_kernel_distribution(sample_xy, grid_xy, weights, debug_mode)
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
if ~exist('grid_xy', 'var') || isempty(grid_xy);
    dist_lims = [min(sample_xy) max(sample_xy)];
    dims = [100 100];
    [xx yy] = meshgrid(...
        linspace(dist_lims(1), dist_lims(2), dims(1)),...
        linspace(dist_lims(3), dist_lims(4), dims(2)));
    grid_xy = [xx(:) yy(:)]; clear xx yy;
end
if ~exist('debug_mode', 'var')
    debug_mode = 0;
end

%Work out the number of points in the breast shape and the number of mass
%centres we have
n_pts = size(grid_xy,1);
N = size(sample_xy, 1);

if ~exist('weights', 'var') || isempty(weights)
    weights = ones(N, 1);
end

%--------------------------------------------------------------------------
%Compute the fixed kernel desnity disribution D

%First compute the optimal sigmas for the fixed Gaussian kernels in the x
%and y directions
xi = sample_xy(:,1);
yi = sample_xy(:,2);
xm = mean(xi);
ym = mean(yi);

sigma_x = 1.06*(N^-.2)*sqrt(sum((xi - xm).^2) / (N - 1));
sigma_y = 1.06*(N^-.2)*sqrt(sum((yi - ym).^2) / (N - 1));

%Compute the density at each point in the breast shape
D = zeros(n_pts,1);
for ii = 1:n_pts
    u = grid_xy(ii,1);
    v = grid_xy(ii,2);
    %D is the sum of densities at each point attributed to the kernel
    %fitted at each mass centre
    D(ii) = sum(weights .* exp(-(u-xi).^2 / (2*sigma_x^2)) .* exp(-(v-yi).^2 / (2*sigma_y^2)));
end

%Compute the normalisation constant
A = sum(D);
D = D / A;

%Save the location model information
kernel_dist.D_f = D;
kernel_dist.x = grid_xy(:,1);
kernel_dist.y = grid_xy(:,2);

if debug_mode

    %--------------------------------------------------------------------------
    %Now compute adaptive kernels

    %First compute the density at each mass centre from the fixed kernel
    %distribution
    D_mass = zeros(N,1);
    for ii = 1:N
        D_mass(ii) = sum(exp(-(xi(ii)-xi).^2 / (2*sigma_x^2)) .* exp(-(yi(ii)-yi).^2 / (2*sigma_y^2)));
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
        u = grid_xy(ii,1);
        v = grid_xy(ii,2);
        D_a(ii) = sum(exp(-(u-xi).^2 ./ (2*(Li.*sigma_x).^2)) .* exp(-(v-yi).^2 ./ (2*(Li.*sigma_y).^2)));
    end
    A_a = sum(D_a);
    D_a = D_a / A_a;
    kernel_dist.D_a = D_a;
    
    %If debugging, display the density distributions in the breast shape

%     figure; 
%     subplot(1,2,1); imagesc(reshape(D, dims(2), dims(1))); axis image; colormap(jet(256));
%     subplot(1,2,2); imagesc(reshape(D_a, dims(2), dims(1))); axis image; colormap(jet(256));
end
