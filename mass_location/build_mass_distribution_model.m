function [location_model] = build_mass_distribution_model(mean_shape, mass_centres, debug_mode)
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
if nargin < 3
    debug_mode = 1;
end

%Get start and end (x,y)-coordinates to bound the mean breast shape
xs = floor(min(mean_shape(:,1)));
xe = ceil(max(mean_shape(:,1)));
ys = floor(min(mean_shape(:,2)));
ye = ceil(max(mean_shape(:,2)));

mask_size(1) = ye - ys + 1;
mask_size(2) = xe - xs + 1;

%Build a BW mask of the shape
[bw] = poly2mask(...
    mean_shape(:,1) - xs + 1,...
    mean_shape(:,2) - ys + 1,...
    mask_size(1),...
    mask_size(2));


%Extract a list of indices belonging to the breast shape in the mask
[in_idx] = find(bw);

%Convert to (x,y) coordinates and translate back to the axes of the mean
%shape - points in the breast shape are (u(i),v(i))
[vv uu] = find(bw); clear bw;
uu = uu - 1 + xs;
vv = vv - 1 + ys;

%Work out the number of points in the breast shape and the number of mass
%centres we have
n_pts = length(uu);
N = size(mass_centres, 1);

%--------------------------------------------------------------------------
%Compute the fixed kernel desnity disribution D

%First compute the optimal sigmas for the fixed Gaussian kernels in the x
%and y directions
xi = mass_centres(:,1);
yi = mass_centres(:,2);
xm = mean(xi);
ym = mean(yi);

sigma_x = 1.06*(N^-.2)*sqrt(sum((xi - xm).^2) / (N - 1));
sigma_y = 1.06*(N^-.2)*sqrt(sum((yi - ym).^2) / (N - 1));

%Compute the density at each point in the breast shape
D = zeros(n_pts,1);
for ii = 1:n_pts
    u = uu(ii);
    v = vv(ii);
    %D is the sum of densities at each point attributed to the kernel
    %fitted at each mass centre
    D(ii) = sum(exp(-(u-xi).^2 / (2*sigma_x^2)) .* exp(-(v-yi).^2 / (2*sigma_y^2)));
end

%Compute the normalisation constant
A = sum(D);
D = D / A;

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
    u = uu(ii);
    v = vv(ii);
    D_a(ii) = sum(exp(-(u-xi).^2 ./ (2*(Li.*sigma_x).^2)) .* exp(-(v-yi).^2 ./ (2*(Li.*sigma_y).^2)));
end
A_a = sum(D_a);
D_a = D_a / A_a;

%If debugging, display the density distributions in the breast shape
if debug_mode
    D_map = zeros(mask_size);
    D_map(in_idx) = D;
    
    D_map_a = zeros(mask_size);
    D_map_a(in_idx) = D_a;
    figure; 
    
    subplot(1,2,1); imagesc(D_map); axis image; colormap(jet(256));
    subplot(1,2,2); imagesc(D_map_a); axis image; colormap(jet(256));
end

%Save the location model information
location_model.D_f = D;
location_model.D_a = D_a;
location_model.x = uu;
location_model.y = vv;
location_model.mean_shape = mean_shape;
