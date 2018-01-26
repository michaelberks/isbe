function [gradient_xyz, gradient_mag] = gradient_3d(data_in)
%GRADIENT_ND compute spatial gradients of 3 dimensional arrays,
%outputing normalised directions and magnitudes
%   [gradient_components, gradient_mag] = gradient_nd(data_in)
%
% Inputs:
%      data_in - ny x nx x nz array
%
%
% Outputs:
%      gradient_xyz - ny x nx x nz x 3 array of spatial gradient components
%      dx, dy, dz, normalised so that sum(gradient_xyz(i,j,k,;).^2,4) = 1
%
%      gradient_mag - magnitude of the gradient
%      (sum(gradient_xyz(i,j,k,;).^2,4) prior to normalising)
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 24-Oct-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
[ny nx nz] = size(data_in);
gradient_xyz = zeros(ny, nx, nz, 3);
[   gradient_xyz(:,:,:,1)...
    gradient_xyz(:,:,:,2)...
    gradient_xyz(:,:,:,3)] = gradient(data_in);
gradient_mag(:,:,:) = sqrt(...
    gradient_xyz(:,:,:,1).^2 + ...
    gradient_xyz(:,:,:,2).^2 +...
    gradient_xyz(:,:,:,3).^2 );

gradient_xyz(:,:,:,1) = gradient_xyz(:,:,:,1) ./ ...
    (gradient_mag(:,:,:) + 1e-6);
gradient_xyz(:,:,:,2) = gradient_xyz(:,:,:,2)./...
    (gradient_mag(:,:,:) + 1e-6);
gradient_xyz(:,:,:,3) = gradient_xyz(:,:,:,3)./...
    (gradient_mag(:,:,:) + 1e-6);
