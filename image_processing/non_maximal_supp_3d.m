function [nms_vol] = non_maximal_supp_3d(vol, grad_xyz, discard_val)
%NON_MAXIMAL_SUPP_3D *Insert a one line summary here*
%   [nms_im] = non_maximal_supp_3d(im, orientations, degrees, normals, discard_val)
%
% Inputs:
%      vol - Input volume, ny * nx * nz
%
%      grad_xyz - 3D gradients of volume, ny * nx * nz * 3
%
%      discard_val - minimum threshold on gradient strength
%
%
% Outputs:
%      nms_vol - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 09-Apr-2018
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

if ~exist('discard_val', 'var') || isempty(discard_val)
    discard_val = 0;
end

[n_y, n_x, n_z] = size(vol);
nms_vol = zeros(n_y, n_x, n_z);

%Loop through each slice
for i_z = 1:n_z
    
    %Make orientation map for slice
    orientations = atan2(grad_xyz(:,:,i_z,2), grad_xyz(:,:,i_z,1));
    
    %Apply nms and save into output volume
    nms_vol(:,:,i_z) = mb_non_maximal_supp(...
        vol(:,:,i_z), orientations, false, true, discard_val);
end

