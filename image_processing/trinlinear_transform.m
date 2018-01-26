function [xyz_t] = trinlinear_transform(xyz, tri_r, tri_t)
%TRINLINEAR_TRANSFORM *Insert a one line summary here*
%   [xyz_t] = trinlinear_transform(xyz, tri_r, tri_t)
%
% Inputs:
%      xyz - *Insert description of input variable here*
%
%      tri_r - *Insert description of input variable here*
%
%      tri_t - *Insert description of input variable here*
%
%
% Outputs:
%      xyz_t - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 09-Mar-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if size(xyz,2) ~= 3
    error('xyz must be n x 3 array of 3-D corrdinates');
end

if size(tri_r,1) ~= 3 || size(tri_r,2) ~= 8
    error('tri_r must be 3 x 8 array');
end

if length(tri_t) ~= 3
    error('tri_t must be 3 element vector');
end

%Shift input to centre of rotation
cx = xyz(:,1) - tri_t(1);
cy = xyz(:,2) - tri_t(2);
cz = xyz(:,3) - tri_t(3);

xyz_t = zeros(size(xyz));
for i_d = 1:3
    xyz_t(:,i_d) =... %Rotation
        cx .* tri_r(i_d,1) + ...
        cy .* tri_r(i_d,2) + ...
        cz .* tri_r(i_d,3) + ...
        cx.*cy .* tri_r(i_d,4) + ...
        cx.*cz .* tri_r(i_d,5) + ...
        cy.*cz .* tri_r(i_d,6) + ...
        cx.*cy.*cz .* tri_r(i_d,7) + ...
        tri_r(i_d,8) + ... %Translation
        tri_t(i_d); %Relocation
end
