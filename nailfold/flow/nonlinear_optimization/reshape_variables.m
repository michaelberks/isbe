function p = reshape_variables(p)
% Reshape the variables to correct sizes.

[m, n, nImages] = size(p.imgStack);

if isfield(p, 'uu')
    p.uu = reshape(p.uu, [m, n]);
    p.vv = reshape(p.vv, [m, n]);
else
    p.F = reshape(p.F, [m, n]);
    p.theta = reshape(p.theta, [m, n]);
end

if isfield(p, 'Q')
    p.Q = reshape(p.Q, [m, n]);
end

if isfield(p, 'f')
    p.f = reshape(p.f, [nImages, 1]);
end

if isfield(p, 'displacements')
    p.displacements = reshape(p.displacements, [nImages-1, 2]);
end