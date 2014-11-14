function [x y z] = mb_ellipsoid(centres, radii, transform, n_pts)

%transform should be orthogonal
if nargin < 4
    n_pts = 20;
end

if nargin < 3
    transform = eye(3);
end
if nargin < 2
    radii = [1 1 1];
end
if nargin < 1
    centres = [0 0 0];
end

[x,y,z] = sphere(n_pts);

[xyz] = transform*[x(:) y(:) z(:)]';
x = reshape(xyz(1,:), 21, 21);
y = reshape(xyz(2,:), 21, 21);
z = reshape(xyz(3,:), 21, 21);

x = radii(1)*x + centres(1);
y = radii(2)*y + centres(2);
z = radii(3)*z + centres(3);


