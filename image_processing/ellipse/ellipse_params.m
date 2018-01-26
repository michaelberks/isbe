function [center, radii, e_axes, eigen_signs] = ellipse_params( v )
%
% Fit an ellispoid/sphere/paraboloid/hyperboloid to a set of xyz data points:
%
%   [center, radii, x_axis, pars ] = ellipsoid_fit( X )
%   [center, radii, x_axis, pars ] = ellipsoid_fit( [x y z] );
%   [center, radii, x_axis, pars ] = ellipsoid_fit( X, 1 );
%   [center, radii, x_axis, pars ] = ellipsoid_fit( X, 2, 'xz' );
%   [center, radii, x_axis, pars ] = ellipsoid_fit( X, 3 );
%
% Parameters:
% * X, [x y z]   - Cartesian data, n x 3 matrix or three n x 1 vectors
% * flag         - '' or empty fits an arbitrary ellipsoid (default),
%                - 'xy' fits a spheroid with x- and y- radii equal
%                - 'xz' fits a spheroid with x- and z- radii equal
%                - 'xyz' fits a sphere
%                - '0' fits an ellipsoid with its axes aligned along [x y z] axes
%                - '0xy' the same with x- and y- radii equal
%                - '0xz' the same with x- and z- radii equal
%
% Output:
% * center    -  ellispoid or other conic center coordinates [xc; yc; zc]
% * radii     -  ellipsoid or other conic radii [a; b; c]
% * x_axis     -  the radii directions as columns of the 3x3 matrix
% * v         -  the 10 parameters describing the ellipsoid / conic algebraically: 
%                Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0
% * chi2      -  residual sum of squared errors (chi^2), this chi2 is in the 
%                coordinate frame in which the ellipsoid is a unit sphere.
%
% Author:
% Yury Petrov, Oculus VR
% Date:
% September, 2015
%

% form the algebraic form of the ellipsoid
A = [ v(1) v(3) v(4); ...
      v(3) v(2) v(5); ...
      v(4) v(5) v(6)];
% find the center of the ellipsoid
center = -A( 1:2, 1:2 ) \ v( 4:5 );
% form the corresponding translation matrix
T = eye( 3 );
T( 3, 1:2 ) = center';
% translate to the center
R = T * A * T';
% solve the eigenproblem
[ evecs, evals ] = eig( R( 1:2, 1:2 ) / -R( 3, 3 ) );
radii = sqrt( 1 ./ diag( abs( evals ) ) );
eigen_signs = sign( diag( evals ) );
radii = radii .* eigen_signs;

%Permute the axes so the first column is most closely aligned with
%x axis, the 2nd with y and so on;
[~, axes1] = max(abs(evecs(1,:)));
axes2 = 3 - axes1;
axes_order = [axes1 axes2];
e_axes = evecs(:, axes_order);
radii = radii(axes_order);
eigen_signs = eigen_signs(axes_order);

%Finally, set sign of axes to be positive in usual sense of direction
e_axes = e_axes * diag(sign(diag(e_axes)));






