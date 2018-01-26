function [v] = ellipse_general( center, radii, x_axis )
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
%
%Convert to ellipse general form using linear algebra matrix multiplication
if all(size(x_axis)==[2 2])
    R = x_axis;
else
    cc = x_axis(1);
    ss = x_axis(2);
    R = [cc ss 0; -ss cc 0; 0 0 1];
end

x0 = -center(1);
y0 = -center(2);
rx = 1 / radii(1)^2;
ry = 1 / radii(2)^2;

T1 = [1 0 0; 0 1 0; x0 y0 1];
T2 = T1';
T2(3,3) = -1;

E = T1*R*diag([rx ry 1])*R'*T2;
v = [E(1,1) E(2,2) E(1,2) E(1,3) E(2,3) E(3,3)]';

% normalize to the more conventional form with constant term = -1
if abs( v(end) ) > 1e-6
    v = -v / v(end); 
else
    v = -sign( v(end) ) * v;
end

%%
% %Brute force method, multiplying out terms
% cc = x_axis(1);
% ss = x_axis(2);
% x0 = center(1)*cc + center(2)*ss;
% y0 = center(2)*cc - center(1)*ss;
% 
% a = 1 / radii(1)^2;
% b = 1 / radii(2)^2;
% d = -x0*a;
% e = -y0*b;
% f = a*x0^2 + b*y0^2 - 1;
% 
% a2 = a*cc^2 + b*ss^2;
% b2 = a*ss^2 + b*cc^2;
% c2 = cc*ss*(a-b);
% d2 = d*cc - e*ss;
% e2 = d*ss + e*cc;
% f2 = f;
% 
% 
% v = [a2 b2 c2 d2 e2 f2]';
% 
% % normalize to the more conventional form with constant term = -1
% if abs( v(end) ) > 1e-6
%     v = -v / v(end); 
% else
%     v = -sign( v(end) ) * v;
% end
