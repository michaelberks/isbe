function ridge = mb_curvature_sign_change(g1r, g2r, orientations, degrees, normals)
% MB_NON_MAXIMAL_SUPP - Suppress non-maximal values
%   @im = non_maximal_supp(im, orientations, degrees)@
%
%   Suppress (set to zero) values in the 2D image @im@ that are not
%   greater than the two values on each side.  The points to compare
%   are the two that are on the normal to the @orientations@ at each
%   point.  The optional parameter @degrees@ says whether the
%   @orientations@ are given in radians (the default) or
%   degrees.
  
if nargin < 4
    degrees = false;
end
if nargin < 5
    normals = false;
end

if degrees
    orientations = orientations * pi / 180;
end

%make meshgrid of (x, y) co-ordinates
[y_lim x_lim] = size(g1r);
xx = repmat(1:x_lim, y_lim, 1);
yy = repmat((1:y_lim)', 1, x_lim);

% Compute sin and cosine of normal orientations (i.e. orientations shifted
%through pi/2)
% If the orientations are already normals we don't need to rotate by pi/2;
if normals
    theta = 0;
else
    theta = pi/2;
end
ori_c = cos(orientations + theta);
ori_s = -sin(orientations + theta); %Note the '-' here because matrix indexing in 'y' direction is inverse to cartesian

%compute (x, y) coordinates along +/- unit normal vectors
x_interp1 = xx + ori_c;
x_interp2 = xx - ori_c;

y_interp1 = yy + ori_s;
y_interp2 = yy - ori_s;

%linearly interpolate values at normal coordinates
z1 = interp2(g1r, x_interp1, y_interp1);
z2 = interp2(g1r, x_interp2, y_interp2);

%discard any point from im that is smaller than either of its two normals 
ridge = .5 * sign(g2r) .* abs(sign(z1) - sign(z2));