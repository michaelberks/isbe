function new_im = mb_non_maximal_supp(im, orientations, degrees, ...
                                      normals, discard_val, max_size)
% MB_NON_MAXIMAL_SUPP - Suppress non-maximal values
%   @im = non_maximal_supp(im, orientations, degrees)@
%
%   Suppress (set to zero) values in the 2D image @im@ that are not
%   greater than the two values on each side.  The points to compare
%   are the two that are on the normal to the @orientations@ at each
%   point.  The optional parameter @degrees@ says whether the
%   @orientations@ are given in radians (the default) or
%   degrees.
  
if ~exist('degrees','var'), degrees = false; end
if ~exist('normals','var'), normals = false; end
if ~exist('discard_val','var'), discard_val = 0; end
if ~exist('max_size','var'), max_size = 512; end

if degrees
    orientations = orientations * pi / 180;
end

%make meshgrid of (x, y) co-ordinates
[y_lim x_lim] = size(im);
xx = repmat(1:x_lim, y_lim, 1);
yy = repmat((1:y_lim)', 1, x_lim);

% Compute sin and cosine of normal orientations (i.e. orientations shifted
% through pi/2)
% If the orientations are already normals we don't need to rotate by pi/2;
if normals
	ori_c = cos(orientations);
	ori_s = sin(orientations);
else
	ori_c = cos(orientations + pi/2);
	ori_s = sin(orientations + pi/2);
end

%compute number of parts we need to break image into
r_parts = ceil(y_lim / max_size);
c_parts = ceil(x_lim / max_size);

z1 = zeros(y_lim, x_lim);
z2 = zeros(y_lim, x_lim);

%Go through each segment computing z1 and z2
for rp = 1:r_parts
    for cp = 1:c_parts
        
        sr = 1 + (rp-1)*max_size;
        er = min(rp*max_size, y_lim);
        sc = 1 + (cp-1)*max_size;
        ec = min(cp*max_size, x_lim);
        
        %linearly interpolate values at normal coordinates
        %Note the '-y' here because matrix indexing in 'y' direction is inverse to cartesian
        z1(sr:er,sc:ec) = interp2(im, xx(sr:er,sc:ec) + ori_c(sr:er,sc:ec), yy(sr:er,sc:ec) - ori_s(sr:er,sc:ec), '*linear');
        z2(sr:er,sc:ec) = interp2(im, xx(sr:er,sc:ec) - ori_c(sr:er,sc:ec), yy(sr:er,sc:ec) + ori_s(sr:er,sc:ec), '*linear');
    end
end

%discard any point from im that is smaller than either of its two normals 
discard = (im < z1) | (im < z2);
new_im = im;
new_im(discard) = discard_val;
