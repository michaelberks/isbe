function y = complex_interp2(c, x_pts, y_pts, w, interpmethod)

% function y = cpxinterp2(c,pts,w,interpmethod)
%
% Complex interpolation routine for the matrix c.
% Uses bandpass interpolation centred on w(1:2) rad/sample centre freq.
% pts specifies the interpolation points between the integer sample
% positions of c. (eg pts = [-1 1]/4  doubles the rate in each direction and 
% pts = [-1 0 1]/3 triples the rate, with uniform sampling.)
% w specifies the approximate expected rotation in radians between consecutive samples.
% w(1) is the rate of rotation down the columns, and w(2) across the rows.
% of c, so that phase differences can be unwrapped correctly.
% interpmethod is a string which specifies the rule used by interp2:
%   'nearest','linear','cubic','spline'  are valid methods in
% increasing order of computational complexity (defaults to 'cubic').
% 'linear' is about 4 times as fast as 'cubic'.
%
% Nick Kingsbury, Cambridge University, October 2004

if nargin < 4, interpmethod = 'cubic'; end

nr = size(c,1);
nc = size(c,2);
nx_pts = length(x_pts);
ny_pts = length(y_pts);
jw = sqrt(-1) * w;

% Set up matrix padding vectors.
z1 = zeros(1,nc+2);
z2 = zeros(nr,1);

% % Create linear ramp matrices for phase wrapping.
thx = repmat(1:nc, nr, 1);
thy = repmat((1:nr)', 1, nc);

% Create matrices of interpolated point locations in the output matrix.
xs = repmat(x_pts(:)', nx_pts, 1); 
ys = repmat(y_pts(:), 1, ny_pts);
% xs = kron(thx,ones(ni,ni)) + kron(ones(size(thx)),ones(ni,1)*pts.');
% ys = kron(thy,ones(ni,ni)) + kron(ones(size(thy)),pts*ones(1,ni));
   
% Get the image and unwrap the phases and extend it by one row/column at each edge.
ye = [z1; z2  c.*exp(-thx*jw(2)-thy*jw(1))  z2; z1];

% Interpolate ye to the new points, specified in (xs,ys).
yi = interp2(ye,xs+1,ys+1,interpmethod);  % Add one to xs,ys to allow for extension of ye.

% Rewrap the phases.
y = yi .* exp(xs*jw(2)+ys*jw(1));

return;
