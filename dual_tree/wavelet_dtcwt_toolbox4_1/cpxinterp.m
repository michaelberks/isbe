function y = cpxinterp(x,pts,w,interpmethod)

% function y = cpxinterp(x,pts,w,interpmethod)
%
% Complex interpolation routine for the matrix x. Perfom 1D interpolation
% along the columns of x
% Uses bandpass interpolation centred on w(1:2) rad/sample centre freq.
% pts specifies the interpolation points between the integer sample
% positions of x. (eg pts = [-1 1]/4  doubles the rate and 
% pts = [-1 0 1]/3 triples the rate, with uniform sampling.)
% w specifies the approximate expected rotation in radians between consecutive samples.
% w is the rate of rotation down the columns
% of x, so that phase differences can be unwrapped correctly.
% interpmethod is a string which specifies the rule used by interp:
%   'nearest','linear','cubic','spline'  are valid methods in
% increasing order of computational complexity (defaults to 'cubic').
% 'linear' is about 4 times as fast as 'cubic'.
%
% Nick Kingsbury, Cambridge University, October 2004

if nargin < 4, interpmethod = 'cubic'; end

[nr nc] = size(x);

pts = pts(:);
ni = length(pts);
jw = sqrt(-1) * w;

% Set up matrix padding vectors.
z1 = zeros(1,nc);

% Create linear ramp matrices for phase wrapping down each column.
thy = repmat((1:nr)', 1, nc);

% Create matrices of interpolated point locations in the output matrix.
ys = repmat(kron((1:nr)', ones(ni,1)) + kron(ones(nr,1), pts), 1, nc);

   
% Get the image and unwrap the phases and extend it by one row at each edge.
ye = [z1; x.*exp(-thy*jw) ; z1];

% Interpolate ye to the new points, specified in xs
yi = interp1(ye, ys(:,1)+1, interpmethod);  % Add one to xs,ys to allow for extension of ye.

% Rewrap the phases.
y = yi .* exp(ys*jw);

return;
