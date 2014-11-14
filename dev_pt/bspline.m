function w = bspline(t, varargin)

if (nargin==0 && nargout==0), test_script(); return; end

w = func(varargin{:});


function w = func(t, knots, order)

for iorder = 0:order    
end


function test_script()
clc;

nKnots = 3;
knots = linspace(0, 1, nKnots);
order = 0;

nPoints = 101;
nSplines = nKnots - order - 1;
t = linspace(0, 1, nPoints);
bt = zeros(nPoints, nSplines);

for i = 1:nPoints
    bt(i,:) = func(t(i), knots, order);
end

plot(bt);
