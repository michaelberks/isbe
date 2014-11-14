function [pathPointArray] = generatepath(varargin)
%CONV Convolution and polynomial multiplication.
%   C = CONV(A, B) convolves vectors A and B.  The resulting
%   vector is length LENGTH(A)+LENGTH(B)-1.
%   If A and B are vectors of polynomial coefficients, convolving
%   them is equivalent to multiplying the two polynomials.
%
%   Class support for inputs A,B: 
%      float: double, single
%
%   See also DECONV, CONV2, CONVN, FILTER and, in the Signal
%   Processing Toolbox, XCORR, CONVMTX.

%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 5.16.4.3 $  $Date: 2006/06/20 20:09:21 $

if (nargin==0 && nargout==0), test(); return; end

[nPoints, sigmaX, roundness] = parsearguments(varargin{:});

xPosition0 = 0;
yPosition0 = 0;

xVelocity0 = 1;
yVelocity0 = 0;

xVelocityDecay = 0.1;
yVelocityDecay = roundness * xVelocityDecay;

pathPointArray = zeros(2*nPoints, 2);

xPosition = xPosition0;
xVelocity = xVelocity0;
yPosition = yPosition0;
yVelocity = yVelocity0;
for i = nPoints:-1:1
    pathPointArray(i,:) = [xPosition yPosition];

    [xPosition, xVelocity] = updatex(xPosition, xVelocity, ...
                                     xVelocityDecay, sigmaX);
    [yPosition, yVelocity] = updatey(yPosition, yVelocity, ...
                                     yVelocityDecay);
end

   
% Reverse initial direction for right-hand limb
xPosition = xPosition0;
xVelocity = -xVelocity0;
yPosition = yPosition0;
yVelocity = -yVelocity0;
for i = nPoints+1:2*nPoints
    pathPointArray(i,:) = [xPosition yPosition];

    [xPosition, xVelocity] = updatex(xPosition, xVelocity, ...
                                     xVelocityDecay, sigmaX);
    [yPosition, yVelocity] = updatey(yPosition, yVelocity, ...
                                     yVelocityDecay);
end



function [nPoints, sigmaX, roundness] = ...
    parsearguments(nPoints, sigmaX, roundness)

if ~exist('nPoints','var') || isempty(nPoints)
    nPoints = 100; 
end

if ~exist('sigmaX','var') || isempty(sigmaX)
    sigmaX = 0.025; 
end

if ~exist('roundness','var') || isempty(roundness)
    roundness = 1.0; 
end



function [xPosition, xVelocity] = updatex(xPosition, xVelocity, ...
                                          xVelocityDecay, sigmaX)

xVelocityTarget = 0;
xVelocity = xVelocity + ...
            xVelocityDecay  * (xVelocityTarget - xVelocity) + ...
            sigmaX * randn;

xPosition = xPosition + xVelocity;



function [yPosition, yVelocity] = updatey(yPosition, yVelocity, ...
                                          yVelocityDecay)
xVelocityTarget = -1;
yVelocity = yVelocity + ...
            yVelocityDecay * (xVelocityTarget - yVelocity);

yPosition = yPosition + yVelocity;



%% Test script
function test()
clc;

pathPointArray = generate_path();
figure(1); clf; hold on;
    plot(pathPointArray(:,1), pathPointArray(:,2), 'k.-');
    plot(pathPointArray(1,1), pathPointArray(1,2), 'go');
    plot(pathPointArray(end,1), pathPointArray(end,2), 'ro');
    axis('equal', 'ij');
    
    