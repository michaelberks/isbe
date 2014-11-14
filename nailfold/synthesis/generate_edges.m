function [widthArray, directionArray, innerArray, outerArray] = ...
    generate_edges(pathPointArray, varargin)

if (nargin==0 && nargout==0), test(); return; end

[sigmaW, widthAtApex, widthRatio] = parsearguments(varargin{:});

widthThick = widthAtApex * sqrt(widthRatio);
widthThin  = widthAtApex / sqrt(widthRatio);

widthVelocityDecay = 0.1;

nPoints = ceil(size(pathPointArray,1) / 2);
widthArray = zeros(2*nPoints, 1);

% Do left hand limb (thicker one)
% Positive velocity makes it grow
width = widthAtApex;
widthVelocity = 0;
for i = nPoints+1:2*nPoints
    widthArray(i) = width;

    widthVelocity = 0.5 * widthVelocity + ...
                    0.5 * widthVelocityDecay * (widthThick - width) + ...
                    sigmaW * randn;
    
    width = width + widthVelocity;
end

% Do right hand limb (thinner one)
% Negative velocity makes it shrink
width = widthAtApex;
widthVelocity = 0;
for i = nPoints:-1:1
    widthArray(i) = width;
                
    width = width + widthVelocity;

    widthVelocity = 0.5 * widthVelocity + ...
                    0.5 * widthVelocityDecay * (widthThin - width) + ...
                    sigmaW * randn;
end

directionArray = vesselDirection(pathPointArray);

[innerArray, outerArray] = ...
    vesselBoundaries(pathPointArray, directionArray, widthArray);



function directionArray = vesselDirection(pathPointArray)

nPoints = size(pathPointArray,1);
directionArray = nan(nPoints, 2);

directionArray(1,:) = (pathPointArray(2,:) - pathPointArray(1,:));
for i = 2:nPoints-1
    directionArray(i,:) = (pathPointArray(i+1,:) - pathPointArray(i-1,:));
end
directionArray(end,:) = (pathPointArray(end,:) - pathPointArray(end-1,:));



function [innerArray, outerArray] = ...
    vesselBoundaries(pathPointArray, directionArray, widthArray)

nPoints = size(pathPointArray,1);
innerArray = nan(nPoints, 2);
outerArray = nan(nPoints, 2);

for i = 1:nPoints
    directionArray(i,:) = directionArray(i,:) / norm(directionArray(i,:));
    normalDirection = [-directionArray(i,2) directionArray(i,1)];
    
    innerArray(i,:) = pathPointArray(i,:) + normalDirection * widthArray(i)/2;
    outerArray(i,:) = pathPointArray(i,:) - normalDirection * widthArray(i)/2;
end



function [sigmaW, widthAtApex, widthRatio] = ...
    parsearguments(sigmaW, widthAtApex, widthRatio)

if ~exist('sigmaW','var') || isempty(sigmaW)
    sigmaW = 0.2; 
end

if ~exist('widthAtApex','var') || isempty(widthAtApex)
    widthAtApex = 8; 
end

if ~exist('widthRatio','var') || isempty(widthRatio)
    widthRatio = 2; 
end



%% Test script
function test()
clc;

pathPoints = generate_path();

[widthArray, directionArray, innerArray, outerArray] = ...
    generate_edges(pathPoints, 0.1, 8, 2.0);

figure(1); clf; hold on;
    plot(pathPoints(:,1), pathPoints(:,2), '.-', 'color', 0.75*[1,1,1]);
    plot(pathPoints(1,1), pathPoints(1,2), 'go');
    plot(pathPoints(end,1), pathPoints(end,2), 'ro');
    
    plot(innerArray(:,1), innerArray(:,2), 'k.-');
    plot(outerArray(:,1), outerArray(:,2), 'k.-');
    axis('equal','ij');
