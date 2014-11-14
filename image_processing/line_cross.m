function [c xy idx1 idx2] = line_cross(line1, line2)
%LINE_CROSS determin whether two lines specified by discrete coordinates intersect
%   [c] = line_cross(line1, line2)
%
% Inputs:
%      line1 - Nx2 array of xy coordinates defining 1st line
%
%      line2 - Mx2 array of xy coordinates defining 2nd line
%
%
% Outputs:
%      c - true if lines cross, false otherwise
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 27-Jan-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

if (size(line1,1) < 2) || (size(line2,1) < 2)
    error('Both lines must have at least 2 points otherwise they don''t define a line');
end

%Get x and y coords of both lines
x1 = line1(:,1); y1 = line1(:,2);
x2 = line2(:,1); y2 = line2(:,2);

%Get direction vectors of each lines segment
a1 = [-diff(y1) diff(x1)]; 
a2 = [-diff(y2) diff(x2)];

b1 = y1(1:end-1).*diff(x1) - x1(1:end-1).*diff(y1);
b2 = y2(1:end-1).*diff(x2) - x2(1:end-1).*diff(y2);

%pre-allocate outputs 
c = false;
idx1 = [];
idx2 = [];
for ii = 1:length(x1) - 1
    for jj = 1:length(x2) - 1
        xy = [a1(ii,:); a2(jj,:)] \ [b1(ii); b2(jj)];
        dxy1 = (xy(1) - x1(ii))/ (x1(ii+1) - x1(ii));
        dxy2 = (xy(1) - x2(jj))/ (x2(jj+1) - x2(jj));
        c = dxy1 >= 0 & dxy1 <= 1 & dxy2 >= 0 & dxy2 <= 1;
        
        if c;
            idx1 = ii;
            idx2 = jj;
            return; 
        end
    end
end
    