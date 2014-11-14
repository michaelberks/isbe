function [auc] = area_under_curve(x, y, xt1, xt2)
%AREA_UNDER_CURVE *Insert a one line summary here*
%   [auc] = area_under_curve(x, y, t1, t2, method)
%
% Inputs:
%      x - *Insert description of input variable here*
%
%      y - *Insert description of input variable here*
%
%      xt1 - *Insert description of input variable here*
%
%      xt2 - *Insert description of input variable here*
%
%      method - *Insert description of input variable here*
%
%
% Outputs:
%      auc - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 29-Nov-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
A = sortrows([x(:) y(:)]);
x = A(:,1);
y = A(:, 2);

i1 = find(x <= xt1, 1, 'last');
j1 = find(x > xt1, 1);
if isempty(i1) || isempty(j1)
    display(['warning: 1st threshold ' num2str(xt1) ' is out of range']);
else
    yt1 = interp1(x([i1 j1]), y([i1 j1]), xt1, 'linear');
    
    x = [xt1; x(j1:end)];
    y = [yt1; y(j1:end)];
end

i2 = find(x <= xt2, 1, 'last');
j2 = find(x > xt2, 1);
if isempty(i2) || isempty(j2)
    display(['warning: 2nd threshold ' num2str(xt2) ' is out of range']);
else
    yt2 = interp1(x([i2 j2]), y([i2 j2]), xt2, 'linear');
    
    x = [x(1:i2); xt2];
    y = [y(1:i2); yt2];
end

auc = trapz(x, y);
