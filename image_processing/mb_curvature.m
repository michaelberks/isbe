function [curvature] = mb_curvature(line_in,d)
%MB_CURVATURE *Insert a one line summary here*
%   [curvature] = mb_curvature(line_in,d)
%
% Inputs:
%      line_in- *Insert description of input variable here*
%
%      d- *Insert description of input variable here*
%
%
% Outputs:
%      curvature- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 16-Nov-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
n_pts = size(line_in,1);
curvature = zeros(n_pts - 2*d, 1);

figure; 
plot(line_in(:,1), line_in(:,2)); hold on; axis equal;

for ii = d+1:n_pts-d
    
    A = line_in(ii,:);
    B = line_in(ii-d,:);
    C = line_in(ii+d,:);
    
    D = C + ((B-C)*(A-C)')*(B-C) / sum((B-C).^2);
%     sgn = sign((C-A)*(B-A)');
    curvature(ii-d) = sum((A-D).^2);
%     
%     if sgn > 0
%         plot([B(1) D(1)], [B(2) D(2)], 'r');
%     else
%         plot([B(1) D(1)], [B(2) D(2)], 'g');
%     end
    
    ABdotAC = ((B-A)*(C-A)') / (sqrt(sum((B-A).^2))*sqrt(sum((C-A).^2)));
%     curvature(ii-d) = acos(ABdotAC);
    if ABdotAC > 0
        plot([A(1) D(1)], [A(2) D(2)], 'r');
    else
        plot([A(1) D(1)], [A(2) D(2)], 'g');
        plot([B(1) C(1)], [B(2) C(2)], 'r');
    end
end
