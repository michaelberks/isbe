function [int_z] = tps_warp(L_inv, source_x, source_y, source_z, int_x, int_y, green)

% source_x, source_y are (x, y) co-ordinates of landmark points
% int_x, int_y are (x, y) co-ordinates of points to be interpolated
% return int_z the interpolated z-displacements
if nargin < 7
    green = 'biharmTPS';
end

source_len = length(source_x); % = length(source_y)
int_len = length(int_x); % = length(int_y)

V_z = [source_z, 0 0 0]';
W_z = L_inv * V_z; clear V_z L_inv;

%
%Code below avoids large loops
%
memlimit = 800000;  %seems to be OK for a 128Mb memory

segments = ceil(source_len*int_len/memlimit);
llx = round(linspace(0, int_len, segments+1));
int_xy = [int_x; int_y]; clear int_x int_y;

for j=1:segments 
    int_z(llx(j)+1:llx(j+1)) = ...
       W_z' * [spline_green([source_x; source_y], int_xy(:,llx(j)+1:llx(j+1)), green) ; ...
           ones(1,size(int_xy(:,llx(j)+1:llx(j+1)),2)); int_xy(:,llx(j)+1:llx(j+1))];
end
    