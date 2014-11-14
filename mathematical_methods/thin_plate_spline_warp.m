function [int_z] = thin_plate_spline_warp(source_x, source_y, source_z, int_x, int_y, green)

% source_x, source_y are (x, y) co-ordinates of landmark points, source_z
% the associated z-displacement
% int_x, int_y are (x, y) co-ordinates of points to be interpolated
% return int_z the interpolated z-displacements

if nargin < 6
    green = 'biharmTPS';
end

int_len = length(int_x); % = length(int_y)

source_len = length(source_x); % = length(source_y)

K = spline_green([source_x; source_y],[], green);

P = [repmat(1, source_len, 1), source_x', source_y'];

L = [[K, P]; [P', zeros(3,3)]]; clear  K P;

V_z = [source_z, 0 0 0]'; clear source_z;

display(['rcond = ', num2str(rcond(L))]);

%try using left divide to calculate weights - most computationally efficient
if rcond(L) > 1e-8,
    W_z = L \ V_z;
else
%if \ doesn't work, use moore-penrose pseudo-inverse (based on SVD)
    display('Matrix L is close to singular, pseudo-inverse used');
    W_z = pinv(L) * V_z;
end
clear L V_z;

memlimit = 800000;  %seems to be OK for a 128Mb memory
segments = ceil(source_len*int_len/memlimit);
llx = round(linspace(0, int_len, segments+1));
%src_xy = [source_x; source_y];
int_xy = [int_x; int_y]; clear int_x int_y;

for j=1:segments 

    int_z(llx(j)+1:llx(j+1)) = ...
       W_z' * [spline_green([source_x; source_y], int_xy(:,llx(j)+1:llx(j+1)), green) ; ...
           ones(1,size(int_xy(:,llx(j)+1:llx(j+1)),2)); int_xy(:,llx(j)+1:llx(j+1))];
end