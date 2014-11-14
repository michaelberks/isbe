function [L_inv] = tps_weights(source_x, source_y, green)

% source_x, source_y are (x, y) co-ordinates of landmark points, source_z
% the associated z-displacement
% int_x, int_y are (x, y) co-ordinates of points to be interpolated
% return int_z the interpolated z-displacements

if nargin < 3
    green = 'biharmTPS';
end

source_len = length(source_x); % = length(source_y)

K = spline_green([source_x; source_y],[], green);

P = [repmat(1, source_len, 1), source_x', source_y'];

L = [[K, P]; [P', zeros(3,3)]]; clear  K P;

% use moore-penrose pseudo-inverse (based on SVD)
L_inv = pinv(L);

clear L;


%
%
%K = zeros(source_len, source_len);
%x_v = repmat(source_x', 1 ,source_len);
%
%y_v = repmat(source_y', 1 ,source_len);
%
%rsq =(x_v-x_v').^2 + (y_v-y_v').^2;
%
%clear x_v y_v;
%
% In order to avoid log(0) errors, all zeroes have been replaced by ones
% This sets the leading diagonal to zero
%
%rsq(rsq < 1e-8) = 1;
%
%K = rsq.*log(rsq);
%
%clear rsq;
