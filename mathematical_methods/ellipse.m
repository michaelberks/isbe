function [x y] = ellipse(r_x, r_y, x0 , y0, x_axis, g_noise, n_pts)

if nargin < 7
    n_pts = 100;
end
if nargin < 6
    g_noise = 0;
end
if nargin < 5
    x_axis = [1 0];
end
if nargin < 3
    x0 = 0;
    y0 = 0;
end
scale = sum(x_axis.^2);
if scale ~= 1
    x_axis = x_axis / scale;
end

co = x_axis(1);
si = x_axis(2);
pts = linspace(0, 2*pi, n_pts);

x_noise = zeros(1, n_pts);
y_noise = zeros(1, n_pts);

if g_noise
    x_noise = g_noise*randn(1, n_pts);
    y_noise = g_noise*randn(1, n_pts);
end

x = (r_x*cos(pts)+x_noise)*co - (si*r_y*sin(pts)+x_noise) + x0;
y = (r_x*cos(pts)+y_noise)*si + (co*r_y*sin(pts)+y_noise) + y0;


