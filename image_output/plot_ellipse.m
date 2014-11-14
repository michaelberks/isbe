function plot_ellipse(r_x, r_y, x_axis, x0 , y0, linestyle)

co = x_axis(1);
si = x_axis(2);
pts = linspace(0, 2*pi, 100);
x = r_x*cos(pts)*co - si*r_y*sin(pts) + x0;
y = r_x*cos(pts)*si + co*r_y*sin(pts) + y0;

plot(x, y, linestyle);
% plot(x, y, 'rx');