function h=plot_warp(CPS)

source_x = [0, 0, -.10, .10];
source_y = [-.10, .10, 0, 0];
source_z = [.005, .005, -.005, -.005];
source_len = 4;%length(source_x);

[xx, yy] = meshgrid(-1:.01:1, -1:.01:1);
int_x = xx(:); int_y = yy(:);

if CPS
    int_z = clamped_plate_spline_warp(source_x, source_y, source_z, int_x, int_y);
else
    L_inv = tps_weights(source_x, source_y);
    int_z = tps_warp(L_inv, source_x, source_y, source_z, int_x, int_y);
end

figure, surf(xx, yy, reshape(int_z, size(xx)));
hold on;
plot3(source_x, source_y, source_z, 'bo');
axis([-1, 1, -1, 1, -.2, .2]);
hold on;