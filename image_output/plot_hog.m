function plot_hog(hog_patch, hog, nx, angle_wrap)

[nxy na] = size(hog);
ny = nxy / nx;

if angle_wrap
    angles = linspace(-pi/2, pi/2, na+1);
else
    angles = linspace(-pi, pi, na+1);
end

colors = lines(na);
cell_sz_x = (size(hog_patch,2)-2)/nx;
cell_sz_y = (size(hog_patch,1)-2)/ny;

figure; imgray(hog_patch);
for i_x = 1:nx
    for i_y = 1:ny
        for i_ori = 1:na
            hr = ny*(i_x-1) + i_y;
            theta = angles(i_ori);
            x = cell_sz_x*i_x - (cell_sz_x-1)/2;
            y = cell_sz_y*i_y - (cell_sz_y-1)/2;
            
            if angle_wrap
                u = 2*cell_sz_x*hog(hr,i_ori)*[-cos(theta) cos(theta)];
                v = 2*cell_sz_y*hog(hr,i_ori)*[-sin(theta) sin(theta)];
            else
                u = 4*cell_sz_x*hog(hr,i_ori)*[0 cos(theta)];
                v = 4*cell_sz_y*hog(hr,i_ori)*[0 sin(theta)];
            end
            plot(x + u, y + v, 'color', colors(i_ori,:));
        end
    end
end