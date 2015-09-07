oh = load('C:\isbe\nailfold\models\vessel\ori_histo.txt');
hog_patch = load('C:\isbe\nailfold\models\vessel\hog_patch.txt');
ori_histo = zeros(8,8,12);

for i_ori = 1:8
    rows = 8*(i_ori-1) + (1:8);
    ori_histo(:,:,i_ori) = oh(rows,:);
end

norm = sum(ori_histo(:)) + eps;
hog = sqrt(ori_histo / norm);

figure; imgray(hog_patch);
angles = linspace(0, pi, 13);
colors = lines(12);

for i_x = 1:8
    for i_y = 1:8
        for i_ori = 1:12
            theta = angles(i_ori);
            x = 8*i_x - 3.5;
            y = 8*i_y - 3.5;
            u = 16*hog(i_y,i_x,i_ori)*[-cos(theta) cos(theta)];
            v = 16*hog(i_y,i_x,i_ori)*[sin(theta) -sin(theta)];
            plot(x + u, y + v, 'color', colors(i_ori,:));
        end
    end
end
%%
nx = 8; ny = 8; na = 9;
angles = linspace(-pi/2, pi/2, na+1);
colors = lines(na);
[hog_mat] = compute_HoG(hog_patch, 'num_ori_bins', na);
hog = hog_mat(:);

figure; imgray(hog_patch);
for i_x = 1:8
    for i_y = 1:8
        for i_ori = 1:9
            i_y0 = i_y-1;
            i_x0 = i_x-1;
            i_o0 = i_ori-1;
            col_c = ny*((i_ori-1)*nx + (i_x-1)) + i_y; %na*(i_y0*nx + i_x0) + i_o0 + 1;
            
            theta = angles(i_ori);
            x = 8*i_x - 3.5;
            y = 8*i_y - 3.5;
            u = 16*hog(col_c)*[-cos(theta) cos(theta)];
            v = 16*hog(col_c)*[-sin(theta) sin(theta)];
            plot(x + u, y + v, 'color', colors(i_ori,:));
        end
    end
end
%%
nx = 8; ny = 8; na = 9;
angles = linspace(-pi/2, pi/2, na+1);
colors = lines(na);
hog = test_data(81,1:576);
%norm = sum(hog(:)) + eps;
%hog = sqrt(hog / norm);

figure; imgray(hog_patch);
for i_x = 1:8
    for i_y = 1:8
        for i_ori = 1:9
            i_y0 = i_y-1;
            i_x0 = i_x-1;
            i_o0 = i_ori-1;
            col_c = ny*((i_ori-1)*nx + (i_x-1)) + i_y; %na*(i_y0*nx + i_x0) + i_o0 + 1;
            
            theta = angles(i_ori);
            x = 8*i_x - 3.5;
            y = 8*i_y - 3.5;
            u = 16*hog(col_c)*[-cos(theta) cos(theta)];
            v = 16*hog(col_c)*[-sin(theta) sin(theta)];
            plot(x + u, y + v, 'color', colors(i_ori,:));
        end
    end
end
%%
