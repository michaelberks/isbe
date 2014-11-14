random_forest = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\305894\random_forest.mat');
random_forest.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
save('C:\isbe\asymmetry_project\data\line_orientation_rfs\305894\random_forest.mat', 'random_forest');

random_forest = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\306060\random_forest.mat');
random_forest.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
save('C:\isbe\asymmetry_project\data\line_orientation_rfs\306060\random_forest.mat', 'random_forest');

random_forest = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\306079\random_forest.mat');
random_forest.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
save('C:\isbe\asymmetry_project\data\line_orientation_rfs\306079\random_forest.mat', 'random_forest');
%%
roi = u_load('C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\mass_roi\024RCC_roi.mat');

rf_3_4 = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\305894\random_forest.mat');
rf_1_4 = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\306060\random_forest.mat');
rf_1_5 = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\306079\random_forest.mat');

args_3_4 = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\305894\sampling_args.mat');
args_1_4 = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\306060\sampling_args.mat');
args_1_5 = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\306079\sampling_args.mat');

args.image_in = roi;
args.forest_type = 'orientation';

args.sampling_args.do_max = 0;
args.sampling_args.rotate = 0;
args.sampling_args.use_nag = 0;
args.sampling_args.feature_shape = 'rect';
args.sampling_args.feature_type = 'conj';

%
args.forest = rf_3_4;
args.sampling_args.num_levels = 4;
args.sampling_args.win_size = 3;
ori_map_3_4 = classify_image(args);
%
args.forest = rf_1_4;
args.sampling_args.num_levels = 4;
args.sampling_args.win_size = 1;
ori_map_1_4 = classify_image(args);
%
args.forest = rf_1_5;
args.sampling_args.num_levels = 5;
args.sampling_args.win_size = 1;
ori_map_1_5 = classify_image(args);
%%
%Extract pectoral region as the box bounded by the start of the breast
%air border and the centorid of the breast
[rows cols] = size(pectoral);

%Filter the pectoral region and discard high/low pixels
pectoral_med = medfilt2(pectoral, [9 9]);
pectoral_mask = pectoral_med > 240 | pectoral_med < 2;
pectoral_med(pectoral_mask) = 0;

%Perform canny edge detection on pectoral region
[pectoral_canny ori] = canny_edge(pectoral_med, [], 1, 0.98, .6);

%Discard edges that aren't approximatelt aligned on an upperleft
%diagonal
ignore_ori = ori > pi/4 | ori < pi/24;
pectoral_ori = pectoral_canny;
pectoral_ori(ignore_ori) = 0;

%Take hough transformation of edges to locate most likely staright line
[line_scores rhos thetas] = hough_line(pectoral_ori);
[max_rhos rho_idx] = max(line_scores);
[dummy theta_idx] = max(max_rhos);

%Convert from rho,theta parametisation to y = mx + c;
theta = thetas(theta_idx);
rho = rhos(rho_idx(theta_idx));
m = tan(theta); 
c = rho / cos(theta);

%Work out points on pectoral edge at top and bottom of image
y1 = 1; x1 = (y1 - c)/m;
y2 = rows; x2 = (y2 - c)/m;
pectoral_edge = [x1 y1; x2 y2];
%%
figure; imagesc(pectoral); axis image; colormap(gray(256));
hold on;
plot([x1 x2], [y1 y2]);

%Get normal points
x_norm = -sin(theta);
y_norm = cos(theta);

len = 100;
norm_width = 100;

normal_p = zeros(len, norm_width);
normal_x = zeros(len, norm_width);
normal_y = zeros(len, norm_width);

for ii = 1:len
    xx = x1 + ii*(x2-x1)/len;
    yy = y1 + ii*(y2-y1)/len;
    
    xn1 = xx - norm_width*x_norm/2;
    xn2 = xx + norm_width*x_norm/2;
    
    yn1 = yy - norm_width*y_norm/2;
    yn2 = yy + norm_width*y_norm/2;
    
    [cx, cy, cp] = improfile(pectoral, [xn1, xn2], [yn1, yn2], norm_width);
    normal_p(ii, :) = cp';
    normal_x(ii, :) = cx';
    normal_y(ii, :) = cy';
    
    plot(xx, yy, 'bx');
    plot(xn1, yn1, 'rx');
    plot(xn2, yn2, 'gx');
end
%%
figure; imagesc(pectoral_canny); axis image; hold on;
[y_edges x_edges] = find(pectoral_canny);
edge_dists = abs((x_edges - x1)*x_norm + (y_edges - y1)*y_norm);
edge_oris = pi/2 - ori(pectoral_canny);

valid_dists = edge_dists < 20;
valid_ori = abs(mb_mod(edge_oris + theta, pi)) < pi/12;
valid_edges = valid_dists & valid_ori;

plot(x_edges(valid_edges), y_edges(valid_edges), 'm.');
%%
r_max = 90;
r_min = 10;
R = 5;
num_angles = 12;
phi = atan(R ./ d);
%%
    [yy xx] = find(pectoral_ori);

    figure; imagesc(pectoral_ori); axis image; hold on;
    N = length(xx);
    x_plot = 1:cols;
    max_score = 0;
    ii = 0;
    while ii < 100
        r_idx = randperm(N);
        x_fit = xx(r_idx(1:3));
        if all(diff(sort(x_fit)))
            ii = ii+1;
            y_fit = yy(r_idx(1:3));
%             polys(ii,:) = polyfit(x_fit,y_fit,2);
%             y_plot = polyval(polys(ii,:), x_plot);

            y_plot = interp1(x_fit, y_fit, x_plot, 'cubic');
            plot(x_plot, y_plot, 'y');

            curr_score = 0;
            for jj = 1:N
                r2 = (xx(jj) - x_plot).^2 + (yy(jj) - y_plot).^2;
                curr_score = curr_score + (min(r2) < 10);
            end
            if curr_score > max_score;
                x_max = x_plot;
                y_max = y_plot;
                max_score = curr_score;
            end
        end
    end
    plot(x_max, y_max, 'r');
    
    y_diff = diff(y_max);
    idx = find(y_diff >= 0, 1);
    x_max(1:idx-1) = [];
    x_max(1:idx-1) = [];
    
    y_diff = diff(y_max);
    idx = find(y_diff < 0, 1);
    x_max(idx:end) = [];
    y_max(idx:end) = [];
    
    plot(x_max, y_max, 'g');

         
%%
ang_res = 2*pi/num_angles;
for ii = 1:num_angles;
    theta_min = ang_res*(ii-1);
    theta_max = ang_res*ii;

    bin_filt = ...
        (d > r_min) & (d < r_max) & (theta > theta_min) & (theta <= theta_max);
    
    figure; imagesc(bin_filt); axis image;
end
