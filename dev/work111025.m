%load in a nailfold for experimenting with
n_files = dir([nailfoldroot 'images\anonymous_oct\bmp\*.bmp']);

nailfold = imread([nailfoldroot 'images\anonymous_oct\bmp\' n_files(4).name]);
nailfold = nailfold(:,:,1);
nailfold_mask = imerode(nailfold > 10 & nailfold < 250, strel('disk', 20));
n_min = min(nailfold(nailfold_mask));
n_max = max(nailfold(nailfold_mask));
%%
%Compute Gaussian first and second derivatives
[g1d_m, g1d_o] = ...
    gaussian_1st_derivative_gradient2(nailfold, [2 4]);
[g2d_m, g2d_o] = ...
    gaussian_clover_line(nailfold, [2 4]);

%Compute line mask of vessels from 2nd deriv.
line_nms = mb_non_maximal_supp(g2d_m, g2d_o);
[line_mask] = bwareaopen(line_nms > 0, 80);
[y_vessels x_vessels] = find(line_mask);
%%
%load in the apex templates
load([nailfoldroot 'data\apex_detection\apex_gd_templates.mat'], 'apex_template*');

%Apply normalised cross correlations using the 2 apex templates
C1 = normxcorr2(apex_template1, g1d_m);
C1 = C1(17:end-16, 17:end-16); %Need to discard floor(patch_size)/2 from start and end

C2 = normxcorr2(apex_template2, g2d_m);
C2 = C2(17:end-16, 17:end-16); %Need to discard floor(patch_size)/2 from start and end

%Look for local maxima in the 2 independent NCC maps
[maxima_pos1, maxima_vals1] = local_image_maxima(C1, 33, nailfold_mask, 0.5);
[maxima_pos2, maxima_vals2] = local_image_maxima(C2, 33, nailfold_mask, 0.5);

%Now multiply the 2 NCC maps to get a single map and again look for local
%maxima
C12 = C1 .* C2;
%C12(~nailfold_mask) = 0;
[maxima_pos12, maxima_vals12] = local_image_maxima(C12, 33, nailfold_mask, 0.2);

%%
figure; imagesc(nailfold); axis image; colormap(gray(256)); hold on; caxis([n_min n_max]);

nailfold_roi = nailfold(251:450, 1301:1400);
figure; imagesc(nailfold_roi); axis image; colormap(gray(256));

gd1_mr = abs(g1d_m(251:450, 1301:1400));
gd1_or = g1d_o(251:450, 1301:1400);
gd2_mr = g2d_m(251:450, 1301:1400);
gd2_or = g2d_o(251:450, 1301:1400);
%%
jetstream(gd1_mr, gd1_or, [60 6], [cos(gd1_or(6,60)) -sin(gd1_or(6,60))],...
    'sigma_theta', 0.1, 'N', 200, 'plot', 1);
jetstream(gd1_mr, gd1_or, [60 6], [-cos(gd1_or(6,60)) sin(gd1_or(6,60))],...
    'sigma_theta', 0.1, 'N', 200, 'plot', 1);

jetstream(gd2_mr, gd2_or, [60 13], [cos(gd2_or(13,60)) -sin(gd2_or(13,60))],...
    'lambda', mean(g2d_m((g2d_m > 0) & (g2d_m < 40))), 'sigma_theta', 0.1, 'plot', 1);

jetstream2(gd1_mr, gd1_or, [60 13], [cos(gd1_or(6,60)) -sin(gd1_or(6,60))], 7,...
    'sigma_theta', 0.1, 'sigma_m', 0.05, 'N', 200, 'plot', 1);
jetstream2(gd1_mr, gd1_or, [60 13], [-cos(gd1_or(6,60)) sin(gd1_or(6,60))], 7,...
    'sigma_theta', 0.1, 'sigma_m', 0.05, 'N', 200, 'plot', 1);
%%
particles1 = jetstream3(gd1_mr, gd2_mr, gd1_or, [60 13], [cos(gd1_or(6,60)) -sin(gd1_or(6,60))], [5 7],...
    'sigma_theta', 0.1, 'sigma_m', 0.2, 'N', 200, 'M', 500, 'plot', 1,...
    'lambda2', mean(g2d_m((g2d_m > 0) & (g2d_m < 40))));
%
particles2 = jetstream3(gd1_mr, gd2_mr, gd1_or, [60 13], [-cos(gd1_or(6,60)) sin(gd1_or(6,60))], [7 5],...
    'sigma_theta', 0.1, 'sigma_m', 0.2, 'N', 200, 'M', 500, 'plot', 1,...
    'lambda2', mean(g2d_m((g2d_m > 0) & (g2d_m < 40))));

figure; imagesc(nailfold_roi); axis image; colormap(gray(256)); hold on;
plot(mean(particles1.xp), mean(particles1.yp), 'r');
plot(mean(particles1.xm), mean(particles1.ym), 'g');
plot(mean(particles2.xp), mean(particles2.yp), 'g');
plot(mean(particles2.xm), mean(particles2.ym), 'r');
%%
profile on;
particles1 = jetstream3(gd1_mr, gd2_mr, gd1_or, [60 13], [cos(gd1_or(6,60)) -sin(gd1_or(6,60))], [5 7],...
    'sigma_theta', 0.1, 'sigma_m', 0.1, 'N', 200, 'plot', 0, 'use_nag', 0,...
    'lambda2', mean(g2d_m((g2d_m > 0) & (g2d_m < 40))));
profile viewer;
%%
%profile on;
figure; imagesc(nailfold); axis image; colormap(gray(256)); hold on; caxis([n_min n_max]);
%figure; imagesc(line_mask); axis image; colormap(gray(256)); hold on;
plot(maxima_pos12(:,1), maxima_pos12(:,2), 'b+', 'markersize', 2);
num_pts = length(maxima_vals12);
for ii = 1:num_pts
    x = maxima_pos12(ii,1);
    y = maxima_pos12(ii,2);
    neighbours = (abs(x_vessels - x) < 10) & (abs(y_vessels - y) < 10);
            
    if ~any(neighbours); continue; end
    
    xn = x_vessels(neighbours);
    yn = y_vessels(neighbours);

    dists = (xn - x).^2 + (yn - y).^2;
    [dummy min_idx] = min(dists);
    x = xn(min_idx);
    y = yn(min_idx);
    
    particles1 = jetstream3(g1d_m, g2d_m, g1d_o, [x y], [cos(g2d_o(y, x)) -sin(g2d_o(y,x))], [5 7],...
        'sigma_theta', 0.1, 'sigma_m', .2, 'N', 200, 'plot', 0, 'M', 500,...
        'lambda2', mean(g2d_m((g2d_m > 0) & (g2d_m < 40))));
    particles2 = jetstream3(g1d_m, g2d_m, g1d_o, [x y], [-cos(g2d_o(y, x)) sin(g2d_o(y,x))], [7 5],...
        'sigma_theta', 0.1, 'sigma_m', .2, 'N', 200, 'plot', 0, 'M', 500,...
        'lambda2', mean(g2d_m((g2d_m > 0) & (g2d_m < 40))));

    plot(mean(particles1.xp), mean(particles1.yp), 'r--');
    plot(mean(particles1.xm), mean(particles1.ym), 'g--');
    plot(mean(particles2.xp), mean(particles2.yp), 'g--');
    plot(mean(particles2.xm), mean(particles2.ym), 'r--');
    
end
%profile viewer;
print_pdf('vessel_track.pdf');
