
borders = zeros(185, 1);
for ii = 1:185
    load(['C:\isbe\dev\annotations\', files(ii).name]);
    borders(ii) = size(mass.mass_outline, 1);
end
[sb inds] = sort(borders);
%%
figure('WindowStyle', 'docked');

for ii = 1:9
    load(['C:\isbe\dev\annotations\', files(inds(ii)).name]);
    subplot(3,3,ii);
    plot(mass.mass_outline([1:end,1], 1), mass.mass_outline([1:end,1], 2));
    hold on;
    [r_x, r_y, x0 , y0, x_axis] = fit_ellipse(mass.mass_outline);
    [pts_x pts_y] = ellipse(r_x, r_y, x0 , y0, x_axis);
    plot(pts_x, pts_y, 'r');
end

%%

%make me radiate!!

%start_vec = [0, 1];
dist = 10;
figure; hold on;
for ii = 0:7
    theta = 2*pi*ii / 8;
    start_vec = [cos(theta) sin(theta)];
    plot(dist*[0 start_vec(1)], dist*[0 start_vec(2)]);
end
%%
[pts_x pts_y] = ellipse(5, 10, 0 , 0);
plot(pts_x, pts_y, 'r');
%%
%loci_im = zeros(200);
%[y_coord x_coord] = ind2sub(size(loci_im), 1:length(loci_im(:)));
xpts = 1:200;
for ii = 1:8
    loci_im = zeros(200);
    theta1 = 2*ii*pi/8;
    theta2 = 2*(ii+1)*pi/8;
    c1 = 100 - 100*tan(theta1);
    c2 = 100 - 100*tan(theta2);
    ypts1 = tan(theta1)*xpts + c1;
    ypts2 = tan(theta2)*xpts + c2;
    inds = (y_coord > (tan(theta2)*x_coord + c2)) & ...
           (y_coord <= (tan(theta1)*x_coord + c1));
    loci_im(inds) = 1;
    figure('WindowStyle', 'docked'); imagesc(loci_im); colormap(jet(256)); 
    axis image; hold on;
    plot(xpts, ypts1, 'g');
    plot(xpts, ypts2, 'r');
end
figure; imagesc(loci_im); axis image;
