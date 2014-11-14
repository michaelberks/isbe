load('C:\isbe\nailfold\images\ncm\annotations\n Tonia MooreV1LD4X3LrgMosaic_vessels.mat');
nailfold = imread('C:\isbe\nailfold\images\ncm\n Tonia MooreV1LD4X3LrgMosaic.bmp');

sr = round(min(vessels{1}(:,2)) - 50);
er = round(max(vessels{1}(:,2)) + 50);
sc = round(min(vessels{1}(:,1)) - 50);
ec = round(max(vessels{1}(:,1)) + 50);

nailfold = nailfold(sr:er, sc:ec);
[mag_2d, ori_2d] = gaussian_clover_line(nailfold, 4);
[mag_1d ori_1d] = gaussian_1st_derivative_gradient(nailfold, 1);

%
dv = diff(vessels{1});
keep = [true; any(dv,2)];
vessel = [vessels{1}(keep,1) - sc vessels{1}(keep,2) - sr];


dists = cumsum([0; sqrt(sum(diff(vessel).^2, 2))]);
len = floor(dists(end));
vessel2 = interp1(dists, vessel, linspace(0, dists(end), len), 'linear');

figure; imagesc(nailfold); axis image; colormap(gray(256)); hold on;
plot(vessel(:,1), vessel(:,2));
plot(vessel2(:,1), vessel2(:,2), 'r.');

norm_width = 40;

% Pre-allocate containers for the normal profiles and the associated
% sampling points
normal_p = zeros(len, norm_width);
normal_x = zeros(len, norm_width);
normal_y = zeros(len, norm_width);
normal_g1 = zeros(len, norm_width);
normal_g2 = zeros(len, norm_width);
%Compute the normal vectors at each point on the inner border
[fx, fy] = gradient(vessel2);

%smooth fy
fy = imfilter(fy, ones(5,1), 'replicate');

%normalise fy
fy = fy ./ [sqrt(sum(fy.^2, 2)), sqrt(sum(fy.^2, 2))];

%Compute normal profiles of the image at every point
for ii = 1:len %= number of rows in skin_air
    
    n1_x = vessel2(ii, 1) - norm_width*fy(ii, 2)/2;
    n1_y = vessel2(ii, 2) + norm_width*fy(ii, 1)/2;
    n2_x = vessel2(ii, 1) + norm_width*fy(ii, 2)/2;
    n2_y = vessel2(ii, 2) - norm_width*fy(ii, 1)/2;
    
    plot([n1_x n2_x], [n1_y n2_y], 'y');
    plot(n1_x, n1_y, 'gx');
    plot(n2_x, n2_y, 'rx');

    [cx, cy, cp] = improfile(nailfold, [n1_x, n2_x], [n1_y, n2_y], norm_width);
    normal_p(ii, :) = cp';
    normal_x(ii, :) = cx';
    normal_y(ii, :) = cy';
    
    normal_g1(ii,:) = improfile(mag_1d, [n1_x, n2_x], [n1_y, n2_y], norm_width);
    normal_g2(ii,:) = improfile(mag_2d, [n1_x, n2_x], [n1_y, n2_y], norm_width);

end

normal_n = bsxfun(@rdivide, normal_p, sum(normal_p,2));
figure; 
subplot(1,4,1); imagesc(normal_p); axis image;
subplot(1,4,2); imagesc(normal_n); axis image;
subplot(1,4,3); imagesc(normal_g1); axis image;
subplot(1,4,4); imagesc(normal_g2); axis image;

%figure; surf(normal_n); axis image;
%%
inner_vessel_mask = roipoly(nailfold, vessel(:,1), vessel(:,2));
outer_vessel_mask = imdilate(inner_vessel_mask, strel('disk', 40)) & ~inner_vessel_mask;

figure; 
subplot(1,2,1); imagesc(inner_vessel_mask); axis image;
subplot(1,2,2); imagesc(outer_vessel_mask); axis image;
%
[mag_2d, ori_2d] = ...
    gaussian_clover_line(nailfold, 4);
nms_2d = bwareaopen(...
    mb_non_maximal_supp(mag_2d, ori_2d) > 0, 0);
[y_2d x_2d] = find(nms_2d);

[mag_1d ori_1d] = gaussian_1st_derivative_gradient2(nailfold, [1 2]);
nms_1d = bwareaopen(...
    mb_non_maximal_supp(mag_1d, ori_1d) > 0, 0);
[y_1di x_1di] = find(nms_1d);% & inner_vessel_mask);
[y_1do x_1do] = find(nms_1d);% & outer_vessel_mask);

figure;
subplot(1,3,1); imagesc(nailfold); axis image; colormap(gray(256)); hold on;
plot(x_1di, y_1di, 'g.', 'markersize', 2);
plot(x_1do, y_1do, 'r.', 'markersize', 2);
plot(x_2d, y_2d, 'y.', 'markersize', 2);
subplot(1,3,2); imagesc(mag_1d); axis image; colormap(gray(256)); hold on;
subplot(1,3,3); imagesc(mag_2d); axis image; colormap(gray(256)); hold on;
%%
figure;
subplot(1,2,1); imagesc(ori_1d); axis image; colormap(hsv(256)); hold on; colorbar; caxis([-pi pi]);
plot(x_2d, y_2d, 'k.', 'markersize', 2);
subplot(1,2,2); imagesc(ori_2d); axis image; colormap(hsv(256)); hold on; colorbar; caxis([-pi pi]);
plot(x_2d, y_2d, 'k.', 'markersize', 2);
%%
figure;
imagesc(nailfold); axis image; colormap(gray(256)); hold on;
for ii = 1:len %= number of rows in skin_air
    
    inner_dist = (x_1di - vessel2(ii, 1)).^2 + (y_1di - vessel2(ii, 2)).^2;
    outer_dist = (x_1do - vessel2(ii, 1)).^2 + (y_1do - vessel2(ii, 2)).^2;
    
    [dummy min_in] = min(inner_dist);
    [dummy min_out] = min(outer_dist);
    nxi = x_1di(min_in);
    nyi = y_1di(min_in);
    nxo = x_1do(min_out);
    nyo = y_1do(min_out);
    
    plot([nxi vessel2(ii, 1) nxo], [nyi vessel2(ii, 2) nyo], 'y');
    plot(nxi, nyi, 'gx');
    plot(nxo, nyo, 'rx');

end
%%
%Compute the normal vectors at each point on the inner border
[fx, fy] = gradient(vessel2);

%normalise fy
fy = fy ./ [sqrt(sum(fy.^2, 2)), sqrt(sum(fy.^2, 2))];

figure;
imagesc(nailfold); axis image; colormap(gray(256)); hold on;

for ii = 1:len
    
    v_to_no = [-fy(ii,2) fy(ii,1)];
    v_to_ni = [fy(ii,2) -fy(ii,1)];
    
    v_to_eo = bsxfun(@minus, [x_1do y_1do], vessel2(ii,:));
    v_to_eo = bsxfun(@rdivide, v_to_eo, sqrt(sum(v_to_eo.^2,2)));
    
    v_to_ei = bsxfun(@minus, [x_1di y_1di], vessel2(ii,:));
    v_to_ei = bsxfun(@rdivide, v_to_ei, sqrt(sum(v_to_ei.^2,2)));

    [dummy min_out] = min(sum(bsxfun(@minus, v_to_eo, v_to_no).^2,2));
    [dummy min_in] = min(sum(bsxfun(@minus, v_to_ei, v_to_ni).^2,2));
    
    nxi = x_1di(min_in);
    nyi = y_1di(min_in);
    nxo = x_1do(min_out);
    nyo = y_1do(min_out);
    
    %plot([nxi vessel2(ii, 1) nxo], [nyi vessel2(ii, 2) nyo], 'y');
    plot([nxi vessel2(ii, 1)], [nyi vessel2(ii, 2)], 'y');
    plot(nxi, nyi, 'rx');
    %plot(nxo, nyo, 'gx');

end
%%
T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
            'transform', 'spline');
        [pts] = geom_transformpoints([i_x; i_y], T);
        
%%
p = 1+7i;
q = 1-3i;
len = 100;
pq = linspace(p, q, len);

n = i*(q-p)/abs(q-p);

npq1 = pq + n;
npq2 = pq - n;

pqt = pq.^2;
npq1t = npq1.^2;
npq2t = npq2.^2;

figure; 
subplot(1,2,1); hold all; axis equal;
plot(real(pq), imag(pq));
plot(real([npq1; npq2]), imag([npq1; npq2]));

subplot(1,2,2); hold all; axis equal;
plot(real(pqt), imag(pqt));

for ii = 1:len
    nn = linspace(npq1(ii), npq2(ii), len);
    plot(real(nn.^2), imag(nn.^2));
end

        
