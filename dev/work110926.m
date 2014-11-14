v1 = u_load('C:\isbe\nailfold\images\ncm\annotations\Andrea Murray_vessels.mat');
v2 = u_load('C:\isbe\nailfold\images\kK_demo\nailfolds\n3_mb_vessels.mat');

v1([1 3]) = [];
%%
min_x1 = min(min(v1{1}(:,1)), min(v1{2}(:,1)));
min_x2 = min(min(v2{1}(:,1)), min(v2{2}(:,1)));

min_y1 = min(min(v1{1}(:,2)), min(v1{2}(:,2)));
min_y2 = min(min(v2{1}(:,2)), min(v2{2}(:,2)));

figure; hold on;
for ii = 1:2
    plot(v1{ii}(:,1)-min_x1, v1{ii}(:,2)-min_y1);
    plot(v2{ii}(:,1)-min_x2, v2{ii}(:,2)-min_y2, 'r');
end
%%
len = 100;

vessels1 = zeros(2*len,2);
vessels2 = zeros(2*len,2);

for ii = 1:2
    keep = [true; any(diff(v1{ii}),2)];
    vessel = [v1{ii}(keep,1) v1{ii}(keep,2)];


    dists = cumsum([0; sqrt(sum(diff(vessel).^2, 2))]);
    vessels1((ii-1)*len+1:ii*len,:) = interp1(dists, vessel, linspace(0, dists(end), len), 'linear');
    
    keep = [true; any(diff(v2{ii}),2)];
    vessel = [v2{ii}(keep,1) v2{ii}(keep,2)];

    dists = cumsum([0; sqrt(sum(diff(vessel).^2, 2))]);
    vessels2((ii-1)*len+1:ii*len,:) = interp1(dists, vessel, linspace(0, dists(end), len), 'linear');
end
%%
[dist, z, trans] = mb_procrustes(vessels1, vessels2);

figure; hold on;
plot(vessels1(:,1),vessels1(:,2), 'rx');
plot(z(:,1),z(:,2), 'gx');
%%
n1 = imread('C:\isbe\nailfold\images\ncm\andrea Murray.bmp');
n2 = imread('C:\isbe\nailfold\images\kk_demo\nailfolds\n3_mb.bmp');
n1 = double(n1(:,:,1));
%%
[line_strength, line_orientation] = ...
    gaussian_clover_line(n1, [2 4]);%gaussian_2nd_derivative_line2

%Apply non-maximal suppression to skeletonise the line strength map
line_nms = mb_non_maximal_supp(line_strength, line_orientation);

%Apply hysterisis to select suitable lines from the NMS map
[line_mask] = hysterisis(line_nms, [], 0.98);
% line_mask = bwareaopen(line_nms > 0, 20);

%Extract (x,y) coordinates of the remaining lines
[y_tgt x_tgt] = find(line_mask);

[line_strength, line_orientation] = ...
    gaussian_clover_line(n2, [2 4]);%gaussian_2nd_derivative_line2

%Apply non-maximal suppression to skeletonise the line strength map
line_nms = mb_non_maximal_supp(line_strength, line_orientation);

%Apply hysterisis to select suitable lines from the NMS map
[line_mask] = hysterisis(line_nms, [], 0.98);
%line_mask = bwareaopen(line_nms > 0, 20);

%Extract (x,y) coordinates of the remaining lines
[y_src x_src] = find(line_mask);

xyt = trans.b * [x_src y_src] * trans.T + repmat(trans.c(1,:), length(x_src), 1);

figure; hold;
plot(x_tgt, y_tgt, 'r.');
plot(xyt(:,1), xyt(:,2), 'b.'); axis ij equal;
plot(vessels1(:,1),vessels1(:,2), 'g.');
plot(z(:,1),z(:,2), 'y.');
%%
v1a = u_load('C:\isbe\nailfold\images\ncm\annotations\Andrea Murray_vessels_all.mat');
v2a = u_load('C:\isbe\nailfold\images\kK_demo\nailfolds\n3_mb_vessels_all.mat');
figure; hold on; axis ij equal;
for ii = 1:length(v1a)
    if size(v1a{ii},1) > 1
        
        keep = [true; any(diff(v1a{ii}),2)];
        vessel = [v1a{ii}(keep,1) v1a{ii}(keep,2)];
        dists = cumsum([0; sqrt(sum(diff(vessel).^2, 2))]);
        vessel = interp1(dists, vessel, linspace(0, dists(end), floor(dists(end))), 'linear');
        
        plot(vessel(:,1), vessel(:,2), 'r.');
    end       
end
for ii = 1:length(v2a)
    if size(v2a{ii},1) > 1
        
        keep = [true; any(diff(v2a{ii}),2)];
        vessel = [v2a{ii}(keep,1) v2a{ii}(keep,2)];
        dists = cumsum([0; sqrt(sum(diff(vessel).^2, 2))]);
        vessel = interp1(dists, vessel, linspace(0, dists(end), floor(dists(end))), 'linear');
        
        vessel =  trans.b * vessel * trans.T + repmat(trans.c(1,:), size(vessel,1), 1);
        plot(vessel(:,1), vessel(:,2), 'b.');
    end
end
%%
%--------------------------------------------------------------------------
%%

% Test if we can locate apexs - first given annotations, looking for the
% apex given the vessel shape. Then we'll try template matching new vessels
% given their appearance

%load in the nailfold
nailfold = imresize(imread('C:\isbe\nailfold\images\kk_demo\nailfolds\n3_mb.bmp'), 0.75, 'bilinear');

%Compute it's G2D transformation
%[mag_2d, ori_2d] = gaussian_clover_line(nailfold, 2);
[mag_2d, ori_2d] = gaussian_clover_line(nailfold, 4);
[mag_1d] = gaussian_1st_derivative_gradient2(nailfold, 2);

%%Convert orientation to a complex form
%ori_2d = complex(cos(ori_2d), sin(ori_2d));

%Load in the annotated vessels
v = u_load('C:\isbe\nailfold\images\kK_demo\nailfolds\n3_mb_vessels_all.mat');
vessels = [];
for ii = 1:length(v)
    if size(v{ii},1) > 1;       
        %Discard duplicate points
        keep = [true; any(diff(v{ii}),2)];
        vessel = 0.75*[v{ii}(keep,1) v{ii}(keep,2)];

        %Sample points evenly along the vessel
        dists = cumsum([0; sqrt(sum(diff(vessel).^2, 2))]);
        vessels{end+1,1} = interp1(dists, vessel, linspace(0, dists(end), floor(dists(end))), 'linear'); %#ok
    end
end
num_vessels = length(vessels);

%Create a structure to hold vessel orientations in
vessel_apex = zeros(num_vessels, 2);
vessel_top = zeros(num_vessels, 2);
%%
figure; hold on; axis equal;

%Loop throught he vessels, getting the vessel orientations and looking for
%the maximum orientation change
for ii = 1:num_vessels
    
    xx = vessels{ii}(:,1);
    yy = vessels{ii}(:,2);
    
    xx = [xx(2); xx(1); xx; xx(end); xx(end-1);]; %#ok
    xt = xx(2:end-3) - xx(4:end-1);
    xtt = xx(1:end-4) - xx(5:end);
    yy = [yy(2); yy(1); yy; yy(end); yy(end-1);]; %#ok
    yt = yy(2:end-3) - yy(4:end-1);
    ytt = yy(1:end-4) - yy(5:end);

    xx([1:2, end-1:end],:) = [];
    yy([1:2, end-1:end],:) = [];
 
    k = (xt.* ytt + yt.*xtt) ./ ((xt.^2 + yt.^2).^(3/2));

    h = fspecial('gaussian', [30 1], 5);
    k_smooth = conv2(k,h, 'same');

    xyn = [-yt xt] ./ [sqrt(xt.^2 + yt.^2) sqrt(xt.^2 + yt.^2)];
    xd = [0; diff(xx)];
    yd = [0; diff(yy)];
    D = sum([xd yd] .* xyn, 2);
    D_smooth = conv2(D,h, 'same');
    s_idx = floor(length(D_smooth)/3);
    e_idx = floor(2*length(D_smooth)/3);
    
    [dummy max_idx] = max(abs(D_smooth(s_idx:e_idx)));
    max_idx = max_idx + s_idx - 1;
    
    plot(vessels{ii}(:,1), vessels{ii}(:,2)); 
    
    %text(xx(min_idx), yy(min_idx), num2str(ii));
    
    vessel_apex(ii,:) = [xx(max_idx) yy(max_idx)];
    
    [y_min min_idx] = min(vessels{ii}(:,2));
    vessel_top(ii,:) = [vessels{ii}(min_idx,1) y_min];
    plot(vessel_apex(ii,1), vessel_apex(ii,2), 'rx');
    plot(vessel_top(ii,1), vessel_top(ii,2), 'gx');
          
end
%%
for ii = [1 2 4 5 21]
    
    %Now do vessel curvature
    xx = vessels{ii}(:,1);
    yy = vessels{ii}(:,2);
    
    xx = [xx(2); xx(1); xx; xx(end); xx(end-1);]; %#ok
    xt = xx(2:end-3) - xx(4:end-1);
    xtt = xx(1:end-4) - xx(5:end);
    yy = [yy(2); yy(1); yy; yy(end); yy(end-1);]; %#ok
    yt = yy(2:end-3) - yy(4:end-1);
    ytt = yy(1:end-4) - yy(5:end);

    xx([1:2, end-1:end],:) = [];
    yy([1:2, end-1:end],:) = [];
 
    k = (xt.* ytt + yt.*xtt) ./ ((xt.^2 + yt.^2).^(3/2));

    h = fspecial('gaussian', [30 1], 5);
    k_smooth = conv2(k,h, 'same');

    xyn = [-yt xt] ./ [sqrt(xt.^2 + yt.^2) sqrt(xt.^2 + yt.^2)];
    xd = [0; diff(xx)];
    yd = [0; diff(yy)];
    D = sum([xd yd] .* xyn, 2);
    D_smooth = conv2(D,h, 'same');
    [dummy min_idx] = max(abs(D_smooth));
    
    figure;
    subplot(1,2,1); hold on;
    plot(vessels{ii}(:,1), vessels{ii}(:,2));
    plot(vessels{ii}(1,1), vessels{ii}(1,2), 'gx');
    plot(xx(min_idx), yy(min_idx), 'rx');
    subplot(1,2,2); plot(1:length(D_smooth), D_smooth);
          
end
%%
%Ok so we've got some apexes (api? apices?) - lets sample from the
%centreline magnitude image
figure; a1 = gca; imagesc(nailfold); axis image; colormap(gray(256)); hold on;
f1 = figure;
for ii = 1:12
    
    %Check we're dealing with a vaild vessel
    plot(a1, vessels{ii}(:,1), vessels{ii}(:,2));
    plot(a1, vessel_top(ii,1), vessel_top(ii,2), 'rx');
    
    mag_patch = sample_window(mag_2d, 33, round(vessel_top(ii,2)), round(vessel_top(ii,1)));
    mag_patch = (mag_patch - mean(mag_patch(:))) / std(mag_patch(:));
    
    figure(f1);
    subplot(3,4,ii); imagesc(mag_patch); axis image; colormap(gray(256));
end
%%
mag_patches1 = zeros(33,33,num_vessels);
mag_patches2 = zeros(33,33,num_vessels);
for ii = 1:num_vessels
    mag_patch = sample_window(mag_1d, 33, round(vessel_top(ii,2)), round(vessel_top(ii,1)));
    mag_patch = (mag_patch - mean(mag_patch(:))) / std(mag_patch(:));
    mag_patches1(:,:,ii) = mag_patch;
    
    mag_patch = sample_window(mag_2d, 33, round(vessel_top(ii,2)), round(vessel_top(ii,1)));
    mag_patch = (mag_patch - mean(mag_patch(:))) / std(mag_patch(:));
    mag_patches2(:,:,ii) = mag_patch;
end

figure; imagesc(mean(mag_patches1,3)); axis image;
figure; imagesc(mean(mag_patches2,3)); axis image;
%%
nailfold = imread('C:\isbe\nailfold\images\ncm\andrea Murray.bmp');
nailfold = nailfold(:,:,1);
nailfold_mask = imerode(nailfold > 10 & nailfold < 250, strel('disk', 20));
[mag_2d] = gaussian_clover_line(nailfold, 4);
[mag_1d] = gaussian_1st_derivative_gradient2(nailfold, 2);

template_apex1 = mean(mag_patches1,3);
C1 = normxcorr2(template_apex1, mag_1d);
C1 = C1(17:end-16, 17:end-16);

template_apex2 = mean(mag_patches2,3);
C2 = normxcorr2(template_apex2, mag_2d);
C2 = C2(17:end-16, 17:end-16);

[maxima_pos1, maxima_vals1] = local_image_maxima(C1, 33, nailfold_mask, 0.5);
[maxima_pos2, maxima_vals2] = local_image_maxima(C2, 33, nailfold_mask, 0.5);

figure; 
subplot(3,1,1); imagesc(nailfold); axis image; hold on; colormap(gray(256));
plot(maxima_pos1(:,1), maxima_pos1(:,2), 'r+');
plot(maxima_pos2(:,1), maxima_pos2(:,2), 'gx');

subplot(3,1,2); imagesc(C1); axis image; hold on;
subplot(3,1,3); imagesc(C2); axis image; hold on;

%%
figure;
subplot(2,1,1);
hold on; axis ij equal;
for ii = 1:num_vessels
    plot(vessels{ii}(:,1), vessels{ii}(:,2));
    plot(vessel_top(ii,1), vessel_top(ii,2), 'rx');       
end
plot(maxima_pos1(:,1), maxima_pos1(:,2), 'g^');
plot(maxima_pos2(:,1), maxima_pos2(:,2), 'gv');

subplot(2,1,2); imagesc(nailfold); axis image; hold on; colormap(gray(256));
plot(maxima_pos1(:,1), maxima_pos1(:,2), 'g^');
plot(maxima_pos2(:,1), maxima_pos2(:,2), 'gv');
%%
C12 = C1 .* C2;
C12(~nailfold_mask) = 0;
[maxima_pos12, maxima_vals12] = local_image_maxima(C12, 33, nailfold_mask, 0.2);

%
figure;
subplot(2,1,1);
hold on; axis ij equal;
for ii = 1:num_vessels
    plot(vessels{ii}(:,1), vessels{ii}(:,2));
    plot(vessel_top(ii,1), vessel_top(ii,2), 'rx');       
end
plot(maxima_pos12(:,1), maxima_pos12(:,2), 'gx');

subplot(2,1,2); imagesc(C12); axis image; hold on; colormap(gray(256));
plot(maxima_pos12(:,1), maxima_pos12(:,2), 'gx');

n_min = min(nailfold(nailfold_mask));
n_max = max(nailfold(nailfold_mask));

colors = hot(256);
color_idx = linspace(min(maxima_vals12), max(maxima_vals12), 256);

figure; imagesc(C12); axis image off; colormap(hot(256)); caxis([0 color_idx(end)]); 
figure; imagesc(nailfold); axis image off; hold on; colormap(gray(256)); caxis([n_min n_max]);
for ii = 1:length(maxima_vals12)
    ms = 1 + 20*(maxima_vals12(ii)-0.2);
    %mc = interp1(color_idx, colors, maxima_vals12(ii));
    %plot(maxima_pos12(ii,1), maxima_pos12(ii,2), 'x', 'MarkerEdgeColor', mc);
    plot(maxima_pos12(ii,1), maxima_pos12(ii,2), 'r+', 'MarkerSize', ms);
end

figure; imagesc(color_idx); colormap(hot(256)); colorbar('location', 'southoutside');
%%
n_files = dir([nailfoldroot 'images\anonymous_oct\bmp\*.bmp']);
do_plot = 0;
for nn = 1:40
    %load in nailfold
    nailfold = imread([nailfoldroot 'images\anonymous_oct\bmp\' n_files(nn).name]);
    nailfold = nailfold(:,:,1);
    nailfold_mask = imerode(nailfold > 10 & nailfold < 250, strel('disk', 20));

    load([nailfoldroot 'images\anonymous_oct\apexes\' n_files(nn).name(1:end-4) '_apexes.mat'], 'maxima_*');
    
    %Plot the maxima in the nailfold image, with marker size scaled according
    %to NCC value
    n_min = min(nailfold(nailfold_mask));
    n_max = max(nailfold(nailfold_mask)); 
    figure; imagesc(nailfold); axis image off; hold on; colormap(gray(256)); caxis([n_min n_max]);
    for ii = 1:length(maxima_vals12)
        ms = 1 + 20*(maxima_vals12(ii)-0.2);
        plot(maxima_pos12(ii,1), maxima_pos12(ii,2), 'r+', 'MarkerSize', ms);
    end
    
end
%%
load([nailfoldroot 'data\apex_detection\apex_gd_templates.mat'], 'apex_template*');

n_files = dir([nailfoldroot 'images\anonymous_oct\bmp\*.bmp']);
do_plot = 1;
for nn = 4
    %load in nailfold
    nailfold = imread([nailfoldroot 'images\anonymous_oct\bmp\' n_files(nn).name]);
    nailfold = nailfold(:,:,1);
    nailfold_mask = imerode(nailfold > 10 & nailfold < 250, strel('disk', 20));

    %Compute the gaussian derivatives
    [mag_2d] = gaussian_clover_line(nailfold, 4);
    [mag_1d] = gaussian_1st_derivative_gradient2(nailfold, 2);

    %Apply normalised cross correlations using the 2 templates
    C1 = normxcorr2(apex_template1, mag_1d);
    C1 = C1(17:end-16, 17:end-16); %Need to discard floor(patch_size)/2 from start and end

    C2 = normxcorr2(apex_template2, mag_2d);
    C2 = C2(17:end-16, 17:end-16); %Need to discard floor(patch_size)/2 from start and end

    %Look for local maxima in the 2 independent NCC maps
    [maxima_pos1, maxima_vals1] = local_image_maxima(C1, 33, nailfold_mask, 0.5);
    [maxima_pos2, maxima_vals2] = local_image_maxima(C2, 33, nailfold_mask, 0.5);

    %Now multiply the 2 NCC maps to get a single map and again look for local
    %maxima
    C12 = C1 .* C2;
    %C12(~nailfold_mask) = 0;
    [maxima_pos12, maxima_vals12] = local_image_maxima(C12, 33, nailfold_mask, 0.2);
    
    if do_plot
        figure; 
        %Display a heat map of the combined tempate matching
        subplot(2,1,1); imagesc(C12); axis image off; colormap(hot(256)); caxis([0 max(maxima_vals12)]);

        %Plot the maxima in the nailfold image, with marker size scaled according
        %to NCC value
        n_min = min(nailfold(nailfold_mask));
        n_max = max(nailfold(nailfold_mask)); 
        subplot(2,1,2); imagesc(nailfold); axis image off; hold on; colormap(gray(256)); caxis([n_min n_max]);
        for ii = 1:length(maxima_vals12)
            ms = 1 + 20*(maxima_vals12(ii)-0.2);
            plot(maxima_pos12(ii,1), maxima_pos12(ii,2), 'r+', 'MarkerSize', ms);
        end
    end
end
%%

c12_bw = C12 > 0.1;
c12_mask = false(size(c12_bw));
for deg = 0:15:165
    c12_mask = c12_mask | imopen(c12_bw, strel('line', 32, deg));
    figure; imagesc(imopen(c12_bw, strel('line', 16, deg))); axis image; colormap gray; hold on;
end
figure; 
subplot(2,1,1); imagesc(c12_mask); axis image;
[r c] = find(c12_mask);
subplot(2,1,2); imagesc(bwselect(c12_bw, c, r,8)); axis image;

[maxima_pos12a, maxima_vals12a] = local_image_maxima(C12, 33, c12_bw & ~bwselect(c12_bw, c, r, 8), 0.2);

figure; imagesc(C12 > 0.1); axis image; colormap gray; hold on;
plot(maxima_pos12(:,1), maxima_pos12(:,2), 'r+');
plot(maxima_pos12a(:,1), maxima_pos12a(:,2), 'gx');

figure; imagesc(nailfold); axis image; colormap gray; hold on;
plot(maxima_pos12a(:,1), maxima_pos12a(:,2), 'gx');



        
