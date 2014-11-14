%--------------------------------------------------------------------------
% Script demonstrating template matching of vessel apexes in nailfold
% images
%--------------------------------------------------------------------------

load('C:\isbe\nailfold\playground\model_experiments\vessel_names.mat', 'train_names');

num_vessels = length(train_names);
gd1_patches = zeros(33,33,num_vessels);
gd2_patches = zeros(33,33,num_vessels);

for ii = 1:num_vessels
       
    %Sample patch from 1st deriv.
    mag_1d = u_load(['C:\isbe\nailfold\playground\model_experiments\images\g1d\vessel_top' train_names{ii} '.mat']);

    %Sample patch from 1st deriv.
    mag_2d = u_load(['C:\isbe\nailfold\playground\model_experiments\images\g2d\vessel_top' train_names{ii} '.mat']);

    %Sample patch from 1st deriv.
    %image_patch = imread(['C:\isbe\nailfold\playground\model_experiments\images\orig\vessel_top' train_names{ii} '.png']);

    %Load in pts
    v_pts = u_load(['C:\isbe\nailfold\playground\model_experiments\points_m\vessel_top' train_names{ii} '.mat']);

    gd1_patches(:,:,ii) = abs(sample_window(mag_1d, 33, round(v_pts(16,2)), round(v_pts(16,1))));
    gd2_patches(:,:,ii) = sample_window(mag_2d, 33, round(v_pts(16,2)), round(v_pts(16,1)));
end

%Compute the mean of the 2 sets of patches to use as templates
apex_template1 = mean(gd1_patches,3);
apex_template2 = mean(gd2_patches,3);

%save the templates for later use
save([nailfoldroot 'data\apex_detection\apex_gd_templates.mat'], 'apex_template*');

figure; imagesc(apex_template1); axis image;
title('Template for 1st derivative patches');
figure; imagesc(apex_template2); axis image;
title('Template for 2nd derivative patches');
%%
load('C:\isbe\nailfold\playground\model_experiments\vessel_names.mat', 'train_names');

num_vessels = length(train_names);
u_shapes = zeros(num_vessels, 62);
areas = zeros(num_vessels, 1);

for ii = 1:num_vessels
    %Load in pts
    v_pts = u_load(['C:\isbe\nailfold\playground\model_experiments\points_m\vessel_top' train_names{ii} '.mat']);
    u_shapes(ii,:) = v_pts(:)';
    areas(ii) = polyarea(v_pts(:,1), v_pts(:,2));
end
[a_shapes, a_scales, mean_target, a_rots, a_trans a_origins] = align_shapes(u_shapes, 'area', mean(areas));
mean_shape = [mean(a_shapes(:,1:31))', mean(a_shapes(:,32:62))'];
theta = mean(acos(squeeze(a_rots(1,1,:))));
theta_rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];


figure; hold on;
plot(a_shapes(:,1:31)', a_shapes(:,32:62)'); axis ij equal;
plot(mean_shape(:,1), mean_shape(:,2), 'k', 'linewidth', 2)

mean_shape = mean_target * theta_rot';
mean_shape = [mean_shape(:,1)-mean_shape(16,1) mean_shape(:,2)-mean_shape(16,2)];
figure; plot(mean_shape(:,1), mean_shape(:,2)); axis ij equal;
%%
x = repmat(-24:24, 49, 1);
y = x';
xy = [x(:) y(:)];

gd1_patches_a = zeros(49,49,num_vessels);
gd2_patches_a = zeros(49,49,num_vessels);
for ii = 1:num_vessels
    %Load in pts
    v_pts = u_load(['C:\isbe\nailfold\playground\model_experiments\points_m\vessel_top' train_names{ii} '.mat']);
    
    %Sample patch from 1st deriv.
    mag_1d = abs(u_load(['C:\isbe\nailfold\playground\model_experiments\images\g1d\vessel_top' train_names{ii} '.mat']));

    %Sample patch from 1st deriv.
    mag_2d = u_load(['C:\isbe\nailfold\playground\model_experiments\images\g2d\vessel_top' train_names{ii} '.mat']);
    
    xya = (xy * a_rots(:,:,ii)' / a_scales(ii))* theta_rot;
    xa = reshape(xya(:,1) + v_pts(16,1), 49, 49);
    ya = reshape(xya(:,2) + v_pts(16,2), 49, 49);
    
    gd1_patches_a(:,:,ii) = interp2(mag_1d, xa, ya, 'bilinear');
    gd2_patches_a(:,:,ii) = interp2(mag_2d, xa, ya, 'bilinear');
    
%     figure;
%     subplot(2,2,1); imgray(mag_1d); plot(xa, ya, 'rx'); plot(v_pts(:,1), v_pts(:,2), 'g'); plot(v_pts(16,1), v_pts(16,2), 'co');
%     subplot(2,2,2); imgray(mag_2d); plot(xa, ya, 'rx'); plot(v_pts(:,1), v_pts(:,2), 'g'); plot(v_pts(16,1), v_pts(16,2), 'co');
%     subplot(2,2,3); imgray(gd1_patches_a(:,:,ii));
%     subplot(2,2,4); imgray(gd2_patches_a(:,:,ii));

end

%Compute the mean of the 2 sets of patches to use as templates
apex_template1_r = zeros(49,49, 11);
apex_template2_r = zeros(49,49, 11);
apex_template1_a = mean(gd1_patches_a,3);
apex_template2_a = mean(gd2_patches_a,3);
circle_mask = x.^2 + y.^2 > 24^2;
%%
figure; 
subplot(1,2,1); imagesc(apex_template1_a); axis image; hold on;
plot(mean_shape(:,1)+24, mean_shape(:,2)+24, 'k', 'linewidth', 2);
title('Template for 1st derivative patches');
subplot(1,2,2); imagesc(apex_template2_a); axis image; hold on;
plot(mean_shape(:,1)+24, mean_shape(:,2)+24, 'k', 'linewidth', 2);
title('Template for 2nd derivative patches');
    
    
for ii = 1:11
    theta = (ii - 6)*6;
    at1a = imrotate(apex_template1_a, theta, 'bilinear', 'crop');
    at1a(circle_mask) = NaN;
    at2a = imrotate(apex_template2_a, theta, 'bilinear', 'crop');
    at2a(circle_mask) = NaN;
    
    apex_template1_r(:,:,ii) = at1a;
    apex_template2_r(:,:,ii) = at2a;
    
    figure; 
    subplot(1,2,1); imagesc(apex_template1_r(:,:,ii)); axis image;
    title('Template for 1st derivative patches');
    subplot(1,2,2); imagesc(apex_template2_r(:,:,ii)); axis image;
    title('Template for 2nd derivative patches');

end

%save the templates for later use
save([nailfoldroot 'data\apex_detection\apex_gd_templates.mat'], 'apex_template*');

%%
load([nailfoldroot 'data\apex_detection\apex_gd_templates.mat'], 'apex_template*');
nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');

for nn = 8:12
    image_path = ['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp'];
    nailfold = imread(image_path);   
    nailfold = nailfold(:,:,1);
    nailfold_mask = imerode(nailfold > 10 & nailfold < 250, strel('disk', 20));

    %Compute the gaussian derivatives
    [mag_2d] = gaussian_2nd_derivative_line(nailfold, 4);
    [mag_1d] = abs(gaussian_1st_derivative_gradient(nailfold, 2));

    %Apply normalised cross correlations using the 2 templates
%     C1 = normxcorr2(apex_template1, mag_1d);
%     C1 = C1(17:end-16, 17:end-16); %Need to discard floor(patch_size)/2 from start and end
% 
%     C2 = normxcorr2(apex_template2, mag_2d);
%     C2 = C2(17:end-16, 17:end-16); %Need to discard floor(patch_size)/2 from start and end

    C1 = -inf(size(nailfold));
    C2 = -inf(size(nailfold));
    
    for ii = 6%1:11
    
        c1 = mb_normxcorr2(apex_template1_r(:,:,ii), mag_1d);
        c2 = mb_normxcorr2(apex_template2_r(:,:,ii), mag_2d);
        
        C1 = max(C1, c1);
        C2 = max(C1, c2);
    end

    %Look for local maxima in the 2 independent NCC maps
%     [maxima_pos1, maxima_vals1] = local_image_maxima(C1, 33, nailfold_mask, 0.5);
%     [maxima_pos2, maxima_vals2] = local_image_maxima(C2, 33, nailfold_mask, 0.5);

    %Now multiply the 2 NCC maps to get a single map and again look for local
    %maxima
    C12 = C1 .* C2;
    C12(~nailfold_mask) = 0;
    [maxima_pos12, maxima_vals12] = local_image_maxima(C12, 33, [], 0.2);
         
    %Display a heat map of the combined tempate matching
    %figure; imagesc(C12); axis image off; colormap(hot(256)); caxis([0 max(maxima_vals12)]);

    %Plot the maxima in the nailfold image, with marker size scaled according
    %to NCC value
    n_min = min(nailfold(nailfold_mask));
    n_max = max(nailfold(nailfold_mask)); 
    figure; imagesc(nailfold); axis image off; hold on; colormap(gray(256)); caxis([n_min n_max]);
%     plot(maxima_pos1(:,1), maxima_pos1(:,2), 'gx');
%     plot(maxima_pos2(:,1), maxima_pos2(:,2), 'co');
%     plot(maxima_pos12(:,1), maxima_pos12(:,2), 'rs');
    for ii = 1:length(maxima_vals12)
        ms = 1 + 20*(maxima_vals12(ii)-0.2);
        plot(maxima_pos12(ii,1), maxima_pos12(ii,2), 'r+', 'MarkerSize', ms);
    end
    
    save([nailfoldroot 'data\apex_detection\' nf_files(nn).name(1:end-4) '_apexes.mat'], 'maxima_*');
%     load([nailfoldroot 'images\anonymous_oct\apexes\' nf_files(nn).name(1:end-11) '_apexes.mat']);
%     for ii = 1:length(maxima_vals12)
%         ms = 1 + 20*(maxima_vals12(ii)-0.2);
%         plot(maxima_pos12(ii,1), maxima_pos12(ii,2), 'go', 'MarkerSize', ms);
%     end
end

%%
load([nailfoldroot 'data\apex_detection\apex_gd_templates.mat'], 'apex_template*');
nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');

nn = 8;
image_path = ['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp'];
nailfold = imread(image_path);   
nailfold = nailfold(:,:,1);
[rows cols] = size(nailfold);
[mag_1d] = abs(gaussian_1st_derivative_gradient(nailfold, 2));
clear nailfold;
    
C1 = -inf(rows, cols);
for jj = 1:11
    for ii = -5:5
        scale = 1 + ii/10;
        mag_1ds = imresize(mag_1d, scale, 'bilinear');
        c1 = imresize(mb_normxcorr2(apex_template1_a(:,:,jj), mag_1ds), [rows cols], 'bilinear');
        C1 = max(C1, c1);
        %figure; imgray(c1); colormap jet; caxis([-1 1]);
    end
end
figure; imgray(C1); colormap jet; caxis([-1 1]);
%%
nf_files = dir('C:\isbe\nailfold\data\apex_detection\centre_rotation\*.mat');

nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');

for nn = 8:12
    %load image
    image_path = ['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp'];
    nailfold = imread(image_path);   
    nailfold = nailfold(:,:,1);
    
    %load potential positions
    load([nailfoldroot 'data\apex_detection\' nf_files(nn).name(1:end-4) '_apexes.mat'], 'maxima_*');
    
    for ii = 1:length(maxima_vals12)
        sr = max(1, floor(min(vessel_top(:,2)) - 50));
        er = min(rows, ceil(max(vessel_top(:,2)) + 50));
        sc = max(1, floor(min(vessel_top(:,1)) - 50));
        ec = min(cols, floor(max(vessel_top(:,1)) + 50));
    end
end

%%
nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');

for nn = 8:12
    image_path = ['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp'];
    nailfold = imread(image_path);   
    nailfold = nailfold(:,:,1);
    nailfold_mask = imerode(nailfold > 10 & nailfold < 250, strel('disk', 20));

    load(['C:\isbe\nailfold\data\apexes\' nf_files(nn).name(1:end-4) '_apexes.mat'], 'apexes');
    
    %Compute the gaussian derivatives
    [mag_2d ori_2d] = gaussian_2nd_derivative_line(nailfold, 4);
    vessel_nms = mb_non_maximal_supp(mag_2d, ori_2d);
    vessel_nms(~nailfold_mask) = 0;
    [yv xv] = find(vessel_nms);
    
    n_min = min(nailfold(nailfold_mask));
    n_max = max(nailfold(nailfold_mask)); 
    figure; imgray(nailfold); caxis([n_min n_max]);
    plot(xv, yv, 'g.', 'markersize', 2);
    
    for ii = 1:length(apexes)
        plot(apexes{ii}(:,1), apexes{ii}(:,2), 'r');
    end
end

%%
if ~isdir([nailfoldroot 'data\apex_detection\']); mkdir([nailfoldroot 'data\apex_detection\']); end
if ~isdir([nailfoldroot 'images\anonymous_oct\apexes\']); mkdir([nailfoldroot 'images\anonymous_oct\apexes\']); end
%% 1) Locate a set of apexes in 1 image given manual annotations
% We'll use MB's markup of the KK demo image (because it just happens to be
% an image I've marked up lots of vessels on)

%load in the nailfold - note we need to resize this by a factor 3/4 to
%match the other nailfolds in our set
nailfold = imresize(imread([nailfoldroot 'images\kk_demo\nailfolds\n3_mb.bmp']), 0.75, 'bilinear');

%Load in the annotated vessels
v = u_load([nailfoldroot 'images\kK_demo\nailfolds\n3_mb_vessels_all.mat']);
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

%Now choose each vessel apex as the point with highest y-val (I know this
%may not strictly be accurate given local deviations, rotations etc - but
%it appears to work more robustly than more complex methods such as
%searching for points of maximum curvature etc)
figure; hold on; axis equal ij;
vessel_apex = zeros(num_vessels, 2);
for ii = 1:num_vessels
    [y_min min_idx] = min(vessels{ii}(:,2));
    vessel_apex(ii,:) = [vessels{ii}(min_idx,1) y_min];
    plot(vessels{ii}(:,1), vessels{ii}(:,2)); 
    plot(vessel_apex(ii,1), vessel_apex(ii,2), 'rx');          
end

%% 2) Now compute Gaussian 1st + 2nd derivatives of the nailfold and sample
% a patch of each about each apex

%Compute the Gaussian derivatives
[mag_2d, ori_2d] = gaussian_clover_line(nailfold, 4);
[mag_1d] = gaussian_1st_derivative_gradient2(nailfold, 2);

f1 = figure;
f2 = figure;

gd1_patches = zeros(33,33,num_vessels);
gd2_patches = zeros(33,33,num_vessels);
for ii = 1:num_vessels
    
    %Sample patch from 1st deriv.
    mag_patch = sample_window(mag_1d, 33, round(vessel_apex(ii,2)), round(vessel_apex(ii,1)));
    mag_patch = (mag_patch - mean(mag_patch(:))) / std(mag_patch(:)); %Normalise patch
    gd1_patches(:,:,ii) = mag_patch;
    
    %Sample patch from 1st deriv.
    mag_patch = sample_window(mag_2d, 33, round(vessel_apex(ii,2)), round(vessel_apex(ii,1)));
    mag_patch = (mag_patch - mean(mag_patch(:))) / std(mag_patch(:));  %Normalise patch
    gd2_patches(:,:,ii) = mag_patch;

    if ii <= 12
        figure(f1);
        subplot(3,4,ii); imagesc(gd1_patches(:,:,ii)); axis image; colormap(gray(256));
        figure(f2);
        subplot(3,4,ii); imagesc(gd2_patches(:,:,ii)); axis image; colormap(gray(256));
    end
end

%Compute the mean of the 2 sets of patches to use as templates
apex_template1 = mean(gd1_patches,3);
apex_template2 = mean(gd2_patches,3);

%save the templates for later use
save([nailfoldroot 'data\apex_detection\apex_gd_templates.mat'], 'apex_template*');

figure; imagesc(apex_template1); axis image;
title('Template for 1st derivative patches');
figure; imagesc(apex_template2); axis image;
title('Template for 2nd derivative patches');

%% 3) Use the template to try and find vessel apexes in a new nailfold

%load in nailfold
nailfold = imread([nailfoldroot 'images\ncm\andrea Murray.bmp']);
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

%Display the maps and plot the local maxima
figure; 
subplot(3,1,1); imagesc(nailfold); axis image; hold on; colormap(gray(256));
plot(maxima_pos1(:,1), maxima_pos1(:,2), 'r+');
plot(maxima_pos2(:,1), maxima_pos2(:,2), 'gx');
subplot(3,1,2); imagesc(C1); axis image; hold on;
subplot(3,1,3); imagesc(C2); axis image; hold on;

%Now multiply the 2 NCC maps to get a single map and again look for local
%maxima
C12 = C1 .* C2;
C12(~nailfold_mask) = 0;
[maxima_pos12, maxima_vals12] = local_image_maxima(C12, 33, [], 0.2);

%Display a heat map of the combined tempate matching
figure; imagesc(C12); axis image off; colormap(hot(256)); caxis([0 max(maxima_vals12)]);

%Plot the maxima in the nailfold image, with marker size scaled according
%to NCC value
n_min = min(nailfold(nailfold_mask));
n_max = max(nailfold(nailfold_mask)); 
figure; imagesc(nailfold); axis image off; hold on; colormap(gray(256)); caxis([n_min n_max]);
for ii = 1:length(maxima_vals12)
    ms = 1 + 20*(maxima_vals12(ii)-0.2);
    plot(maxima_pos12(ii,1), maxima_pos12(ii,2), 'r+', 'MarkerSize', ms);
end
%%
%Repeat the above test for each nailfold in the anonymous OCT set

load([nailfoldroot 'data\apex_detection\apex_gd_templates.mat'], 'apex_template*');

n_files = dir([nailfoldroot 'images\anonymous_oct\bmp\*.bmp']);
do_plot = 0;
for nn = 15
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
    C12(~nailfold_mask) = 0;
    [maxima_pos12, maxima_vals12] = local_image_maxima(C12, 33, [], 0.2);
    
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
    
    %%
    save([nailfoldroot 'images\anonymous_oct\apexes\' n_files(nn).name(1:end-4) '_apexes.mat'], 'maxima_*');
end
%%