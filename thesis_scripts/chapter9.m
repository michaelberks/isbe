%--------------------------------------------------------------------------
%------------- Chapter 9 script -------------------------------------------
%--------------------------------------------------------------------------

%%
%--------------------------------------------------------------------------
%Build mass distribution model - CC

%Get list of CC segmentations
cc_list = dir('C:\isbe\dev\segmentation\breast_borders\*CC*.mat');
mass_cc_list = dir('C:\isbe\dev\annotations\*CC*.mat');

%Get standardised CC shapes from segmentations
[cc_shapes, cc_areas] = get_cc_shape(cc_list, 50, 0);

%Compute the mean CC breast shape
[cc_mean_shape] = compute_cc_mean(cc_shapes, cc_areas, 0);
%%
% Now we have a target mean CC shape, we can map the position of mass into
% this shape
[mass_centres] = warp_mass_centres(cc_mean_shape, mass_cc_list, 0, 0);

%Plot the mass centres in the mean shape
figure; plot(cc_mean_shape(:,1), cc_mean_shape(:,2)); axis ij equal; hold on;
plot(mass_centres(:,1), mass_centres(:,2), 'r+');

%Compute the location distribution model
[location_model] = build_mass_distribution_model(cc_mean_shape, mass_centres, 1);

%Show some samples from the model
for ii = 1:9
    thresh = (ii/10);
    sample_mass_location(location_model,[], thresh, 1);
end
%
mkdir C:\isbe\dev\location
save C:\isbe\dev\location\mean_CC_shape cc_mean_shape
save C:\isbe\dev\location\CC_shapes cc_shapes
save C:\isbe\dev\location\CC_centres mass_centres
save C:\isbe\dev\location\CC_location_model location_model
%%
%--------------------------------------------------------------------------
% Build mass distribution model - MLO
%

%Get list of ML segmentations
ml_list = dir('C:\isbe\dev\segmentation\breast_borders\*ML*.mat');
mass_ml_list = dir('C:\isbe\dev\annotations\*ML*.mat');

%Get standardised ML shapes from segmentations
[ml_shapes, ml_areas] = get_mlo_shape(ml_list, 50, 0);
%
%Compute the mean ML breast shape
[ml_mean_shape] = compute_cc_mean(ml_shapes, ml_areas, 0);

% Now we have a target mean ML shape, we can map the position of mass into
% this shape
[mass_centres] = warp_mass_centres(ml_mean_shape, mass_ml_list, 1, 0);

%Plot the mass centres in the mean shape
figure; plot(ml_mean_shape(:,1), ml_mean_shape(:,2)); axis ij equal; hold on;
plot(mass_centres(:,1), mass_centres(:,2), 'r+');

%Compute the location distribution model
[location_model] = build_mass_distribution_model(ml_mean_shape, mass_centres, 1);
%
%Show some samples from the model
for ii = 1:9
    thresh = (ii/10);
    sample_mass_location(location_model,[], thresh, 1);
end
%%
save C:\isbe\dev\location\mean_ML_shape ml_mean_shape
save C:\isbe\dev\location\ML_shapes ml_shapes
save C:\isbe\dev\location\ML_centres mass_centres
save C:\isbe\dev\location\ML_location_model location_model

%%
%--------------------------------------------------------------------------
% Show figure of mass centres in mean shape
%--------------------------------------------------------------------------
load C:\isbe\dev\location\mean_CC_shape
load C:\isbe\dev\location\CC_centres

apos = [40 40 max(cc_mean_shape(:,1)) - min(cc_mean_shape(:,1)) + 20, max(cc_mean_shape(:,2)) - min(cc_mean_shape(:,2)) + 20];
xlims = [min(cc_mean_shape(:,1)) - 10 max(cc_mean_shape(:,1)) + 10];
ylims = [min(cc_mean_shape(:,2)) - 10 max(cc_mean_shape(:,2)) + 10];
fpos = [0 0 apos(3) + 80 apos(4) + 80];

f1 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', fpos,...
    'PaperPositionMode','auto');
axes(...
    'Units', 'pixels',...
    'position', apos,...
    'xlim', xlims,...
    'ylim', ylims); 

axis ij equal off; hold on;

plot([cc_mean_shape(:,1); cc_mean_shape(1,1)] , [cc_mean_shape(:,2); cc_mean_shape(1,2)], 'b', 'LineWidth', 2);
plot(cc_mean_shape(51,1), cc_mean_shape(51,2), 'k.');
plot(mass_centres(:,1), mass_centres(:,2), 'r+', 'MarkerSize', 10);
print('-dtiff', '-noui', '-painters', f1, '-r300', 'C:\isbe\thesis\figures\10\mass_centres_cc.tif');
%%
load C:\isbe\dev\location\mean_ML_shape
load C:\isbe\dev\location\ML_centres

apos = [40 40 max(ml_mean_shape(:,1)) - min(ml_mean_shape(:,1)) + 20, max(ml_mean_shape(:,2)) - min(ml_mean_shape(:,2)) + 20];
xlims = [min(ml_mean_shape(:,1)) - 10 max(ml_mean_shape(:,1)) + 10];
ylims = [min(ml_mean_shape(:,2)) - 10 max(ml_mean_shape(:,2)) + 10];
fpos = [0 0 apos(3) + 80 apos(4) + 80];

f2 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', fpos,...
    'PaperPositionMode','auto');
axes(...
    'Units', 'pixels',...
    'position', apos,...
    'xlim', xlims,...
    'ylim', ylims); 

axis ij equal off; hold on;

plot([ml_mean_shape(:,1); ml_mean_shape(1,1)] , [ml_mean_shape(:,2); ml_mean_shape(1,2)], 'b', 'LineWidth', 2);
plot(ml_mean_shape(51,1), ml_mean_shape(51,2), 'k.');
plot(mass_centres(:,1), mass_centres(:,2), 'r+', 'MarkerSize', 10);
print('-dtiff', '-noui', '-painters', f2, '-r300', 'C:\isbe\thesis\figures\10\mass_centres_mlo.tif');

%%
plot(ml_mean_shape(:,1), ml_mean_shape(:,2)); 
plot(mass_centres(:,1), mass_centres(:,2), 'r+');

load C:\isbe\dev\location\mean_ML_shape
load C:\isbe\dev\location\ML_shapes
load C:\isbe\dev\location\ML_centres
load C:\isbe\dev\location\ML_location_model
%%
%--------------------------------------------------------------------------
% Show the adaptive and fixed kernel distributions of mass location
%--------------------------------------------------------------------------
load C:\isbe\dev\location\ML_location_model
%load C:\isbe\dev\location\CC_location_model %repeat for CC

%Get start and end (x,y)-coordinates to bound the mean breast shape
xs = floor(min(location_model.mean_shape(:,1)));
xe = ceil(max(location_model.mean_shape(:,1)));
ys = floor(min(location_model.mean_shape(:,2)));
ye = ceil(max(location_model.mean_shape(:,2)));

mask_size(1) = ye - ys + 1;
mask_size(2) = xe - xs + 1;

%Build a BW mask of the shape
[bw] = poly2mask(...
    location_model.mean_shape(:,1) - xs + 1,...
    location_model.mean_shape(:,2) - ys + 1,...
    mask_size(1),...
    mask_size(2));


%Extract a list of indices belonging to the breast shape in the mask
[in_idx] = find(bw);
%
D_map = nan(mask_size);
D_map(in_idx) = location_model.D_f;

D_map_a = nan(mask_size);
D_map_a(in_idx) = location_model.D_a;

D_map = padarray(D_map, [50 50], nan);
D_map_a = padarray(D_map_a, [50 50], nan);

figure;    
subplot(1,2,1); imagesc(D_map); axis image; colormap(jet(256));
subplot(1,2,2); imagesc(D_map_a); axis image; colormap(jet(256));

write_im_from_colormap(D_map, 'C:\isbe\thesis\figures\10\D_f_mlo.jpg', jet(256), [], [0 0 0]);
write_im_from_colormap(D_map_a, 'C:\isbe\thesis\figures\10\D_a_mlo.jpg', jet(256), [], [0 0 0]);

% write_im_from_colormap(D_map, 'C:\isbe\thesis\figures\10\D_f_mlo.jpg', jet(256), [], [0 0 0]); %repeat for CC
% write_im_from_colormap(D_map_a, 'C:\isbe\thesis\figures\10\D_a_mlo.jpg', jet(256), [], [0 0 0]); %repeat for CC

%%
%---------------------------------------------------------------------------
%%
mlo_list = dir('C:\isbe\dev\segmentation\breast_borders\*CC*.mat');
mkdir C:\isbe\dev\segmentation\mammograms
for ii = 1:length(mlo_list)
    segmentation = u_load(['C:\isbe\dev\segmentation\breast_borders\', mlo_list(ii).name]);
    mam_name = [mlo_list(ii).name(1:end-17), '.bmp'];
    mam = imread(['C:\isbe\mammograms\new_CAD\BMP_2004\', mam_name]);
    
    mam_small = imresize(mam, segmentation.size, 'bilinear');
    clear mam;
    save(['C:\isbe\dev\segmentation\mammograms\', mam_name(1:end-4), '.mat'], 'mam_small');
    clear mam_small;
end
%%
%--------------------------------------------------------------------------
% Produce mean histogram of intensities for all mammograms
%--------------------------------------------------------------------------
mam_list = dir('C:\isbe\dev\segmentation\mammograms\*CC*.mat');
all_counts = zeros(256,1);
for ii = 1: length(mam_list)
    mam = u_load(['C:\isbe\dev\segmentation\mammograms\', mam_list(ii).name]);
    counts = hist(mam(:), 0:255);
    all_counts = all_counts + counts(:);
    clear mam counts;
end
%%
f = figure(...
        'windowstyle', 'normal',...
        'Units', 'pixels',...
        'position', [100 100 600 600],...
        'PaperPositionMode','auto'); 
bar(0:255, all_counts, 'r'); set(gca, 'yticklabel', [], 'units', 'pixels', 'position', [10 30 580 560], 'xlim', [0 255]);
print('-dtiff', '-noui', '-painters', f, '-r300', 'C:\isbe\thesis\figures\10\mammogram_histogram_mean.tif');
%%
%--------------------------------------------------------------------------
% Show intensity histograms for specific mammograms
%--------------------------------------------------------------------------
% Use CC: 006LCC and MLO: 012ML
mam_cc = u_load('C:\isbe\dev\segmentation\mammograms\o04_006LCC.mat');
mam_ml = u_load('C:\isbe\dev\segmentation\mammograms\o04_012LML.mat');

f = figure(...
        'windowstyle', 'normal',...
        'Units', 'pixels',...
        'position', [100 100 600 600],...
        'PaperPositionMode','auto'); 
hist(mam_cc(:), 0:255); set(gca, 'yticklabel', [], 'units', 'pixels', 'position', [10 30 580 560], 'xlim', [0 255]);
print('-dtiff', '-noui', '-painters', f, '-r300', 'C:\isbe\thesis\figures\10\mammogram_histogram_cc.tif');

f = figure(...
        'windowstyle', 'normal',...
        'Units', 'pixels',...
        'position', [100 100 600 600],...
        'PaperPositionMode','auto'); 
hist(mam_ml(:), 0:255); set(gca, 'yticklabel', [], 'units', 'pixels', 'position', [10 30 580 560], 'xlim', [0 255]);
print('-dtiff', '-noui', '-painters', f, '-r300', 'C:\isbe\thesis\figures\10\mammogram_histogram_ml.tif');

%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------
% BW mask for finding inner breast edge
%--------------------------------------------------------------------------
% Use 014LML
%
mam = double(u_load('C:\isbe\dev\segmentation\mammograms\o04_014LML.mat'));
[rows, cols] = size(mam);
[counts, x_range] = hist(mam(:), 128);

%Find the index of the maximal bin in the first half of the intensity range
[dummy max_idx] = max(counts(1:64));

%As the intensity increases, the bin counts decrease - find the point at
%which they reach a minimum before increasing again, and use this as an
%upper threshold
counts(1:max_idx-1) = [];
x_range(1:max_idx-1) = [];
min_idx = find(diff(counts) >= 0, 1);
upper_threshold = ceil(x_range(min_idx));

%Get binary masks based on the threshold
upper_mask = mam > upper_threshold;

figure; imagesc(upper_mask); axis image; colormap gray
imwrite(upper_mask, 'C:\isbe\thesis\figures\10\bw_mask_initial.bmp');

[ll, no_objects] = bwlabel(upper_mask, 4); %#ok
[nn] = hist(ll(:), no_objects + 1);
[mm ind] = max(nn(2:end));
upper_mask = ll == ind;

%Close the mask (to fill in any holes inside the main region)
upper_mask = imdilate(upper_mask, strel('disk', 11));
upper_mask = imerode(upper_mask, strel('disk', 11));

%Open the mask to remove the straight line border artifacts
upper_mask = imerode(upper_mask, strel('disk', 75));
upper_mask = imdilate(upper_mask, strel('disk', 75));

%Again, throw away any regions we have now disconnected - note the breast
%region may no longer be the largest region, so we assume the point
%half-way up the chestwall will always belong to the main breast region
[ll] = bwlabel(upper_mask, 4);
upper_mask = ll == ll(round(rows/2), 1);

figure; imagesc(upper_mask); axis image; colormap gray
imwrite(upper_mask, 'C:\isbe\thesis\figures\10\bw_mask_processed.bmp');
%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------
% Show normal profiles and edge energy profiles
%--------------------------------------------------------------------------
% Use 014LML
mam = u_load('C:\isbe\dev\segmentation\mammograms\o04_014LML.mat');
[segmentation normal] = segment_breast('image', mam, 'mlo', 1, 'right', 0, 'plot', 1);

%Get edge energies
edge = imfilter(normal.profiles, [1 2 4 0 -4 -2 -1], 'replicate');
edge_n2 = imfilter(normal.profiles, [1 2 4 -2 -4 -2 -1], 'replicate');
edge_n2(edge_n2 < -150) = -150;

edge = (edge - min(edge(:))) / (max(edge(:)) - min(edge(:)));
edge_n2 = (edge_n2 - min(edge_n2(:))) / (max(edge_n2(:)) - min(edge_n2(:)));

%Display edge energies
figure; 
subplot(1,2,1); imagesc(edge); colormap(jet(256));
subplot(1,2,2); imagesc(edge_n2); colormap(jet(256));

%Save edge energies as 3D plots
f1 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [100 100 800 800],...
    'PaperPositionMode','auto'); 
surf(edge); set(gca, 'yticklabel', [], 'CameraPosition', [-76.3898 -145.5122 8.4115]);
print('-dtiff', '-noui', '-painters', f1, '-r300', 'C:\isbe\thesis\figures\10\initial_edge_energy.tif');

f2 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [100 100 800 800],...
    'PaperPositionMode','auto');
surf(edge_n2); set(gca, 'yticklabel', [], 'CameraPosition', [-76.3898 -145.5122 8.4115]);
print('-dtiff', '-noui', '-painters', f2, '-r300', 'C:\isbe\thesis\figures\10\final_edge_energy.tif');

%%
%--------------------------------------------------------------------------
% Show segmentation results for 1 CC and 1 MLO mammogram
%--------------------------------------------------------------------------
%%
% Use 014LCC and 014LML
mam = u_load('C:\isbe\dev\segmentation\mammograms\o04_014LCC.mat');
segmentation = u_load('C:\isbe\dev\segmentation\breast_borders\o04_014LCC_segmentation.mat');

f3 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [100 100 size(mam,2)/2 size(mam,1)/2],...
    'PaperPositionMode','auto');
axes('units', 'pixels', 'position', [0 0 size(mam,2)/2 size(mam,1)/2]);
imagesc(mam); colormap(gray(256));  axis image; hold on; 
plot(segmentation.breast_border(segmentation.breast_air,1), segmentation.breast_border(segmentation.breast_air,2), 'r', 'LineWidth', 1.5);
set(gca, 'xticklabel', [], 'yticklabel', []);

print('-dtiff', '-noui', '-painters', f3, '-r300', 'C:\isbe\thesis\figures\10\segmented_mammogram_cc.tif');
%%
mam = u_load('C:\isbe\dev\segmentation\mammograms\o04_014LML.mat');
segmentation = u_load('C:\isbe\dev\segmentation\breast_borders\o04_014LML_segmentation.mat');

f5 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [100 100 size(mam,2)/2 size(mam,1)/2],...
    'PaperPositionMode','auto');
axes('units', 'pixels', 'position', [0 0 size(mam,2)/2 size(mam,1)/2]);
imagesc(mam); colormap(gray(256));  axis image; hold on; 
plot(segmentation.breast_border(segmentation.breast_air,1), segmentation.breast_border(segmentation.breast_air,2), 'r', 'LineWidth', 1.5);
set(gca, 'xticklabel', [], 'yticklabel', []);

print('-dtiff', '-noui', '-painters', f5, '-r300', 'C:\isbe\thesis\figures\10\segmented_mammogram_mlo.tif');
%%
% Use 014LCC and 014LML
mam = u_load('C:\isbe\dev\segmentation\mammograms\o04_031RCC.mat');
segmentation = u_load('C:\isbe\dev\segmentation\breast_borders\o04_031RCC_segmentation.mat');

f3 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [100 100 size(mam,2)/2 size(mam,1)/2],...
    'PaperPositionMode','auto');
axes('units', 'pixels', 'position', [0 0 size(mam,2)/2 size(mam,1)/2]);
imagesc(mam); colormap(gray(256));  axis image; hold on; 
plot(segmentation.breast_border(segmentation.breast_air,1),...
     segmentation.breast_border(segmentation.breast_air,2), 'g--', 'LineWidth', 1.5);
set(gca, 'xticklabel', [], 'yticklabel', []);

print('-dtiff', '-noui', '-painters', f3, '-r300', 'C:\isbe\thesis\figures\10\intro_mammogram_cc.tif');
%%
mam = u_load('C:\isbe\dev\segmentation\mammograms\o04_031LML.mat');
segmentation = u_load('C:\isbe\dev\segmentation\breast_borders\o04_031LML_segmentation.mat');

f5 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [100 100 size(mam,2)/2 size(mam,1)/2],...
    'PaperPositionMode','auto');
axes('units', 'pixels', 'position', [0 0 size(mam,2)/2 size(mam,1)/2]);
imagesc(mam); colormap(gray(256));  axis image; hold on; 
plot(segmentation.breast_border(segmentation.breast_air(1:end-10),1),...
     segmentation.breast_border(segmentation.breast_air(1:end-10),2), 'g--', 'LineWidth', 1.5);
set(gca, 'xticklabel', [], 'yticklabel', []);

print('-dtiff', '-noui', '-painters', f5, '-r300', 'C:\isbe\thesis\figures\10\intro_mammogram_mlo.tif');
%%
%--------------------------------------------------------------------------
% Show pectoral fitting result
%--------------------------------------------------------------------------
% Choose o04_008LML_segmentation.mat
mam = u_load('C:\isbe\dev\segmentation\mammograms\o04_008LML.mat');
[segmentation normal] = segment_breast('image', mam, 'mlo', 1, 'right', 0, 'plot', 1);
f6 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [100 100 size(mam,2)/2 size(mam,1)/2],...
    'PaperPositionMode','auto');
axes('units', 'pixels', 'position', [0 0 size(mam,2)/2 size(mam,1)/2]);
imagesc(mam); colormap(gray(256));  axis image; hold on; 
plot(segmentation.breast_border(segmentation.breast_air,1), segmentation.breast_border(segmentation.breast_air,2), 'r', 'LineWidth', 1.5);
plot(segmentation.pectoral_edge(:,1), segmentation.pectoral_edge(:,2), 'r', 'LineWidth', 1.5);
plot(250, 500, 'gx', 'MarkerSize', 10);
plot([1 segmentation.breast_border(segmentation.breast_air(1),1)],...
     [500 500], 'b:', 'LineWidth', 2);
 plot([segmentation.breast_border(segmentation.breast_air(1),1) segmentation.breast_border(segmentation.breast_air(1),1)],...
     [1 500], 'b:', 'LineWidth', 2);
print('-dtiff', '-noui', '-painters', f6, '-r300', 'C:\isbe\thesis\figures\10\segmented_pectoral.tif');
%%
%--------------------------------------------------------------------------
% Show how we get a standard MLO shape
%--------------------------------------------------------------------------
% Use 014LML
ml_list = dir('C:\isbe\dev\segmentation\breast_borders\*014LML*.mat');
[ml_shapes, ml_areas] = get_mlo_shape(ml_list, 50, 1);
%%
mam = u_load('C:\isbe\dev\segmentation\mammograms\o04_014LML.mat');
segmentation = u_load('C:\isbe\dev\segmentation\breast_borders\o04_014LML_segmentation.mat');

breast_border = segmentation.breast_border;
breast_air = segmentation.breast_air;
pectoral_edge = segmentation.pectoral_edge;
rows = segmentation.size(1);
cols = segmentation.size(1);

bw1 = ~poly2mask([1; pectoral_edge(:,1)], [1; pectoral_edge(:,2)], rows, cols);
bw2 = poly2mask(breast_border(:,1), breast_border(:,2), rows, cols);
%figure; imagesc(bw1 & bw2); axis image; colormap gray
[yy xx] = find(bw1 & bw2);
%breast_centre = [mean(xx) mean(yy)];
    
x1 = pectoral_edge(1,1);
x2 = pectoral_edge(2,1);
y1 = pectoral_edge(1,2);
y2 = pectoral_edge(2,2);

m = (y1 - y2) / (x1 - x2);
c = y1 - m*x1;

x3 = 1;
y3 = m*x3 + c;

x4 = mean(xx);
y4 = mean(yy);

mn = -1 / m;
cn = y4 - mn*x4;

x5 = (cn - c) / (m - mn);
y5 = mn*x5 + cn;

apos = [40 40 max(breast_border(:,1)) + 20, max(breast_border(:,2)) + 20];
xlims = [-10 max(breast_border(:,1)) + 10];
ylims = [-10 max(breast_border(:,2)) + 10];
fpos = [0 0 apos(3) + 80 apos(4) + 80];

f1 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', fpos,...
    'PaperPositionMode','auto');
axes(...
    'Units', 'pixels',...
    'position', apos,...
    'xlim', xlims,...
    'ylim', ylims);

axis ij equal off; hold on;
plot(breast_border(:,1), breast_border(:,2), 'b--');
plot(breast_border(breast_air,1),...
     breast_border(breast_air,2), 'b', 'LineWidth', 2);
plot([x1 x3], [y1 y3], 'r', 'LineWidth', 2);
plot(x4, y4, 'rx', 'MarkerSize', 10);

print('-dtiff', '-noui', '-painters', f1, '-r300', 'C:\isbe\thesis\figures\10\mlo_shape_unrotated.tif');
%
theta = atan2(diff(pectoral_edge(:,2)), diff(pectoral_edge(:,1)));
    
%Compute the rotation needed so the pectoral edge lies vertically in
%the image
rot_mat = [cos(pi/2-theta) sin(pi/2-theta); -sin(pi/2-theta) cos(pi/2-theta)];
    
%Rotate the pectoral edge and compute the translation required to move
%this to the y-axis
pectoral_edge_rot = pectoral_edge * rot_mat;
t = pectoral_edge_rot(1,:);
    
pectoral_edge_rot = pectoral_edge_rot - [t;t];
breast_centre_rot = [x4 y4] * rot_mat - t;

%Rotate and translate the breast border
breast_border_rot = breast_border * rot_mat - repmat(t, size(breast_border,1),1);

xx = breast_border_rot(breast_air,1);
yy = breast_border_rot(breast_air,2);

xx = [xx(2); xx(1); xx; xx(end); xx(end-1);]; %#ok
xt = xx(2:end-3) - xx(4:end-1);
xtt = xx(1:end-4) - xx(5:end);
yy = [yy(2); yy(1); yy; yy(end); yy(end-1);]; %#ok
yt = yy(2:end-3) - yy(4:end-1);
ytt = yy(1:end-4) - yy(5:end);

xx([1:2, end-1:end],:) = [];
yy([1:2, end-1:end],:) = [];
xyn = [-yt xt] ./ [sqrt(xt.^2 + yt.^2) sqrt(xt.^2 + yt.^2)];
xd = [0; diff(xx)];
yd = [0; diff(yy)];
D = sum([xd yd] .* xyn, 2);
%
D_smooth = conv2(D,fspecial('gaussian', [30 1], 5), 'same');

halfway = round(length(D)/2);
[min_up min_up_idx] = min(D_smooth(1:halfway));
[min_do min_do_idx] = min(D_smooth(halfway+1:end));
min_do_idx = min_do_idx + halfway;

keep = true(length(breast_air),1);
keep(1:min_up_idx-1) = 0;
keep(min_do_idx+1:end) = 0;
    
x6 = xx(min_up_idx);
y6 = yy(min_up_idx);
x7 = xx(min_do_idx);
y7 = yy(min_do_idx);

breast_air_pts = breast_border_rot(breast_air(keep),:);

%Find the point on breast air closest to the intersection of a line
%running normal from the chest wall through the centroid
%The line by construction is horizontal, so just need to find the
%breast air point with closest matching y to centre_y
[dummy idx] = min((breast_air_pts(:,2)-breast_centre_rot(2)).^2);

%Split breast border into 3 segemnts - upper breast air, lower breast
%air and chest wall and space 50 points equally along each
segment1 = breast_air_pts(1:idx,:);
segment2 = breast_air_pts(idx:end,:);
n_pts = 20;

diff1 = diff(segment1);
cum_dist1 = [0; cumsum(sqrt(sum(diff1.^2,2)))];
segment1 = interp1(cum_dist1, segment1, linspace(0, cum_dist1(end), n_pts), 'linear');

diff2 = diff(segment2);
cum_dist2 = [0; cumsum(sqrt(sum(diff2.^2,2)))];
segment2 = interp1(cum_dist2, segment2, linspace(0, cum_dist2(end), n_pts), 'linear');
%%
segment3 = [ones(10,1) linspace(segment1(1,2), segment2(1,2), 10)'];
segment4 = [ones(10,1) linspace(segment2(1,2), segment2(end,2), 10)'];    

%
apos = [40 40 max(breast_border_rot(:,1)) + 20, max(breast_border_rot(:,2)) - min(breast_border_rot(:,2)) + 20];
xlims = [-10 max(breast_border_rot(:,1)) + 10];
ylims = [min(breast_border_rot(:,2)) - 10 max(breast_border_rot(:,2)) + 10];
fpos = [0 0 apos(3) + 80 apos(4) + 80];

f2 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', fpos,...
    'PaperPositionMode','auto');
axes(...
    'Units', 'pixels',...
    'position', apos,...
    'xlim', xlims,...
    'ylim', ylims);

axis ij equal off; hold on;
%plot(breast_border_rot(:,1), breast_border_rot(:,2), 'b--');
plot(breast_border_rot(breast_air,1), breast_border_rot(breast_air,2), 'b', 'LineWidth', 2);
plot(breast_centre_rot(1), breast_centre_rot(2)+4, 'rx', 'MarkerSize', 10);
plot([1 1], [breast_border_rot(breast_air(1),2) breast_border_rot(breast_air(end),2)], 'r', 'LineWidth', 2);
plot(segment1(:,1), segment1(:,2), 'ro', 'MarkerSize', 10);
plot(segment2(:,1), segment2(:,2), 'ro', 'MarkerSize', 10);
plot(segment3(:,1), segment3(:,2), 'ro', 'MarkerSize', 10);
plot(segment4(:,1), segment4(:,2), 'ro', 'MarkerSize', 10);
plot(segment2(1,1), segment2(1,2), 'go', 'MarkerSize', 10);
plot(segment4(1,1), segment4(1,2), 'go', 'MarkerSize', 10);
plot(x6, y6, 'go', 'MarkerSize', 10);
plot(1, y6, 'go', 'MarkerSize', 10);
plot(x7, y7, 'go', 'MarkerSize', 10);
plot(1, y7, 'go', 'MarkerSize', 10);


print('-dtiff', '-noui', '-painters', f2, '-r300', 'C:\isbe\thesis\figures\10\mlo_shape_rotated.tif');

%%
%--------------------------------------------------------------------------
% Compute probability densities at the skin-air border
%--------------------------------------------------------------------------
%load C:\isbe\dev\location\ML_location_model
load C:\isbe\dev\location\CC_location_model %repeat for CC

%Get start and end (x,y)-coordinates to bound the mean breast shape
xs = floor(min(location_model.mean_shape(:,1)));
xe = ceil(max(location_model.mean_shape(:,1)));
ys = floor(min(location_model.mean_shape(:,2)));
ye = ceil(max(location_model.mean_shape(:,2)));

mask_size(1) = ye - ys + 1;
mask_size(2) = xe - xs + 1;

%Build a BW mask of the shape
[bw] = poly2mask(...
    location_model.mean_shape(:,1) - xs + 1,...
    location_model.mean_shape(:,2) - ys + 1,...
    mask_size(1),...
    mask_size(2));

[in_idx] = find(bw);
D_map_a = zeros(mask_size);
D_map_a(in_idx) = location_model.D_a;


bw2 = bw;
bw2(:,1) = 1;
bw2 = imerode(bw2, strel('disk', 100));

skin_region = bw & ~bw2;

D_map_mask = D_map_a;
D_map_mask(skin_region) = 0;
D_map_mask = D_map_mask / sum(D_map_mask(:));

figure;    
subplot(1,2,1); imagesc(skin_region); axis image; colormap(jet(256));
subplot(1,2,2); imagesc(D_map_mask); axis image; colormap(jet(256));

hold on;
plot(location_model.mean_shape(:,1) - xs + 1, location_model.mean_shape(:,2) - ys + 1,'y');

display(sum(D_map_a(:)));
display(sum(D_map_a(skin_region)));
%%
%
%%
%--------------------------------------------------------------------------
% Show MLO and CC mammogram with breast border
%--------------------------------------------------------------------------

%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-------- Some useful other auxilliary code -------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
mlo_list = dir('C:\isbe\dev\segmentation\breast_borders\*ML*.mat');
mlo_list2 = dir('C:\isbe\dev\segmentation\mammograms\*ML*.mat');
for ii = 12
    segmentation = u_load(['C:\isbe\dev\segmentation\breast_borders\', mlo_list(ii).name]);
    mam = u_load(['C:\isbe\dev\segmentation\mammograms\', mlo_list2(ii).name]);
    figure; imagesc(mam); axis image; colormap(gray(256)); hold on;
    plot(segmentation.breast_border(:,1), segmentation.breast_border(:,2), 'g', 'LineWidth', 1.5);
    plot(segmentation.breast_border(segmentation.breast_air,1), segmentation.breast_border(segmentation.breast_air,2), 'r', 'LineWidth', 1.5);
    plot(segmentation.pectoral_edge(:,1), segmentation.pectoral_edge(:,2), 'y', 'LineWidth', 1.5);
end