%script for Janine

%1 Sample from density distribution (a raleigh distribution with sigma = 25
%density_sample = sort(25*sqrt(-2*log(rand(15,1))));
%save D:\isbe\dev\blobs\density_sample.mat density_sample
load D:\isbe\dev\blobs\density_sample.mat

%%
%Display PDF of distribution
x = 0:100;
pdf = (x/25^2) .* exp(-0.5*(x / 25).^2);
figure; plot(x, pdf, 'linewidth', 2);
xlabel('X', 'fontsize', 24);
ylabel('f(x)', 'fontsize', 24);
title('PDF of Raleigh distribution, \sigma = 25', 'fontsize', 24);
set(gca, 'fontsize', 16);
%%

load C:\isbe\mammograms\new_CAD\bmP_2004_half\o04_013LCC.mat
load D:\isbe\dev\segmentation\breast_borders\o04_013LCC_segmentation.mat
%%
mammogram = imresize(mammogram, segmentation.size, 'bilinear');

load D:\isbe\dev\mass_model\models\model_w500_50K.mat

len = 500;
dists = cumsum([0; sqrt(sum(diff(segmentation.breast_border).^2, 2))]);
%
breast_shape = interp1(dists, segmentation.breast_border, linspace(0, dists(end), len), 'linear');
breast_mask = roipoly(mammogram, breast_shape(:,1), breast_shape(:,2));

[dummy nipple_idx] = max(breast_shape(:,1));
start_point = breast_shape(nipple_idx,:) - [25 0];
breast_area = polyarea(breast_shape(:,1), breast_shape(:,2));

[d, aligned_shape,transform] = mb_procrustes(mass_model.mean_target,breast_shape);
keep = 10;
b_shape = mass_model.P_shape' * (aligned_shape(:)'-mass_model.mean_shape)';
%
for dd = 1:15

    b_shape(keep+1:end) = (randn(length(mass_model.L_shape)-keep, 1) .* sqrt(mass_model.L_shape(keep+1:end)));

    blob_shape = reshape(mass_model.P_shape*b_shape + mass_model.mean_shape', [], 2);

%     figure; hold on; axis equal;
%     plot(blob_shape(:,1), blob_shape(:,2), 'r');
%     plot(aligned_shape(:,1), aligned_shape(:,2), 'g');

    blob_area = polyarea(blob_shape(:,1), blob_shape(:,2));
    
    scale_factor = sqrt(density_sample(dd)*breast_area/(blob_area*100)); %#ok

    blob_shape = scale_factor * blob_shape* inv(transform.T);
    blob_shape(:,1) = blob_shape(:,1) - blob_shape(nipple_idx,1) + start_point(1);
    blob_shape(:,2) = blob_shape(:,2) - blob_shape(nipple_idx,2) + start_point(2);
    
    
    blob_areas(dd) = polyarea(blob_shape(:,1), blob_shape(:,2));  %#ok
    blob_mask = roipoly(mammogram, blob_shape(:,1), blob_shape(:,2));
    
    blob_mammogram = mammogram;
    blob_mammogram(breast_mask) = mean(mammogram(breast_mask));
    blob_mammogram(blob_mask) = blob_mammogram(blob_mask)*1.4;
    
    %write_im_from_colormap(blob_mammogram, ['D:\isbe\dev\blobs\blobs\blob' zerostr(dd,2) '.bmp'], gray(256));
    
    figure; 
    subplot(1,2,1); imagesc(mammogram); axis image; colormap(gray(256)); hold on;

    plot(...
        [1 segmentation.size(2) segmentation.size(2) 1 1],...
        [1 1 segmentation.size(1) segmentation.size(1) 1]);
    plot(segmentation.breast_border(:,1), segmentation.breast_border(:,2));
    plot(blob_shape(:,1), blob_shape(:,2), 'r');
     
    subplot(1,2,2); imagesc(blob_mammogram); axis image; colormap(gray(256));
end
%%
%--------------------------------------------------------------------------
% Need to workout what size blob we should use when varying breast size
seg_list = dir('D:\isbe\dev\segmentation\breast_borders\*CC*segmentation.mat');
breast_areas = zeros(length(seg_list),1);
for ii = 1:length(seg_list)
    load(['D:\isbe\dev\segmentation\breast_borders\' seg_list(ii).name]);
    breast_areas(ii) = polyarea(segmentation.breast_border(:,1), segmentation.breast_border(:,2));
end
%%
figure; hold all; ylim([0 1]);
for dd = 1:15; 
    plot(1:length(seg_list), blob_areas(dd) ./ sort(breast_areas), '-'); 
end
%%
for bb = 1:15
    display(['Blob ', num2str(bb)]);
    for dd = 1:15
        [dummy breast_idx(dd)] = min(abs(density_sample(dd) - 100*blob_areas(bb) ./ breast_areas)); %#ok
        display(['min dist = ', num2str(dummy)]);
    end
    display('-*-');
end
% outcome -> choose blob 8
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Now do blobs with varying breast size
load D:\isbe\dev\mass_model\models\model_w500_50K.mat
for dd = 1:15
    [dummy breast_idx(dd)] = min(abs(density_sample(dd) - 100*blob_areas(8) ./ breast_areas)); %#ok
    new_density_sample(dd) = 100*blob_areas(8) ./ breast_areas(breast_idx(dd)); %#ok
    
    load(['D:\isbe\mammograms\new_CAD\bmP_2004_half\' seg_list(breast_idx(dd)).name(1:10) '.mat']);
    load(['D:\isbe\dev\segmentation\breast_borders\' seg_list(breast_idx(dd)).name]);
    
    mammogram = imresize(mammogram, segmentation.size, 'bilinear');
    
    len = 500;
    dists = cumsum([0; sqrt(sum(diff(segmentation.breast_border).^2, 2))]);
    %
    breast_shape = interp1(dists, segmentation.breast_border, linspace(0, dists(end), len), 'linear');
    breast_mask = roipoly(mammogram, breast_shape(:,1), breast_shape(:,2));

    if isempty(strfind(seg_list(breast_idx(dd)).name(1:10), 'R'))
        [dummy nipple_idx] = max(breast_shape(:,1));
        start_point = breast_shape(nipple_idx,:) - [30 0];
    else
        [dummy nipple_idx] = min(breast_shape(:,1));
        start_point = breast_shape(nipple_idx,:) + [30 0];
    end
    
    breast_area = polyarea(breast_shape(:,1), breast_shape(:,2));

    [d, aligned_shape,transform] = mb_procrustes(mass_model.mean_target,breast_shape);
    keep = 10;
    b_shape = mass_model.P_shape' * (aligned_shape(:)'-mass_model.mean_shape)';
    
    b_shape(keep+1:end) = (randn(length(mass_model.L_shape)-keep, 1) .* sqrt(mass_model.L_shape(keep+1:end)));

    blob_shape = reshape(mass_model.P_shape*b_shape + mass_model.mean_shape', [], 2);

    blob_area = polyarea(blob_shape(:,1), blob_shape(:,2));
    
    scale_factor = sqrt(new_density_sample(dd)*breast_area/(blob_area*100));

    blob_shape = scale_factor * blob_shape* inv(transform.T);
    blob_shape(:,1) = blob_shape(:,1) - blob_shape(nipple_idx,1) + start_point(1);
    blob_shape(:,2) = blob_shape(:,2) - blob_shape(nipple_idx,2) + start_point(2);
    display(['blob area = ' num2str(polyarea(blob_shape(:,1), blob_shape(:,2)))]);
    
    blob_mask = roipoly(mammogram, blob_shape(:,1), blob_shape(:,2));
    
    blob_mammogram = mammogram;
    blob_mammogram(breast_mask) = mean(mammogram(breast_mask));
    blob_mammogram(blob_mask) = blob_mammogram(blob_mask)*1.4;
    
    write_im_from_colormap(blob_mammogram, ['D:\isbe\dev\blobs\breasts\mammo' zerostr(dd,2) '.bmp'], gray(256));
    
    figure; 
    subplot(1,2,1); imagesc(mammogram); axis image; colormap(gray(256)); hold on;

    plot(...
        [1 segmentation.size(2) segmentation.size(2) 1 1],...
        [1 1 segmentation.size(1) segmentation.size(1) 1]);
    plot(segmentation.breast_border(:,1), segmentation.breast_border(:,2));
    plot(blob_shape(:,1), blob_shape(:,2), 'r');
     
    subplot(1,2,2); imagesc(blob_mammogram); axis image; colormap(gray(256));
end
save D:\isbe\dev\blobs\density_sample.mat density_sample new_density_sample
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Now repeat the varying blobs with the new density sample
load D:\isbe\mammograms\new_CAD\bmP_2004_half\o04_013LCC.mat
load D:\isbe\dev\segmentation\breast_borders\o04_013LCC_segmentation.mat
mammogram = imresize(mammogram, segmentation.size, 'bilinear');

load D:\isbe\dev\mass_model\models\model_w500_50K.mat
len = 500;
dists = cumsum([0; sqrt(sum(diff(segmentation.breast_border).^2, 2))]);
%
breast_shape = interp1(dists, segmentation.breast_border, linspace(0, dists(end), len), 'linear');
breast_mask = roipoly(mammogram, breast_shape(:,1), breast_shape(:,2));

[dummy nipple_idx] = max(breast_shape(:,1));
start_point = breast_shape(nipple_idx,:) - [25 0];
breast_area = polyarea(breast_shape(:,1), breast_shape(:,2));

[d, aligned_shape,transform] = mb_procrustes(mass_model.mean_target,breast_shape);
keep = 10;
b_shape = mass_model.P_shape' * (aligned_shape(:)'-mass_model.mean_shape)';
%
for dd = 1:15

    b_shape(keep+1:end) = (randn(length(mass_model.L_shape)-keep, 1) .* sqrt(mass_model.L_shape(keep+1:end)));

    blob_shape = reshape(mass_model.P_shape*b_shape + mass_model.mean_shape', [], 2);

    blob_area = polyarea(blob_shape(:,1), blob_shape(:,2));
    
    scale_factor = sqrt(new_density_sample(dd)*breast_area/(blob_area*100));

    blob_shape = scale_factor * blob_shape* inv(transform.T);
    blob_shape(:,1) = blob_shape(:,1) - blob_shape(nipple_idx,1) + start_point(1);
    blob_shape(:,2) = blob_shape(:,2) - blob_shape(nipple_idx,2) + start_point(2);
    
    blob_mask = roipoly(mammogram, blob_shape(:,1), blob_shape(:,2));
    
    blob_mammogram = mammogram;
    blob_mammogram(breast_mask) = mean(mammogram(breast_mask));
    blob_mammogram(blob_mask) = blob_mammogram(blob_mask)*1.4;
    
    write_im_from_colormap(blob_mammogram, ['D:\isbe\dev\blobs\blobs\blob' zerostr(dd,2) '.bmp'], gray(256));
    
    figure; 
    subplot(1,2,1); imagesc(mammogram); axis image; colormap(gray(256)); hold on;

    plot(...
        [1 segmentation.size(2) segmentation.size(2) 1 1],...
        [1 1 segmentation.size(1) segmentation.size(1) 1]);
    plot(segmentation.breast_border(:,1), segmentation.breast_border(:,2));
    plot(blob_shape(:,1), blob_shape(:,2), 'r');
     
    subplot(1,2,2); imagesc(blob_mammogram); axis image; colormap(gray(256));
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Create streaky mammograms
mkdir D:\isbe\dev\blobs\lines\
load D:\isbe\mammograms\new_CAD\bmP_2004_half\o04_013LCC.mat
load D:\isbe\dev\segmentation\breast_borders\o04_013LCC_segmentation.mat
mammogram = imresize(mammogram, segmentation.size, 'bilinear');

len = 500;
dists = cumsum([0; sqrt(sum(diff(segmentation.breast_border).^2, 2))]);
%
breast_shape = interp1(dists, segmentation.breast_border, linspace(0, dists(end), len), 'linear');
breast_mask = roipoly(mammogram, breast_shape(:,1), breast_shape(:,2));

[dummy nipple_idx] = max(breast_shape(:,1));
start_point = round(breast_shape(nipple_idx,:) - [1 0]);
breast_area = polyarea(breast_shape(:,1), breast_shape(:,2));

for dd = 1:15
    line_sum = breast_area * density_sample(dd) / 100;
    line_mask = false(size(breast_mask));
    while sum(line_mask(:)) < line_sum
        line_mask = generate_line_tree(breast_mask, round(start_point), 3*pi*rand/4 - 3*pi/8, line_mask, line_sum);
    end
    
    line_mammogram = mammogram;
    line_mammogram(breast_mask) = mean(mammogram(breast_mask));
    line_mammogram(line_mask) = line_mammogram(line_mask)*1.4;
    figure; 
    subplot(1,2,1); imagesc(line_mask); axis image; colormap(gray(256));
    subplot(1,2,2); imagesc(line_mammogram); axis image; colormap(gray(256));
    
    write_im_from_colormap(line_mammogram, ['D:\isbe\dev\blobs\lines\line' zerostr(dd,2) '.bmp'], gray(256));
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Create multiblobs
mkdir D:\isbe\dev\blobs\multi\
load D:\isbe\dev\blobs\density_sample.mat
load C:\isbe\mammograms\new_CAD\bmP_2004_half\o04_013LCC.mat
load D:\isbe\dev\segmentation\breast_borders\o04_013LCC_segmentation.mat
mammogram = imresize(mammogram, segmentation.size, 'bilinear');

load D:\isbe\dev\mass_model\models\model_w500_50K.mat
len = 500;
dists = cumsum([0; sqrt(sum(diff(segmentation.breast_border).^2, 2))]);
%
breast_shape = interp1(dists, segmentation.breast_border, linspace(0, dists(end), len), 'linear');
breast_mask = roipoly(mammogram, breast_shape(:,1), breast_shape(:,2));
%breast_mask(:,1:2) = true;
breast_area = polyarea(breast_shape(:,1), breast_shape(:,2));

%--------------------------------------------------------------------------
mean_shape  = mass_model.mean_shape;
P_shape     = mass_model.P_shape;
L_shape     = mass_model.L_shape;
mean_scale  = mass_model.mean_scale;
P_scale     = mass_model.P_scale;
mean_tex    = mass_model.mean_tex;
P_tex       = mass_model.P_tex;
L_tex       = mass_model.L_tex;

P_com         = mass_model.P_com;
L_com         = mass_model.L_com;

W_shape     = mass_model.W_shape;
W_tex       = mass_model.W_tex;
W_scale     = mass_model.W_scale;

mean_shape_pl = mass_model.mean_shape_pl;
size_shape_vec = length(mean_shape) / 2;

k_shape     = length(L_shape);
k_tex       = length(L_tex);

scaled_mean_shape = mass_model.scaled_mean_shape;
  
%Define source points for TPS - as row vectors
s_x = scaled_mean_shape(1:size_shape_vec);
s_y = scaled_mean_shape(size_shape_vec+1:end);

%Define points to be interpolated by TPS - as row vectors
i_x = mean_shape_pl(:,1)';
i_y = mean_shape_pl(:,2)';
%--------------------------------------------------------------------------

%
for dd = 1:14

    blob_sum = breast_area * density_sample(dd) / 100;
    blob_mask = false(size(breast_mask));
    blob_mammogram = zeros(size(breast_mask));
    
    while sum(blob_mask(:)) < blob_sum
        
        %Generate random rotation of shape
        theta = 2*pi*rand;
        theta_m = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        % Compute new shape vector, texture vector and scale
        blob_com = (randn(length(L_com), 1) .* sqrt(L_com));

        Q_shape = P_com(1:k_shape,:); 
        Q_tex = P_com(k_shape+1:k_shape + k_tex,:);
        Q_scale = P_com(end, :);

        %Use sampled parameters to generate texture
        blob_tex = mean_tex + (P_tex*inv(W_tex)*Q_tex*blob_com)';
        blob_tex = blob_tex * 1.2;
        
        %Use sampled parameters to generate shape
        blob_shape = mean_shape + (P_shape*inv(W_shape)*Q_shape*blob_com)';
        blob_shape = rand*reshape(blob_shape, [], 2)* theta_m;

        %Workout where we can stick this blob
        blob_size = round(max(max(blob_shape) - min(blob_shape)));
        blob_map = imerode(breast_mask, strel('disk', blob_size)) & ~blob_mask;
        [pts_y pts_x] = find(blob_map);
        
        %Are there any viable pts...?
        if ~isempty(pts_x)
            %If yes, select one at random
            centre_idx = ceil(length(pts_x)*rand);

            %Centre the new blow shape here
            blob_shape(:,1) = blob_shape(:,1) - mean(blob_shape(:,1)) + pts_x(centre_idx);
            blob_shape(:,2) = blob_shape(:,2) - mean(blob_shape(:,2)) + pts_y(centre_idx);
            
            %Get a mask of the new blob, and compute the indices of pixels
            %belonging to the blob
            blob_bw = roipoly(blob_mask, blob_shape(:,1), blob_shape(:,2));
            blob_idx = find(blob_bw);
            [blob_y blob_x] = ind2sub(size(blob_mask), blob_idx);

            % Compute TPS warp to texture from mean to new blob
            %%%%

            %Define displacement to target points
            z_x = blob_shape(:,1)';
            z_y = blob_shape(:,2)';

            %Compute displacement of interpolated points        
            T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
                'transform', 'spline');
            [f_xy] = geom_transformpoints([i_x; i_y], T);

            %Interpolate uneven grid of sampled texture points to pixel
            %lattice
            blob_shape_tex = griddata(f_xy(1,:), f_xy(2,:), blob_tex,...
                blob_x, blob_y);
            
            %Remove any unrealistic values
            blob_shape_tex(isnan(blob_shape_tex)) = 0;
            blob_shape_tex(blob_shape_tex < 0) = 0;
            
            %Add the new blob to the blob mask and blob mammogram
            blob_mask = blob_mask | blob_bw;
            blob_mammogram(blob_bw) =  blob_mammogram(blob_bw)+blob_shape_tex;
            
            
%             figure; 
%             subplot(1,2,1); imagesc(blob_map); axis image; colormap(gray(256)); hold on;
%             subplot(1,2,2); 
%             imagesc(blob_mask); axis image; colormap(gray(256)); hold on;
%             plot(segmentation.breast_border(:,1), segmentation.breast_border(:,2));

        end
    
    end
    
    %blob_mammogram = mammogram;
    %blob_mammogram(breast_mask) = mean(mammogram(breast_mask));
    %blob_mammogram(blob_mask) = blob_mammogram(blob_mask)*1.4;
    
    %write_im_from_colormap(blob_mammogram, ['D:\isbe\dev\blobs\multi\blob' zerostr(dd,2) '.bmp'], gray(256));
    
    figure; 
    subplot(1,2,1); imagesc(blob_mask); axis image; colormap(gray(256)); hold on;
    subplot(1,2,2); imagesc(blob_mammogram); axis image; colormap(gray(256));
end