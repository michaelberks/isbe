
%load image
nailfold = imread('C:\isbe\nailfold\anonymous_bmap_files_oct\OCT013V1LD4X3LrgMosaic.bmp');
nailfold(:,2001:end) = [];
mask = nailfold < 255;
[r c] = size(nailfold);

%Compute line strength and orientation of full size image
[line_strength, line_ori] = gaussian_2nd_derivative_line(nailfold, 4);

%Compute line strength and orientation of half size images then upscale to
%full size
[line_strength2, line_ori2] =...
    gaussian_2nd_derivative_line(imresize(nailfold, 0.5, 'bilinear'), 2);
%
%Apply NMS to full size maps and rescaled half-size maps
vessel_nms = mb_non_maximal_supp(line_strength, line_ori);
vessel_nms2 = mb_non_maximal_supp(imresize(line_strength2, [r c], 'bilinear'), imresize(line_ori2, [r c], 'nearest'));
%%
%Get the non-zero points from both
[y_and x_and] = find(vessel_nms & vessel_nms2 & mask);
[y x] = find(vessel_nms & ~vessel_nms2 & mask);
[y2 x2] = find(vessel_nms2 & ~vessel_nms & mask);

%Display the nailfold and plot the 2 sets of points
figure; imagesc(nailfold); axis image; colormap(gray(256)); hold on;
plot(x_and, y_and, 'g.', 'markersize', 2);
plot(x, y, 'r.', 'markersize', 2);
plot(x2, y2, 'b.', 'markersize', 2);
%%
mkdir C:\isbe\nailfold\anonymous_bmap_files_enhanced
n_list = dir('C:\isbe\nailfold\anonymous_bmap_files_oct\*.bmp');
for ii = 1:length(n_list)
    
    nailfold = imread(['C:\isbe\nailfold\anonymous_bmap_files_oct\' n_list(ii).name]);
    if size(nailfold,3)==3
        nailfold = nailfold(:,:,1);
    end
    
    nailfold = double(nailfold);
    mask = nailfold < 255;
    nailfold(~mask) = 0;
    win_size = 63;

    local_n = imfilter(double(mask), ones(win_size));
    local_sum = imfilter(nailfold, ones(win_size));
    local_sum2 = imfilter(nailfold.^2, ones(win_size));

    local_mean = local_sum ./ local_n;
    local_mean2 = local_sum2 ./ local_n;
    local_mean(~mask) = 0;
    local_mean2(~mask) = 0;

    local_std = sqrt(local_mean2 - local_mean.^2);
    local_norm = (nailfold - local_mean) ./ local_std;

    nailfold = uint8(255*(local_norm - min(local_norm(:))) / (max(local_norm(:)) - min(local_norm(:))));
    imwrite(nailfold, ['C:\isbe\nailfold\anonymous_bmap_files_enhanced\' n_list(ii).name]);
    
end
%%
mkdir C:\isbe\nailfold\anonymous_bmap_files\line\g2d
mkdir C:\isbe\nailfold\anonymous_bmap_files\ori\g2d
n_list = dir('C:\isbe\nailfold\anonymous_bmap_files_oct\*.bmp');
for ii = 1:length(n_list)
    
    nailfold = imread(['C:\isbe\nailfold\anonymous_bmap_files_oct\' n_list(ii).name]);
    if size(nailfold,3)==3
        nailfold = nailfold(:,:,1);
    end
    
    nailfold = double(nailfold);
    
    %Compute line strength and orientation of full size image
    [line_strength, line_ori] = gaussian_2nd_derivative_line(nailfold, [1 2 4 8]);
    
    %Save both line and orientation maps
    save(['C:\isbe\nailfold\anonymous_bmap_files\line\g2d' n_list(ii).name(1:end-3) '_line.mat'], 'line_strength');
    save(['C:\isbe\nailfold\anonymous_bmap_files\ori\g2d' n_list(ii).name(1:end-3) '_ori.mat'], 'line_ori');
    clear nailfold line_strength line_ori;
end
%%
n_list = dir('C:\isbe\nailfold\anonymous_bmap_files_enhanced\*.bmp');
for ii = 37:length(n_list)
    
    nailfold = imread(['C:\isbe\nailfold\anonymous_bmap_files_enhanced\' n_list(ii).name]);
    figure; imagesc(nailfold); axis image; colormap(gray(256)); hold on;
    clear nailfold;
    
    %Save both line and orientation maps
    line_strength = u_load(['C:\isbe\nailfold\anonymous_bmap_files\line\g2d' n_list(ii).name(1:end-3) '_line.mat']);
    line_ori = u_load(['C:\isbe\nailfold\anonymous_bmap_files\ori\g2d' n_list(ii).name(1:end-3) '_ori.mat']);
    
    col = floor(size(line_strength,2)/2) ;
    
    %Apply NMS to full size maps and rescaled half-size maps
    vessel_nms1 = mb_non_maximal_supp(line_strength(:,1:col), line_ori(:,1:col));

    %Get the non-zero points from both
    [y1 x1] = find(bwareaopen(vessel_nms1 > 0, 20)); clear vessel_nms1;
    
    %Apply NMS to full size maps and rescaled half-size maps
    vessel_nms2 = mb_non_maximal_supp(line_strength(:,col+1:end), line_ori(:,col+1:end));

    %Get the non-zero points from both
    [y2 x2] = find(bwareaopen(vessel_nms2 > 0, 20)); clear vessel_nms2;
    
    plot(x1, y1, 'g.', 'markersize', 2);
    plot(x2+col, y2, 'g.', 'markersize', 2);


    clear line_strength line_ori;
    pack;
end

%%    
nailfold = imread('C:\isbe\nailfold\ncm\n Tonia MooreV1LD4X3LrgMosaic.bmp');
nailfold(:,2001:end) = [];
mask = nailfold < 255;
[r c] = size(nailfold);

%Compute line strength and orientation of full size image
%[line_strength, line_ori] = gaussian_2nd_derivative_line(nailfold, 4);
load C:\isbe\nailfold\annotations\'n Tonia MooreV1LD4X3LrgMosaic_vessels.mat'
%%
%Compute line strength and orientation of half size images then upscale to
%full size
for scale = 2:8
    [line_strength, line_ori] =...
        gaussian_2nd_derivative_line(imresize(nailfold, 1 / scale, 'bilinear'), 4 / scale);
    
    vessel_nms = mb_non_maximal_supp(...
        imresize(line_strength, [r c], 'bilinear'),...
        imresize(line_ori, [r c], 'nearest'));
    
    %vessel_nms = mb_non_maximal_supp(line_strength, line_ori);
    
    [y_vessels x_vessels] = find(vessel_nms);
    %vessels_scale = vessels;
    
    %snap the [x y] point to nearest potential vessel
    %figure; imagesc(nailfold); axis image; colormap(gray(256)); hold on;
    
    display(['Scale ' num2str(scale)]);
    for ii = 1:length(vessels)
        %plot(vessels{ii}(:,1), vessels{ii}(:,2), 'y');
        %plot(vessels{ii}(:,1), vessels{ii}(:,2), 'gx');
        missing_pts = 0;
        num_pts = size(vessels{ii},1);
        for jj = 1:num_pts
            x = vessels{ii}(jj,1);
            y = vessels{ii}(jj,2);
            
            neighbours = (abs(x_vessels - x) < 5) & (abs(y_vessels - y) < 5);

            if any(neighbours)
                xn = x_vessels(neighbours);
                yn = y_vessels(neighbours);

                dists = (xn - x).^2 + (yn - y).^2;
                [dummy min_idx] = min(dists);
                xn = xn(min_idx);
                yn = yn(min_idx);
                %plot([x xn], [y yn], 'b');
                %plot(xn, yn, 'm.');
            else
                %plot(x, y, 'rx');
                missing_pts = missing_pts + 1;
            end
            %vessels_scale{ii}(jj,:) = [x y];
        end
        
        display(['Vessel ' num2str(ii) ' missing ' num2str(missing_pts) ' out of ' num2str(num_pts) ' (' num2str(100*missing_pts/num_pts) '%)']);
        
        %plot(vessels{ii}(:,1), vessels{ii}(:,2), 'g');
        %plot(vessels_scale{ii}(:,1), vessels_scale{ii}(:,2), 'r');
    end        
        %figure; imagesc(line_strength); axis image; colormap(gray(256)); hold on;
        %plot(x, y, 'g.');
end
%%
scale = 2;
for noise = 1:8
    test_nailfold = imresize(nailfold, 1 / scale, 'bilinear');
    test_nailfold = imnoise(test_nailfold, 'gaussian', 0, 10^(noise-10));
    [line_strength, line_ori] =...
        gaussian_2nd_derivative_line(test_nailfold, 4 / scale);
    
%     vessel_nms = mb_non_maximal_supp(...
%         imresize(line_strength, [r c], 'bilinear'),...
%         imresize(line_ori, [r c], 'nearest'));
    
    vessel_nms = mb_non_maximal_supp(line_strength, line_ori);
    
    [y_vessels x_vessels] = find(vessel_nms);
    %vessels_scale = vessels;
    
    %snap the [x y] point to nearest potential vessel
    figure; imagesc(test_nailfold); axis image; colormap(gray(256)); hold on;
    plot(x_vessels, y_vessels, 'm.', 'markersize', 2);
    
%     display(['Noise ' num2str(noise)]);
%     for ii = 1:length(vessels)
%         %plot(vessels{ii}(:,1), vessels{ii}(:,2), 'y');
%         %plot(vessels{ii}(:,1), vessels{ii}(:,2), 'gx');
%         missing_pts = 0;
%         num_pts = size(vessels{ii},1);
%         for jj = 1:num_pts
%             x = vessels{ii}(jj,1);
%             y = vessels{ii}(jj,2);
%             
%             neighbours = (abs(x_vessels - x) < 5) & (abs(y_vessels - y) < 5);
% 
%             if any(neighbours)
%                 xn = x_vessels(neighbours);
%                 yn = y_vessels(neighbours);
% 
%                 dists = (xn - x).^2 + (yn - y).^2;
%                 [dummy min_idx] = min(dists);
%                 xn = xn(min_idx);
%                 yn = yn(min_idx);
%                 %plot([x xn], [y yn], 'b');
%                 %plot(xn, yn, 'm.');
%             else
%                 %plot(x, y, 'rx');
%                 missing_pts = missing_pts + 1;
%             end
%             %vessels_scale{ii}(jj,:) = [x y];
%         end
%         
%         display(['Vessel ' num2str(ii) ' missing ' num2str(missing_pts) ' out of ' num2str(num_pts) ' (' num2str(100*missing_pts/num_pts) '%)']);
%         
%         %plot(vessels{ii}(:,1), vessels{ii}(:,2), 'g');
%         %plot(vessels_scale{ii}(:,1), vessels_scale{ii}(:,2), 'r');
%     end        
        %figure; imagesc(line_strength); axis image; colormap(gray(256)); hold on;
        %plot(x, y, 'g.');
end
%%
% Why are some loops missed?
nailfold = imread('C:\isbe\nailfold\anonymous_bmap_files_oct\OCT004rptV2LD4X3LrgMosaic.bmp');
figure; imagesc(nailfold); axis image; colormap(gray(256)); hold on;
    
%Save both line and orientation maps
line_strength = u_load('C:\isbe\nailfold\anonymous_bmap_files\line\g2dOCT004rptV2LD4X3LrgMosaic._line.mat');
line_ori = u_load('C:\isbe\nailfold\anonymous_bmap_files\ori\g2dOCT004rptV2LD4X3LrgMosaic._ori.mat');

col = floor(size(line_strength,2)/2) ;

%Apply NMS to full size maps and rescaled half-size maps
vessel_nms1 = mb_non_maximal_supp(line_strength(:,1:col), line_ori(:,1:col));

%Get the non-zero points from both
[y1 x1] = find(bwareaopen(vessel_nms1 > 0, 20)); clear vessel_nms1;

%Apply NMS to full size maps and rescaled half-size maps
vessel_nms2 = mb_non_maximal_supp(line_strength(:,col+1:end), line_ori(:,col+1:end));

%Get the non-zero points from both
[y2 x2] = find(bwareaopen(vessel_nms2 > 0, 20)); clear vessel_nms2;

plot(x1, y1, 'g.', 'markersize', 2);
plot(x2+col, y2, 'g.', 'markersize', 2);
%%
nailfold_small = nailfold(201:700,1001:1800);
line_strength_small = line_strength(201:700,1001:1800);
line_ori_small = line_ori(201:700,1001:1800);
clear line_strength line_ori

vessel_nms = mb_non_maximal_supp(line_strength_small, line_ori_small);
[y_vessels x_vessels] = find(bwareaopen(vessel_nms > 0, 20));
    
figure; imagesc(nailfold_small); axis image; colormap(gray(256)); hold on;
plot(x_vessels, y_vessels, 'r.', 'markersize', 2);
%%
tic; [line_strength_small, line_ori_small, line_scale_small] =...
        gaussian_2nd_derivative_line(nailfold_small, [4 8]); toc;
tic; [line_strength_small2, line_ori_small2, line_scale_small2] =...
        gaussian_2nd_derivative_line2(nailfold_small, [4 8]); toc;
    
vessel_nms = mb_non_maximal_supp(line_strength_small, line_ori_small);
[y_vessels x_vessels] = find(bwareaopen(vessel_nms > 0, 20));


vessel_nms2 = mb_non_maximal_supp(line_strength_small2, line_ori_small2);
[y_vessels2 x_vessels2] = find(bwareaopen(vessel_nms2 > 0, 0));

figure; imagesc(nailfold_small); axis image; colormap(gray(256)); hold on;
plot(x_vessels, y_vessels, 'r.', 'markersize', 2);
figure; imagesc(nailfold_small); axis image; colormap(gray(256)); hold on;
plot(x_vessels2, y_vessels2, 'g.', 'markersize', 2);

figure; imagesc(log(line_scale_small)); axis image; colormap([0 0 0; 1 0 0; 0 1 0; 0 0 1]); hold on;
figure; imagesc(log(line_scale_small2)); axis image; colormap([0 0 0; 1 0 0; 0 1 0; 0 0 1]); hold on;
%%
figure; imagesc(nailfold_small); axis image; colormap(gray(256)); hold on;
plot(x_vessels, y_vessels, 'r.', 'markersize', 2);
colors([1 2 4 8]) = 'rgbm';
%%
figure; imagesc(nailfold_small); axis image; colormap(gray(256)); hold on;
vessel_mask = false(size(nailfold_small));
vessel_strength = zeros(size(nailfold_small));
[line_strength1] = gaussian_2nd_derivative_line(nailfold_small, 1);
for scale = [2 4 8]
    [line_strength, line_ori] =...
        gaussian_2nd_derivative_line2(nailfold_small, scale);
    
    vessel_nms = mb_non_maximal_supp(line_strength, line_ori);
    %[y_vessels x_vessels] = find(bwareaopen(vessel_nms & line_scale==scale, 20));
    if scale == 22
        scale_mask = (line_scale_small2==scale) & (line_strength > line_strength1);
    else
        scale_mask = line_scale_small2==scale;
    end
    
    [y_vessels x_vessels] = find(vessel_nms);
    figure; imagesc(line_strength); axis image; colormap(gray(256)); hold on; colorbar;
    plot(x_vessels, y_vessels, [colors(scale) '.'], 'markersize', 2);
    
    [y_vessels x_vessels] = find(vessel_nms & scale_mask);
    vessel_mask = vessel_mask | (vessel_nms & scale_mask);
    vessel_strength(scale_mask) = vessel_strength(scale_mask) + vessel_nms(scale_mask);
end
vessel_strength(vessel_strength < 0) = 0;
[y_vessels x_vessels] = find(bwareaopen(vessel_mask, 5, 8));

figure; imagesc(nailfold_small); axis image; colormap(gray(256)); hold on;
plot(x_vessels, y_vessels, 'r.', 'markersize', 2);
figure; imagesc(vessel_strength); axis image; colormap(gray(256)); hold on;
%%
[line_strength_small4, line_ori_small4, line_scale_small4] =...
        gaussian_2nd_derivative_line2(nailfold_small, 4);    
vessel_nms4 = mb_non_maximal_supp(line_strength_small4, line_ori_small4);

[line_strength_small8, line_ori_small8, line_scale_small8] =...
        gaussian_2nd_derivative_line2(nailfold_small, 8);    
vessel_nms8 = mb_non_maximal_supp(line_strength_small8, line_ori_small8);

%mask4 = line_strength_small4 >= line_strength_small8;
mask4 = vessel_nms4 >= vessel_nms8;

vessel_nms = zeros(size(nailfold_small));
vessel_nms(mask4) = vessel_nms4(mask4);
vessel_nms(~mask4) = vessel_nms8(~mask4);

[y_vessels x_vessels] = find(bwareaopen(vessel_nms > 0, 10));



figure; imagesc(nailfold_small); axis image; colormap(gray(256)); hold on;
plot(x_vessels, y_vessels, 'r.', 'markersize', 2);
%%
for ii = 0:15:165
    bar = create_gauss_bar(4, 1, ii, 256, 256, 128, 128);
    [line_strength, line_ori] =...
        gaussian_2nd_derivative_line2(bar, 4);
    theta = mod(180*line_ori(128,128)/pi,180);
    display(['True ori = ' num2str(ii) ', predicted ori = ' num2str(theta)]);
end