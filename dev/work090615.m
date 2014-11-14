im1 = double(rot90(imresize(imread('C:\isbe\density\mammograms\018LCC1824.tif'),0.25, 'bilinear'),-1));
figure; imagesc(im1); axis image; colormap(gray(256));
mask = roipoly;
dt1 = dtwavexfm2(im1, 6);
[ilp1 icp1] = mb_dual_tree_transform(dt1);
%%
for lev = 2:6;
    
    lev_mask = imresize(mask, size(dt1{lev}(:,:,1)));
    
    max_dt = max(dt1{lev}, [], 3);
    max_dt(~lev_mask) = 0;
    figure; image(complex2rgb(max_dt, [-pi, pi])); axis image;
    title(['Dual-tree, level ', num2str(lev)]);
    
    max_ilp = max(ilp1{lev}, [], 3);
    figure; imagesc(abs(angle(max_ilp)) < pi/6); axis image;
    
    max_ilp(~lev_mask) = -1;
    %figure; image(complex2rgb(max_ilp)); axis image;
    figure; imagesc(angle(max_ilp)); colormap hsv; caxis([-pi pi]); axis image;
    title(['ILP, level ', num2str(lev)]);
    
    max_icp = max(icp1{lev}, [], 3);
    max_icp(~lev_mask) = 0;
    %figure; image(complex2rgb(max_icp)); axis image;
    figure; imagesc(angle(max_icp)); colormap hsv; axis image;
    title(['ICP, level ', num2str(lev)]);
    
end
%%
mask5 = imresize(mask, size(dt1{5}(:,:,1)));
max_dt5 = max(dt1{5}, [], 3);
max_ilp5 = max(ilp1{5}, [], 3);
max_icp5 = max(icp1{5}, [], 3);
%%
figure; imagesc(abs(angle(max_ilp5)) < pi/6); axis image;
%%
figure; imagesc(im1); axis image; colormap(gray(256)); hold on;
[mask, xi, yi] = roipoly;
%
xi(end) = []; yi(end) = [];
d = [0; cumsum(sum(diff([xi yi]).^2,2))];
[xy] = interp1(d, [xi yi], linspace(d(1), d(end), 200), 'spline');
hold on; plot(xy(:,1), xy(:,2), 'g');
%%
[im_r im_c] = size(im1);

%Now take normal profiles along the xy border
[fx, gxy] = gradient(xy); clear fx;

%normalise gxy
gxy = gxy ./ [sqrt(sum(gxy.^2, 2)), sqrt(sum(gxy.^2, 2))];
%%
for ii = 1:length(gxy(:,1)) %= number of rows in xy
    n1_x = xy(ii, 1) - 80*gxy(ii, 2);
    n1_y = xy(ii, 2) + 80*gxy(ii, 1);
    n2_x = xy(ii, 1) + 20*gxy(ii, 2);
    n2_y = xy(ii, 2) - 20*gxy(ii, 1);

    if (n1_x >= 1 && n1_x < im_c && n2_x >= 1 && n2_x < im_c && n1_y >= 1 && n1_y < im_r && n2_y >= 1 && n2_y < im_r)
        [cx, cy, cp] = improfile(im1, [n1_x, n2_x], [n1_y, n2_y], 100);
        normal_p(ii, :) = cp';
        normal_x(ii, :) = cx';
        normal_y(ii, :) = cy';

        %[hist_n, hist_x] = hist(cp, 16);
        %hist_mode = hist_x(find(hist_n == max(hist_n(1:2)), 1) + 1);
        %border_idx = find(cp <= hist_mode, 1);
        %xy(ii,1) = cx(border_idx);
        %xy(ii,2) = cy(border_idx);
    end

    if 1%not(rem(ii, 10))
        plot(cx, cy, 'b');
        plot(n1_x, n1_y, 'gx');
        plot(n2_x, n2_y, 'rx');
    end
end
%%
%%
mam_dir = 'C:\isbe\mammograms\new_CAD\bmP_2004\';
mam_list = dir('C:\isbe\mammograms\new_CAD\bmP_2004\*.bmp');
%

for ii = 1:2%length(mam_list);
    mlo = ~isempty(strfind(mam_list(ii).name, 'ML'));
    
    mam = double(imresize(imread([mam_dir, mam_list(ii).name]), [1024 NaN], 'bilinear'));
    if ~isempty(strfind(mam_list(ii).name, 'R'))
        mam = fliplr(mam);
    end
    [segmentation] = segment_breast('image', mam, 'mlo', mlo, 'right', right, 'plot', 0);
    segmentation.size = size(mam);
    %save(['C:\isbe\dev\segmentation\breast_borders\', mam_list(ii).name(1:end-4), '_segmentation.mat'], 'mam');
    
    set(gcf, 'Name', mam_list(ii).name);
    %saveas(gcf, ['C:\isbe\dev\segmentation\figures\', mam_list(ii).name(1:end-3), 'fig']);
    %close(gcf);
end
%%
%%
mam_dir = 'C:\isbe\mammograms\new_CAD\bmP_2004\';
mam_list = dir('C:\isbe\mammograms\new_CAD\bmP_2004\*ML*.bmp');
for ii = 1:10%length(mam_list); %[1:302 304:length(mam_list)];
    try
        mlo = ~isempty(strfind(mam_list(ii).name, 'ML'));
        right = ~isempty(strfind(mam_list(ii).name, 'R'));
        
        mam = double(imresize(imread([mam_dir, mam_list(ii).name]), [1024 NaN], 'bilinear'));

        [segmentation] = segment_breast('image', mam, 'mlo', mlo, 'right', right, 'plot', 1);
        %save(['C:\isbe\dev\segmentation\breast_borders\', mam_list(ii).name(1:end-4), '_segmentation.mat'], 'segmentation');
    catch
        %display(['Skipping ', mam_list(ii).name]);
    end
%     set(gcf, 'Name', mam_list(ii).name);
    %saveas(gcf, ['C:\isbe\dev\segmentation\figures\', mam_list(ii).name(1:end-3), 'fig']);
    %close(gcf);
end
%%
for ii = 11:20%[1:302 304:length(mam_list)];

    mam = double(imresize(imread([mam_dir, mam_list(ii).name]), [1024 NaN], 'bilinear'));
    figure; subplot(1,2,1); imagesc(mam); axis image; colormap(gray(256));
    subplot(1,2,2); imagesc(-mam); axis image; colormap(gray(256));
    set(gcf, 'Name', mam_list(ii).name);
    
end
    
%%
mam_dir = 'F:\';
mam_list = dir('F:\*.tif');
for ii = 1:length(mam_list);
    try
        mlo = ~isempty(strfind(mam_list(ii).name, 'ML'));

        mam = double(imresize(imread([mam_dir, mam_list(ii).name]) / 256, [1024 NaN], 'bilinear'));
        if isempty(strfind(mam_list(ii).name, '2430'))
            if ~isempty(strfind(mam_list(ii).name, 'R'))
                mam = rot90(mam, 1);
                mam = fliplr(mam);
            else
                mam = rot90(mam, -1);
            end
        else
            if ~isempty(strfind(mam_list(ii).name, 'R'))
                mam = fliplr(mam);
            else
                mam = rot90(mam, 2);
            end
        end
        segmentation.breast_border = segment_breast('image', mam, 'mlo', mlo, 'plot', 0);
        segmentation.size = size(mam);
        save(['C:\isbe\dev\segmentation\', mam_list(ii).name(1:end-4), '_segmentation.mat'], 'segmentation');
    catch
        display(['Skipping ', mam_list(ii).name]);

mam_dir = 'C:\isbe\density\mammograms\';
mam_list = dir('C:\isbe\density\mammograms\*2430*.tif');
for ii = 1:length(mam_list) 
    mlo = ~isempty(strfind(mam_list(ii).name, 'ML'));
    right = ~isempty(strfind(mam_list(ii).name, 'R'));
    size2430 = ~isempty(strfind(mam_list(ii).name, '2430'));
    mam = double(imresize(imread([mam_dir, mam_list(ii).name]) / 256, [1024 NaN], 'bilinear'));
    if size2430
        if ~right
            mam = rot90(mam, 2);
        end
    else
        if right
            mam = rot90(mam, 1);
        else
            mam = rot90(mam, -1);
        end

    end

end
%%
mam_dir = 'F:\';
mam_list = dir('F:\*.tif');
for ii = 21:40;
    try
        mlo = ~isempty(strfind(mam_list(ii).name, 'ML'));

        mam = double(imresize(imread([mam_dir, mam_list(ii).name]) / 256, [1024 NaN], 'bilinear'));
        if isempty(strfind(mam_list(ii).name, '2430'))
            if ~isempty(strfind(mam_list(ii).name, 'R'))
                mam = rot90(mam, 1);
                mam = fliplr(mam);
            else
                mam = rot90(mam, -1);
            end
        else
            if ~isempty(strfind(mam_list(ii).name, 'R'))
                mam = fliplr(mam);
            else
                mam = rot90(mam, 2);
            end
        end
        
        load(['C:\isbe\dev\segmentation\', mam_list(ii).name(1:end-4), '_segmentation.mat'], 'segmentation');
        figure; imagesc(mam); axis image; hold on;
        plot(segmentation.breast_border(:,1), segmentation.breast_border(:,2), 'r-');
        plot(segmentation.breast_border(:,1), segmentation.breast_border(:,2), 'rx');
    catch
        display(['Skipping ', mam_list(ii).name]);
    end
end
    [breast_border breast_air] = segment_breast('image', mam, 'mlo', mlo, 'right', right, 'plot', 1);
    segmentation.breast_border = breast_border;
    segmentation.breast_air = breast_air;
    segmentation.size = size(mam);
        save(['C:\isbe\density\segmentation\', mam_list(ii).name(1:end-4), '_segmentation.mat'], 'segmentation');
end

