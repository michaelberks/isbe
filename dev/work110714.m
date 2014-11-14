caps_list = dir('C:\isbe\nailfold\ncm\*.bmp');
for ii = 1:length(caps_list)
    
    cap = imread(['C:\isbe\nailfold\ncm\' caps_list(ii).name]);
    if size(cap,3)==3
        cap = cap(:,:,1);
    end
    %[line_strength, line_orientation, line_scale] = gaussian_2nd_derivative_line(cap, [4 8 16]);
    
    figure; imagesc(cap); axis image; colormap(gray(256));
    %subplot(2,1,1); 
    %subplot(2,1,2); imagesc(line_strength); axis image; colormap(gray(256));
end
%%
cap = double(cap);
mask = cap < 255;
cap(~mask) = 0;
%figure; imagesc(mask); axis image;
win_size = 63;

local_n = imfilter(double(mask), ones(win_size));
local_sum = imfilter(cap, ones(win_size));
local_sum2 = imfilter(cap.^2, ones(win_size));

figure; imagesc(local_n); axis image;
figure; imagesc(local_sum); axis image;
figure; imagesc(local_sum2); axis image;
%
local_mean = local_sum ./ local_n;
local_mean2 = local_sum2 ./ local_n;
local_mean(~mask) = 0;
local_mean2(~mask) = 0;

local_std = sqrt(local_mean2 - local_mean.^2);
local_norm = (cap - local_mean) ./ local_std;

figure; imagesc(local_mean); axis image;
figure; imagesc(local_std); axis image;
%%
figure; 
subplot(2,1,1); imagesc(cap - local_mean); axis image; colormap(gray(256));
subplot(2,1,2); imagesc(local_norm); axis image; colormap(gray(256));
%%
figure; 
subplot(2,1,1); imagesc(cap); axis image; colormap(gray(256));
subplot(2,1,2); imagesc(local_norm); axis image; caxis([-2 2]); colormap(gray(256));
%%
local_med = medfilt2(local_norm, [15 15]);
figure; 
subplot(2,1,1); imagesc(cap); axis image; colormap(gray(256));
subplot(2,1,2); imagesc(local_med); axis image; caxis([-2 2]); colormap(gray(256));
%%
[cap_strength, cap_orientation] = gaussian_2nd_derivative_line(cap, [2 4 8]);
[norm_strength, norm_orientation] = gaussian_2nd_derivative_line(local_norm, [2 4 8]);
figure; 
subplot(2,1,1); imagesc(cap_strength); axis image;
subplot(2,1,2); imagesc(norm_strength); axis image;
%
cap_nms = mb_non_maximal_supp(cap_strength, cap_orientation);
norm_nms = mb_non_maximal_supp(norm_strength, norm_orientation);
figure; 
subplot(2,1,1); imagesc(cap_nms); axis image;
subplot(2,1,2); imagesc(norm_nms); axis image;
%%
figure; 
subplot(2,1,1); imagesc(cap_nms > 0); axis image;
subplot(2,1,2); imagesc(norm_nms > 0); axis image;
%%
cap_nms2 = bwareaopen(cap_nms > 0, 20);
norm_nms2 = bwareaopen(norm_nms > 0, 20);
figure; 
subplot(2,1,1); imagesc(cap_nms2); axis image;
subplot(2,1,2); imagesc(norm_nms2); axis image;
%%
[yc xc] = find(cap_nms2);
[yn xn] = find(norm_nms2);
figure; imagesc(cap); axis image; hold on; colormap(gray(256));
plot(xc, yc, 'r.', 'markersize', 2);
figure; imagesc(local_norm); axis image; colormap(gray(256)); hold on;
plot(xn, yn, 'r.', 'markersize', 2);
%%
cap = imread(['C:\isbe\nailfold\ncm for MB and PT\' caps_list(end).name]);
if size(cap,3)==3
    cap = cap(:,:,1);
end
mask = mb_pad(cap, 1, 255) == 255;
mask = bwselect(mask, 1, 1, 4);
mask = mask(2:end-1, 2:end-1);

mask_inner = imdilate(mask, strel('disk', 20)) & ~mask;
%%
[y_int x_int] = find(mask);
[y x] = find(mask_inner);
z = double(cap(mask_inner));
z_int = griddata(x, y, z, x_int, y_int, 'nearest');

cap_filled = cap;
cap_filled(mask) = z_int;
figure; imagesc(cap_filled); axis image;
%%
%%
%5) Fill image sector by sector with average of internal sector

%Compute orientation of vector from each pixel to centre of image
[row col] = size(cap);
[x, y] = meshgrid(1:col, 1:row);
x = x - col/2;
y = y - row/2;
tan_yx = atan2(y, x);

%Select how many sectors to use - try varying this number, currently set at
%360
n_sectors = 360;
ori = (linspace(-pi,pi, n_sectors));

cap_filled2 = cap;

%Work through each sector, filling the outside pixels with the mean of
%background pixels from the corresponding sector in the inner ring
for ii = 1:n_sectors-1

    sector = (tan_yx >  ori(ii)) & (tan_yx <=  ori(ii+1));
    mean_val = mean(cap(sector & mask_inner));
    
    cap_filled2(sector & ~mask_inner) = mean_val;

end
cap_filled2 = imfilter(cap_filled2, fspecial('gaussian', 11, 2), 'replicate');
cap_filled2(~mask) = cap(~mask);
figure; imagesc(cap_filled2); axis image;


