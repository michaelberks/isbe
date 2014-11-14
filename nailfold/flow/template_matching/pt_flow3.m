clc;
clear;
% close all;
timebar('closeall');

imgroot = 'U:\projects\nailfold\capture\2012_10_22\Left\Digit4\x300\';
imgpath = fullfile(imgroot, 'seq1\preprocessed\registered_g1d\masked');
imgdir = dir(fullfile(imgpath,'*.png'));

patch_sz = 10;
search_sz = 10;
match_type = 'normxcorr2';

rng = -search_sz:search_sz;
[xx,yy] = meshgrid(rng, rng);
ind0 = find(xx==0 & yy==0);

img1 = mean(imread(fullfile(imgpath, imgdir(1).name)), 3);
xmin = 1; xmax = size(img1,2);
ymin = 1; ymax = size(img1,1);

step = 8;
x_rng = xmin:step:xmax;
y_rng = ymin:step:ymax;

x_rng = 390:500;
y_rng = 216;

n_frames = 20;
dxy = cell(ymax, xmax);
for i = 2:n_frames
    img2 = mean(imread(fullfile(imgpath, imgdir(i).name)), 3);
    
    for y = y_rng
        for x = x_rng
            rng = -(patch_sz+search_sz):(patch_sz+search_sz);
%             if (any((x+rng)<1 | (x+rng)>size(img2,2) | ...
%                     (y+rng)<1 | (y+rng)>size(img2,1)))
            if (any((x+rng)<1 | (x+rng)>size(img2,2)))
                continue;
            end
%             search_region = img2(y+rng, x+rng);
            search_region = img2(y, x+rng);

            rng = -patch_sz:patch_sz;
%             template = img1(y+rng, x+rng);
            template = img1(y, x+rng);

            if all(template(:) == template(1))
                continue;
            end

            xc = normxcorr2(double(template), double(search_region));
            
            c = (size(xc,2)+1)/2;
            rng = -search_sz:search_sz;
%             xc_valid = xc(c(1)+rng, c(2)+rng);
            xc_valid = xc(1, c+rng);

            figure(2);
                subplot(3,1,1); plot(-patch_sz:patch_sz, template);
                subplot(3,1,2); plot(-(patch_sz+search_sz):(patch_sz+search_sz), search_region);
                subplot(3,1,3); plot(xx(1,:), xc_valid);
            continue;

            % Accumulate correlation scores over time
            if isempty(dxy{y,x})
                dxy{y,x} = xc_valid;
            else
                dxy{y,x} = dxy{y,x} + xc_valid;
            end                
        end
    end

    img1 = img2;
end

dxy = dxy(y_rng, x_rng);

dx_img = nan(length(y_rng), length(x_rng));
dx_wt = zeros(size(dx_img));

for y = 1:length(y_rng)
    for x = 1:length(x_rng)
        if ~isempty(dxy{y, x})
            [xc_max, index] = max(dxy{y, x}(:));
            dx_img(y,x) = xx(index);
            dx_wt(y,x) = xc_max / (n_frames-1);
        end
    end
end

alpha = 0.5;
underlay = 255 * normim(img1(y_rng, x_rng, :));
underlay = repmat(alpha * double(underlay), [1,1,3]);

% dx_img = dx_img(y_rng, x_rng);
dx_img(1) = -1; dx_img(end) = 1;
dx2 = 255 * normim(dx_img, 'stretch_fixed');
dx_rgb = 255 * ind2rgb(ceil(dx2), redgreen());

figure(1); colormap(gray(256));
%     image(uint8(underlay)); axis('image');
    image(uint8(img1)); axis('image');
figure(3); 
    dx_overlay = (1-alpha) * dx_rgb .* repmat(dx_wt, [1,1,3]);
    image(uint8(underlay + dx_overlay)); axis('image');
figure(4); colormap(gray(256));
    imagesc(dx_wt); axis('image');
  
outpath = 'U:\projects\nailfold\optic_flow';
imwrite(uint8(255*normim(dx_wt)), ...
        fullfile(outpath, 'weights.png'));
imwrite(uint8(underlay + dx_overlay), ...
        fullfile(outpath, 'overlaid.png'));