clc;
clear;
% close all;
timebar('closeall');

imgroot = 'U:\projects\nailfold\capture\2012_10_22\Left\Digit4\x300\';
imgpath = fullfile(imgroot, 'seq1\preprocessed\registered_g1d\masked');
imgdir = dir(fullfile(imgpath,'*.png'));

patch_sz = 10; t_rng = -patch_sz:patch_sz;
search_sz = 10; s_rng = -search_sz:search_sz;
im_rng = -(patch_sz+search_sz):(patch_sz+search_sz);

match_type = 'normxcorr2';

[xx,yy] = meshgrid(s_rng, s_rng);
ind0 = find(xx==0 & yy==0);

img1 = mean(imread(fullfile(imgpath, imgdir(1).name)), 3);
xmin = 1; xmax = size(img1,2);
ymin = 1; ymax = size(img1,1);

step = 8;
x_rng = xmin:step:xmax;
y_rng = ymin:step:ymax;

x_rng = 300:500;
y_rng = 216;

n_frames = 20;
dx = zeros(size(x_rng));
for i = 2:n_frames
    img2 = mean(imread(fullfile(imgpath, imgdir(i).name)), 3);
%     img2 = img1;
    
    for y = y_rng
        for x = x_rng
            if (any((x+im_rng)<1 | ...
                    (x+im_rng)>size(img2,2)))
                continue;
            end
            search_region = img2(y, x+im_rng);

            template = img1(y, x+t_rng);

            if all(template(:) == template(1))
                continue;
            end

            xc = normxcorr2(double(template), double(search_region));
            
            c = (size(xc,2)+1)/2;
            xc_valid = xc(1, c+s_rng);
            
            [ignore, index] = max(xc_valid);
            dx(find(x_rng == x)) = xx(1, index);

%             figure(1); clf; hold on; colormap(gray(256));
%                 imagesc(img1);
%                 rectangle('edgecolor', 'r', ...
%                           'position', [x-patch_sz-search_sz, y-patch_sz, ...
%                                        2*(patch_sz+search_sz)+1, 2*patch_sz+1]);
%                 rectangle('edgecolor', 'g', ...
%                           'position', [x-patch_sz, y-patch_sz, ...
%                                        2*patch_sz+1, 2*patch_sz+1]);
%                 
%             figure(2);
%                 n = 4;
%                 subplot(n,1,1); hold on; 
%                     plot(img1(y_rng, x_rng), 'r-');
%                     plot(img2(y_rng, x_rng), 'g-');
%                 subplot(n,1,2); plot(-patch_sz:patch_sz, template);
%                     ax = axis(); ax(1:2) = [-(patch_sz+search_sz) (patch_sz+search_sz)]; axis(ax);
%                 subplot(n,1,3); plot(-(patch_sz+search_sz):(patch_sz+search_sz), search_region);
%                 subplot(n,1,4); plot(xx(1,:), xc_valid); 
%                     axis([-patch_sz patch_sz -1 1]);
        end
    end

    figure(2); clf; hold on;
        mn = mean(img1(y_rng, x_rng));
        plot(x_rng, img1(y_rng, x_rng)-mn, 'r-');
        plot(x_rng, img2(y_rng, x_rng)-mn, 'g-');
        plot(x_rng, dx, 'k-');
        plot(x_rng, 0, 'k:');

    img1 = img2;
end
