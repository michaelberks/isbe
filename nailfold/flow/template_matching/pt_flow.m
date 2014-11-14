clc;
clear;
% close all;
timebar('closeall');

imgroot = 'U:\projects\nailfold\capture\2012_10_22\Left\Digit4\x300\';
imgpath = fullfile(imgroot, 'seq1\preprocessed\registered_g1d\masked');

if (false)
    lineim = 128*ones(480,640);
    lineim(:,300:320) = 64;
    img1 = lineim(:,1:end-3);
    img2 = lineim(:,3:end);
else
    img1 = mean(imread(fullfile(imgpath, 'frame_0001.png')), 3);
    img2 = mean(imread(fullfile(imgpath, 'frame_0003.png')), 3);
%     img2 = img1(1:end-2, 1:end-2);
%     img1 = img1(3:end, 3:end);
end

sigma = 0;
img1 = img1 + round(randn(size(img1))*sigma);
img2 = img2 + round(randn(size(img2))*sigma);

% figure(1); colormap(gray(256)); clf; hold on;
%     imagesc(img1); axis('image', 'ij');
% figure(2); colormap(gray(256)); clf; hold on;
%     imagesc(img2); axis('image', 'ij');
% return

patch_sz = 15;
search_sz = 15;
match_type = 'sse';

rng = -search_sz:search_sz;
[xx,yy] = meshgrid(rng, rng);
ind0 = find(xx==0 & yy==0);

xmin = 1; xmax = size(img1,2);
ymin = 1; ymax = size(img1,1);

step = 8;
x_rng = xmin:step:xmax;
y_rng = ymin:step:ymax;

dx_img = zeros(size(img1));
dy_img = zeros(size(img1));

tic;
% n_pts = length(x_rng) * length(y_rng);
x_vec = []; u_vec = [];
y_vec = []; v_vec = [];
for y = y_rng
    for x = x_rng
        rng = -(patch_sz+search_sz):(patch_sz+search_sz);
        if (any((x+rng)<1 | (x+rng)>size(img2,2) | ...
                (y+rng)<1 | (y+rng)>size(img2,1)))
            continue;
        end
        search_region = img2(y+rng, x+rng);
        
        rng = -patch_sz:patch_sz;
        template = img1(y+rng, x+rng);
        
        if all(template(:) == template(1))
            continue;
        end
        
        xc = match_cost(double(template), double(search_region), match_type);
        switch match_type
            case {'normxcorr2'},    xc_min_thresh = -0.5;
            case {'filter2'},       xc_min_thresh = -1e5;
            case {'phasecorr'},     xc_min_thresh = -0.05; % ???
        end
        xc_min_thresh = inf;
                
        c = (size(xc)+1)/2;
        rng = -search_sz:search_sz;
        xc_valid = xc(c(1)+rng, c(2)+rng);
%         xc_valid([1,end],:) = -inf;
%         xc_valid(:,[1,end]) = -inf;

        xc0 = xc_valid(ind0);
        [xc_min, ind] = min(xc_valid(:));
        
        if (xc_min > xc_min_thresh)
            % No strong correlation at all
            dxy = [0, 0];
            col = 'g';
%         elseif (xc_min > xc0*1.1)
%             % No significant difference between displaced and original 
%             % positions
%             dxy = [0, 0];
%             col = 'r';
        else
            % Well-defined movement between patches
            dxy = [xx(ind), yy(ind)];
            col = 'b';
        end
        
        dxy = [xx(ind), yy(ind)];
        dx_img(y,x) = dxy(1);
        dy_img(y,x) = dxy(2);
        
        x_vec(end+1) = x;
        y_vec(end+1) = y;
        u_vec(end+1) = dxy(1);
        v_vec(end+1) = dxy(2);
    end
end
toc;

figure(10); clf; hold on; colormap(gray(256));
    imagesc(img1);
    quiver(x_vec, y_vec, u_vec, v_vec);
    axis('image','ij');

alpha = 0.8;
underlay = 255 * normim(img1(y_rng, x_rng, :));
underlay = repmat(alpha * double(underlay), [1,1,3]);
figure(1); colormap(gray(256)); maximize();
    image(uint8(underlay)); axis('image');

dx = dx_img(y_rng, x_rng);
dx(1) = -1; dx(end) = 1;
dx2 = 255 * normim(dx, 'stretch_fixed');
dx_rgb = 255 * ind2rgb(ceil(dx2), redgreen());
figure(3); maximize(); 
    dx_overlay = (1-alpha) * dx_rgb;
    image(uint8(underlay + dx_overlay)); axis('image');

dy = dy_img(y_rng, x_rng);
dy(1) = -1; dy(end) = 1;
dy2 = 255 * normim(dy, 'stretch_fixed');
dy_rgb = 255 * ind2rgb(ceil(dy2), redgreen());
figure(4); maximize();
    dy_overlay = (1-alpha) * dy_rgb;
    image(uint8(underlay + dy_overlay)); axis('image');
