clear all; clc;

imgroot = 'U:\projects\nailfold\capture\PatientX\2013_10_02.1';
d = dir(fullfile(imgroot, '*.png'));
d = d(1:300);

img = imread(fullfile(imgroot,d(1).name));
img = double(img(:,:,1));
for i = 1:100
    img_cropped = img(:, i:i+500);
    m(i) = mean(img_cropped(:));
    v(i) = std(img_cropped(:));
end

max(v)/min(v)

figure(1); clf;
subplot(2,1,1); 
    plot(m);
subplot(2,1,2); 
    plot(v);
    
return


gain(1) = 1.0;
offset(1) = 0.0;
for i = 1:length(d)
    img = imread(fullfile(imgroot, d(i).name));
    img = double(img(:,:,1));
    
    mn(i) = mean(img(:));
    vr(i) = var(img(:));
    
    if (i > 1)
%         scale = vcl_sqrt(dest_var/src_var);
%         offset = dest_mean - scale*src_mean;

        gain_i = sqrt(vr(i-1)/vr(i));
        offset_i = mn(i-1) - gain_i*mn(i);

        offset(i) = offset(i-1) + gain(i-1)*offset_i;
        gain(i) = gain(i-1) * gain_i;
    end
end

figure(1); clf;
subplot(4,1,1);
    plot(mn);
subplot(4,1,2);
    plot(vr);
subplot(4,1,3);
    plot(gain);
subplot(4,1,4);
    plot(offset);

return;

n = 100;
t = linspace(0, 16*pi, n);
x = 0.1*sin(4*t) + 1;

scl = 1;
for i = 2:n
    scl(i) = scl(i-1) * x(i)/x(i-1);
end

figure(1); hold on;
plot(x, 'b.-');
plot(scl, 'r.-');
axis([0,n, 0,2]);

return;

imgroot = 'U:\projects\nailfold\capture\PatientX\2013_10_02.1';
d = dir(fullfile(imgroot, 'frame*.png'));

count = 0;
Imag = [];
for i = 1:100:length(d)
    img = imread(fullfile(imgroot, d(i).name));

    for j = 1:8
        img = conv2([1 2 1]/4, [1 2 1]/4, img, 'same');
    end
    img = img(2:2:end-1, 2:2:end-1);

    Ix = conv2([1 2 1]/4, [1 0 -1]/2, img, 'same');
    Iy = conv2([1 0 -1]/2, [1 2 1]/4, img, 'same');
    if (isempty(Imag))
        Imag = zeros(size(Ix));
    end
    
    Imag = Imag + sqrt(Ix.*Ix + Iy.*Iy);
    count = count + 1;
end

Imag = Imag / count;
threshold = mean(Imag(:));
mask = (Imag > threshold);

Imag = sqrt(Ix.*Ix + Iy.*Iy);
Imag(mask) = 0;
peaks = ones(size(Imag)-2);
% peaks = (Imag(2:end-1,2:end-1) > Imag(2:end-1,1:end-2)) & ...
%         (Imag(2:end-1,2:end-1) > Imag(2:end-1,3:end  )) & ...
%         (Imag(2:end-1,2:end-1) > Imag(1:end-2,2:end-1)) & ...
%         (Imag(2:end-1,2:end-1) > Imag(3:end  ,2:end-1));
peaks = peaks & ...
        (Imag(2:end-1,2:end-1) > threshold);

figure(1); clf; colormap(gray(256));
% subplot(2,1,1);exit

%     imagesc(mask);
%     axis('image');
% subplot(2,1,2);
    imagesc(peaks);
    axis('image');

return

outroot = 'S:\projects\mammography\matlab\papers\9999 journal (synthesis)';

outpath = fullfile(outroot, 'fig\generate_edge');
delete(fullfile(outpath,'*.*'));

for n = linspace(0, 0.4, 7)
    if ~exist('pathPoints0','var')
        pathPoints0 = generate_path();
    end

    [widthArray, directionArray, innerArray, outerArray] = ...
        generate_edges(pathPoints0, n);

    % Convert from arbitrary units to pixels
    [pathPoints, widthArray, innerArray, outerArray, imageSize] = ...
        scale_vessel(pathPoints0, widthArray, innerArray, outerArray, 160);

    figure(1); clf; hold on;
        plot(pathPoints(:,1), pathPoints(:,2), '.-', 'color', 0.75*[1,1,1]);
        plot(pathPoints(1,1), pathPoints(1,2), 'go');
        plot(pathPoints(end,1), pathPoints(end,2), 'ro');
        plot(innerArray(:,1), innerArray(:,2), 'k.-');
        plot(outerArray(:,1), outerArray(:,2), 'k.-');
        axis('equal', 'ij', 'off', [0 imageSize(2) 0 imageSize(1)]);

    d = dir(fullfile(outpath, '*.pdf'));
    filename = sprintf('path%02d', length(d)+1);
    exportfig(fullfile(outpath, filename));    
end

return

for n = linspace(0, 0.06, 7)
    pathPoints = generate_path([], n);

    % Convert from arbitrary units to pixels
    [pathPoints, widths, inner, outer, imageSize] = ...
        scale_vessel(pathPoints, [], pathPoints, pathPoints, 160);

    figure(1); clf; hold on;
        plot(pathPoints(:,1), pathPoints(:,2), 'k.-');
    %     plot(pathPoints(1,1), pathPoints(1,2), 'go');
    %     plot(pathPoints(end,1), pathPoints(end,2), 'ro');
        axis('equal', 'ij', 'off', [0 imageSize(2) 0 imageSize(1)]);

    outroot = 'S:\projects\mammography\matlab\papers\9999 journal (synthesis)';
    outpath = fullfile(outroot, 'fig\generate_path');

    d = dir(fullfile(outpath, '*.pdf'));
    filename = sprintf('path%02d', length(d)+1);
    exportfig(fullfile(outpath, filename));
end


return

imgpath = 'U:\projects\nailfold\synthesis\20130724T110328\halfsize';
d = dir(fullfile(imgpath,'frame*.png'));

for i = 1:length(d)
    img = imread(fullfile(imgpath, d(i).name));
    figure(1); clf;
        image(uint8(img(40:60,40:60))); axis('image');
end


return

img = zeros(51,51);
k = 1;
img(26-k:26+k, 26-k:26+k) = 1;

flt = fspecial('gaussian', [51,51], 2);

img = 1 - conv2(img, flt, 'same');
img = img / max(img(:));

bg = 0.01 * randn(51,50);
imgStack(:,:,1) = img(:, 1:50) + bg;
imgStack(:,:,2) = img(:, 2:51) + bg;

[Ix, Iy, It] = image_derivatives(imgStack);

step = [];
lambda = 0.00000;
tau = 0.1;
% output = pt_flow_derivs(img1, img2, [], [], 3.0);
output = pt_flow_horn_schunk(Ix(:,:,1), Iy(:,:,1), It(:,:,1), ...
                             [], lambda, tau);

figure(2); colormap(gray(256));
    imagesc(imgStack(:,:,1));
    axis('image');
figure(3); clf;
    image(flow_image([output.u zeros(51,1) output.v]));
    axis('image');


return

imgroot = 'U:\projects\nailfold\synthesis\';
d = dir(fullfile(imgroot, '*T*'));
imgroot = fullfile(imgroot, d(end).name);
imgroot = fullfile(imgroot, 'halfsize');

nFrames = 64;

imgStack = load_image_stack(imgroot, nFrames);
meanIm = mean(imgStack, 3);

imgStack = reshape(imgStack - repmat(meanIm, [1,1,nFrames]), [], nFrames);
[u,s,v] = svd(imgStack, 0);

svec = cumsum(diag(s).^2);
svec = svec / svec(end);

count = sum(svec < 0.99) + 1;

dx = linspace(-1, 1, 101);
for c = 1:count
    for i = 1:length(dx)
        img = reshape( meanIm(:) + dx(i)*s(c,c)*u(:,c), size(meanIm) );
        img = uint8(255 * normim(img, 'stretch_fixed'));
        
        figure(1); clf;
%             colormap(gray(256));
            image(img); axis('image');
            title(sprintf('%i / %i\n', c, count));
        drawnow;
    end
    pause_if(c < count);
end


return


a = sort(randn(2,10000000));

x = linspace(-4,4,101);
figure(1); clf;
    subplot(3,1,1); hist(a(:),x);
    ax = axis;
    subplot(3,1,2); hist(a(1,:),x);
    axis(ax);
    subplot(3,1,3); hist(abs(a(2,:)-a(1,:)),x);
    axis(ax);

return

x = linspace(-5, 5, 101);
s = 1;
g = exp(-0.5 * x.*x / (s*s));

s = 2;
g2 = exp(-0.5 * x.*x / (s*s));

g3 = 0.5*g + 0.5*g2; 

g4 = 0.99*g + 0.01;

figure(5); clf; 
subplot(2,1,1); hold on;
    plot(x, g);
    plot(x, g2, 'r-');
    plot(x, g3, 'g-');
    plot(x, g4, 'c-');
subplot(2,1,s); hold on;
%     plot(x, -log(g));
    plot(x, -log(g2), 'r-');
    plot(x, -log(g3), 'g-');
    plot(x, -log(g4), 'c-');


return

clc;
clear;
% close all;

imgroot = 'U:\projects\nailfold\remote_folder\images\study';

mkproot = 'U:\projects\nailfold\remote_folder\markup\ppyrkotsch';
d = dir(fullfile(mkproot,'*c_markup.txt'));

for i = 1:length(d)
	markup = read_markup_from(fullfile(mkproot, d(i).name));
    
    nVessels = length(markup.vessels);
    if (nVessels == 0)
        [s,t] = strtok(d(i).name, '#');
        [s,t] = strtok(t(2:end), '_');
        filename = fullfile(imgroot, [s, '.png']);

        if exist(filename, 'file')
            img = imread(filename);
            figure(1); clf; colormap(gray(256));
                imagesc(img);
                axis('image');
            disp(markup.image_grade);
        end
    end 
end


return

ggrpath = fullfile(imgpath, 'ggr\Stage0001');

nframes = 30;

for iframe = 1:nframes
    d = dir(fullfile(imgpath, sprintf('*%04d.png', iframe)));
    imgfile = fullfile(imgpath, d(1).name);
    ptsfile = fullfile(ggrpath, [d(1).name,'.pts']);

    img = imread(imgfile);
    pts = read_pts(ptsfile);
    
    x = pts(1,1):2:pts(end,1);
    y = pts(1,2):2:pts(end,2);
    [xx,yy] = meshgrid(x,y);

    if (iframe == 1)
        intimg = interp2(double(img), xx, yy, '*linear');
    else
        intimg = intimg + interp2(double(img), xx, yy, '*linear');
    end
end

intimg = intimg / nframes;
figure(1); clf; colormap(gray(256));
imagesc(intimg);

