clc;
clear;
% close all;

datroot = 'U:\projects\nailfold\capture';
imgpath = fullfile(datroot, '\2012_10_18\Left\Digit4\x300\seq1');

nframes = 30;

for iframe = 2 % 1:nframes
    d = dir(fullfile(imgpath, sprintf('*%04d.png', iframe)));
    imgfile = fullfile(imgpath, d(1).name);
    img = mean(double(imread(imgfile)),3);
    img = img(:,1:480);

    imsz = 2*round(size(img)/2);

    if 0
        figure(2); clf; colormap(gray(256));
            imagesc(img);
            [x,y] = ginput(2);
            close();
    else
        x = [393,425];
        y = [290,289];
    end
        
    wiener = zeros(size(img));
%     dx = x(2)-x(1);
%     dy = y(2)-y(1);
%     for xi = round(x(1)):round(x(2))
%         yi = round(y(1) + (xi-x(1)) * dy/dx);
%         wiener(yi,xi) = 1;
%     end
%     wiener = wiener / dx; % kernel should sum to unity?

%     cx = size(wiener,2)/2;
%     cy = size(wiener,1)/2;
%     wiener(cy, (cx-16):(cx+17)) = 1;
    wiener = zeros(32,32);
    wiener(16,:) = 1;
    wiener = wiener / sum(wiener(:));

    [img2,wiener] = deconvblind(img,wiener,10, 0.1);
    newimg = conv2(img2, wiener, 'same');
    
    figure(1); clf; colormap(gray(256));
        subplot(2,2,1); imagesc(img); axis('image');
        subplot(2,2,2); imagesc(img2); axis('image');
        subplot(2,2,3); imagesc(newimg); axis('image');
    
    return
    
%     img = img - mean(img(:)); % zero-mean image
    
    f1 = fft2(img);
    
    % Deconvolve image
    h = fft2(wiener);
    k = 1e0;
    f1 = (f1.*conj(h)) ./ (h*h' + k);
    
    mag1 = abs(f1);
    phase1 = angle(f1);
    recon = real(ifft2(mag1 .* complex(cos(phase1), sin(phase1))));
%     recon = fftshift(recon);

    f1 = fftshift(f1);
    mag1 = abs(f1);
    phase1 = angle(f1);
    
%     % Find central frequency (=0)
%     [cy,cx] = find(mag1 < 1e-3);
%     w = 5;
%     window = mag1(cy-w:cy+w, cx-w:cx+w);
%     power = sum(window(:))
    
    figure(1); clf; colormap(gray(256)); 
        subplot(2,2,1); imagesc(img); axis('image');
        subplot(2,2,2); imagesc(recon); axis('image');
        subplot(2,2,3); imagesc(log(mag1)); axis('image');
        subplot(2,2,4); imagesc(phase1); axis('image');

%     figure(1); clf; colormap(gray(256));
%         imagesc(img);

%     pause;
end

return

datroot = 'U:\projects\nailfold\capture';
imgpath = fullfile(datroot, '\2012_10_18\Left\Digit4\x300\seq2');
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

