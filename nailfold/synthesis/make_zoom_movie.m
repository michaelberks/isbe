clc;

imgroot = 'U:\projects\nailfold\synthesis\_movies\zoom';
outpath = fullfile(imgroot, 'output');
if ~exist(outpath, 'dir')
    mkdir(outpath);
else
    delete(fullfile(outpath, 'frame*.png'));
end
outf = 1;

img1 = double( imread(fullfile(imgroot,'fullnail.png')) );
img2 = double( imread(fullfile(imgroot,'dermo.png')) );

[m,n] = size(img1(:,:,1));

imw = 720;
imh = imw / aratio;
aratio = 1.5;

wx = linspace(0, 1, n/2);
wy = linspace(0, 1, m/2)';
imwts = [wy; wy(end:-1:1)] * [wx wx(end:-1:1)];

% Define the dermoscopy region in the full nail image
x0 = 723;
x1 = x0+180;
y0 = 492; 
y1 = y0 + (x1-x0)/aratio;

nBlend = 100;
nPadding = 250;

s = sigmoid(linspace(-5, 5, 2*nBlend));
s(1:nBlend) = 2*s(1:nBlend);
s(nBlend+1:2*nBlend) = 2*s(nBlend+1:2*nBlend) - 1;
s = sqrt(s);
s = s / max(s(:));

for alpha = s(1:nBlend)
% for alpha = sqrt(linspace(0, 1, 100))
    dx0 = x0 - 1;
    dx1 = x1 - n;
    dy0 = y0 - 1;
    dy1 = y1 - m;
    
    % Define the zoomed region of the full nail image
    x0a = 1 + alpha*dx0;
    x1a = n + alpha*dx1;
    y0a = 1 + alpha*dy0;
    y1a = m + alpha*dy1; 

    [xi,yi] = meshgrid(linspace(x0a, x1a, imw), linspace(y0a, y1a, imh));

    % Sample the zoomed region of the full nail image
    ri = interp2(double(img1(:,:,1)), xi, yi, '*linear');
    gi = interp2(double(img1(:,:,2)), xi, yi, '*linear');
    bi = interp2(double(img1(:,:,3)), xi, yi, '*linear');
    zi = cat(3, ri,gi,bi);

    img_zoomed = zi;

    % Define the dermoscopy region in the sampled (zoomed) image
    n2 = x1a-x0a;
    m2 = y1a-y0a;
    scl = (imw / n2);
    
    x0b = 1   + scl*(1-alpha)*dx0;
    x1b = imw + scl*(1-alpha)*dx1;
    y0b = 1   + scl*(1-alpha)*dy0;
    y1b = imh + scl*(1-alpha)*dy1;
    
    [xi,yi] = meshgrid(ceil(x0b):floor(x1b), ceil(y0b):floor(y1b));
    ri = interp2(linspace(x0b, x1b, n), ...
                 linspace(y0b, y1b, m)', ...
                 double(img2(:,:,1)), xi, yi, '*linear');
    gi = interp2(linspace(x0b, x1b, n), ...
                 linspace(y0b, y1b, m)', ...
                 double(img2(:,:,2)), xi, yi, '*linear');
    bi = interp2(linspace(x0b, x1b, n), ...
                 linspace(y0b, y1b, m)', ...
                 double(img2(:,:,3)), xi, yi, '*linear');
    zi = cat(3, ri,gi,bi);
    
    wi = interp2(linspace(x0b, x1b, n), ...
                 linspace(y0b, y1b, m)', ...
                 double(imwts), xi, yi, '*linear');
    
    rrng = ceil(y0b):floor(y1b);
    crng = ceil(x0b):floor(x1b);

    wi = repmat(wi, [1,1,3]) * alpha;
    img_zoomed(rrng, crng, :) = (1-wi) .* img_zoomed(rrng, crng, :) + ...
                                   wi  .* zi;
    
%     figure(2); clf; hold on;
%         image(uint8(img_zoomed));
%         axis('image','ij');
%     drawnow;
        
    filename = sprintf('frame_%04d.png', outf);
    imwrite(uint8(img_zoomed), fullfile(outpath, filename));
    outf = outf + 1;
end

% Sample a new img1
[xi,yi] = meshgrid(linspace(x0a, x1a, n), linspace(y0a, y1a, m));
ri = interp2(double(img1(:,:,1)), xi, yi, '*linear');
gi = interp2(double(img1(:,:,2)), xi, yi, '*linear');
bi = interp2(double(img1(:,:,3)), xi, yi, '*linear');

wi = repmat(imwts, [1,1,3]);
img_blended = (1-wi) .* cat(3, ri,gi,bi) + ...
                 wi  .* double(img2);

cap_frame = 1;
caproot = fullfile(imgroot, 'ncm_video/looped/');
d = dir(fullfile(caproot, 'frame*.png'));
n_cap_frames = length(d);

for alpha = s(nBlend+1:2*nBlend)
    dx0 = x0 - 1;
    dx1 = x1 - n;
    dy0 = y0 - 1;
    dy1 = y1 - m;
    
    % Define the zoomed region of the full nail image
    x0a = 1 + alpha*dx0;
    x1a = n + alpha*dx1;
    y0a = 1 + alpha*dy0;
    y1a = m + alpha*dy1; 

    [xi,yi] = meshgrid(linspace(x0a, x1a, imw), linspace(y0a, y1a, imh));

    imgbg = (1-alpha) * img_blended + ...
               alpha  * img2;
           
    % Sample the zoomed region of the full nail image
    ri = interp2(double(imgbg(:,:,1)), xi, yi, '*linear');
    gi = interp2(double(imgbg(:,:,2)), xi, yi, '*linear');
    bi = interp2(double(imgbg(:,:,3)), xi, yi, '*linear');
    zi = cat(3, ri,gi,bi);
    
    img_zoomed = zi;
    
    % Define the capillaroscopy region in the sampled (dermoscopy) image
    n2 = x1a-x0a;
    m2 = y1a-y0a;
    scl = (imw / n2);
    
    x0b = 1   + scl*(1-alpha)*dx0;
    x1b = imw + scl*(1-alpha)*dx1;
    y0b = 1   + scl*(1-alpha)*dy0;
    y1b = imh + scl*(1-alpha)*dy1;
    
    filename = sprintf('frame_%04d.png', cap_frame);
    img3 = double( imread(fullfile(caproot,filename)) );
    img3 = img3(:,:,1);
    [mc,nc] = size(img3);
    cap_frame = mod(cap_frame, n_cap_frames) + 1;    

    [xi,yi] = meshgrid(ceil(x0b):floor(x1b), ceil(y0b):floor(y1b));
    ri = interp2(linspace(x0b, x1b, nc), ...
                 linspace(y0b, y1b, mc)', ...
                 double(img3(:,:,1)), xi, yi, '*linear');
    zi = cat(3, ri,ri,ri);
    
    wi = interp2(linspace(x0b, x1b, n), ...
                 linspace(y0b, y1b, m)', ...
                 double(imwts), xi, yi, '*linear');
    
    rrng = ceil(y0b):floor(y1b);
    crng = ceil(x0b):floor(x1b);

    wi = repmat(wi, [1,1,3]) * alpha;
    img_zoomed(rrng, crng, :) = (1-wi) .* img_zoomed(rrng, crng, :) + ...
                                   wi  .* zi;
    
%     figure(2); clf; hold on;
%         image(uint8(img_zoomed));
%         axis('image','ij');
%     drawnow;

    filename = sprintf('frame_%04d.png', outf);
    imwrite(uint8(img_zoomed), fullfile(outpath, filename));
    outf = outf + 1;
end

% Get dermoscopy portion
[xi,yi] = meshgrid(linspace(x0, x1, imw), linspace(y0, y1, imh));

% Sample the zoomed region of the full nail image
ri = interp2(img2(:,:,1), xi, yi, '*linear');
gi = interp2(img2(:,:,2), xi, yi, '*linear');
bi = interp2(img2(:,:,3), xi, yi, '*linear');
derm = cat(3, ri,gi,bi);

wx = linspace(0, 1, imw/2);
wy = linspace(0, 1, imh/2)';
imwts = [wy; wy(end:-1:1)] * [wx wx(end:-1:1)];

% Now fade to capillaroscopy images only
for alpha = linspace(0, 1, ceil(nBlend/2))
    % Skip first frame (same as previous)
    if (alpha == 0)
        continue;
    end
    
    filename = sprintf('frame_%04d.png', cap_frame);
    img3 = double( imread(fullfile(caproot,filename)) );
    img3 = img3(:,:,1);
    [mc,nc] = size(img3);
    cap_frame = mod(cap_frame, n_cap_frames) + 1;
    
    [xi,yi] = meshgrid(linspace(1,nc,imw), linspace(1,mc,imh));
    
    zi = interp2(double(img3(:,:,1)), xi, yi, '*linear');
    zi = cat(3, zi,zi,zi);
    
    wi = alpha + repmat(imwts, [1,1,3]) * (1-alpha);
    img_zoomed = (1-wi) .* derm + ...
                    wi  .* zi;
    
%     figure(2); clf; hold on;
%         image(uint8(img_zoomed));
%         axis('image','ij');
%     drawnow;

    filename = sprintf('frame_%04d.png', outf);
    imwrite(uint8(img_zoomed), fullfile(outpath, filename));
    outf = outf + 1;
end

%% Pad the sequence with capillaroscopy frames
for i = 1:nPadding
    filein = sprintf('frame_%04d.png', cap_frame);
    fileout = sprintf('frame_%04d.png', outf);
    copyfile(fullfile(caproot, sprintf('frame_%04d.png', cap_frame)), ...
             fullfile(outpath, sprintf('frame_%04d.png', outf)));
         
    cap_frame = mod(cap_frame, n_cap_frames) + 1;
    outf = outf + 1;
end

