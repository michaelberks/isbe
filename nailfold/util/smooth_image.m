clc;
clear;
% close all;
timebar('closeall');

imgroot = 'U:\projects\nailfold\capture\2012_10_22\Left.Digit4.x300\';
imgpath = fullfile(imgroot, 'seq2\corrected\registered_g1d');

mosaic = double(imread(fullfile(imgpath,'mosaic.png')));
mosaic = mosaic(:,:,1);

% hw = 50; fw = (2*hw)+1;
% mosaic = zeros(fw,fw);
% mosaic(hw+1,hw+1) = 255;

[nj, ni] = size(mosaic);

figure(1); clf; colormap(gray(256));
    subplot(3,2,1); imagesc(mosaic); axis('image');

%% Smooth with large Gaussian
sigma = 16;
g = isbe_fspecial('gaussian', (5*sigma+1)*ones(1,2), sigma);
gv = g(:,1);

disp('Separable convolution');
tic;
imgout0 = conv2(gv', gv', mosaic, 'same');
toc;

disp('Nonseparable convolution');
tic;
imgout0 = conv2(mosaic, g, 'same');
toc;

subplot(3,2,2); imagesc(imgout0); axis('image');

%% Smooth with repeated [1,1,1] convolution
n = round(sigma*sigma / (2/3));
imgin = mb_pad(mosaic, [1, 1]);
imgout = imgin;

disp('Repeated [1,1,1]');
tic;
for it = 1:n
    % Horizontal filter
    for j = 2:nj+1
        for i = 2:ni+1
            imgout(j, i) = imgin(j, i-1) + ...
                           imgin(j, i  ) + ...
                           imgin(j, i+1);
        end
    end
    imgin = imgout / 3;
    
    % Vertical filter
    for i = 2:ni+1
        for j = 2:nj+1
            imgout(j, i) = ...
                    imgin(j-1, i) + ...
                    imgin(j,   i) + ...
                    imgin(j+1, i);
        end
    end
    imgin = imgout / 3;
end
toc;
imgout = imgout(2:nj+1, 2:ni+1);
subplot(3,2,3); imagesc(imgout); axis('image');

%% Smooth with repeated convolutions, implemented with integral image
n = 6;
w0 = sqrt((12*sigma*sigma/n) + 1);
wL = floor(w0); 
if (mod(wL,2)==0)
    wL = wL - 1; 
end
m = (12*sigma*sigma - n*wL*wL - 4*n*wL - 3*n) / ...
    (-4*wL - 4);
m = round(m);

nf = wL*wL;
fhw = (wL-1)/2;
offset = fhw+2;

imgin = mb_pad(mosaic, [offset, offset]);
imgout = imgin;

done_m = false;
varTot = 0;

disp('Repeated ones(n,n) via intimg');
tic;
for it = 1:n
    intimgin = cumsum(cumsum(imgin,1),2);
    
    % Increase filter halfwidth after m passes
    if (it > m) && (done_m == false)
        done_m = true;
        fhw = fhw + 1;
        wL = wL + 2;
        nf = wL*wL;
    end
    
    for j = (1:nj)+offset
        for i = (1:ni)+offset
            imgout(j, i) = ...
                    intimgin(j+fhw,   i+fhw) - ...
                    intimgin(j+fhw,   i-fhw-1) - ...
                    intimgin(j-fhw-1, i+fhw) + ...
                    intimgin(j-fhw-1, i-fhw-1);
        end
    end
    imgin = imgout / nf;
    
    varTot = varTot + (wL*wL - 1) / 12;
end
toc;
imgout = imgin((1:nj)+offset, (1:ni)+offset);
subplot(3,2,4); imagesc(imgout); axis('image');
% disp(sqrt(varTot));


%% Smooth with repeated [1,2,1] convolution
n = round(sigma*sigma / (1/2));
imgin = mb_pad(mosaic, [1, 1]);
imgout = imgin;

disp('Repeated [1,2,1]');
tic;
for it = 1:n
    % Horizontal filter
    for j = 2:nj+1
        for i = 2:ni+1
            imgout(j, i) = imgin(j, i-1) + ...
                           imgin(j, i  )*2 + ...
                           imgin(j, i+1);
        end
    end
    imgout = imgout / 4;
    imgin = imgout;
    
    % Vertical filter
    for i = 2:ni+1
        for j = 2:nj+1
            imgout(j, i) = imgin(j-1, i) + ...
                           imgin(j  , i)*2 + ...
                           imgin(j+1, i);
        end
    end
    imgout = imgout / 4;
    imgin = imgout;
end
toc;
imgout = imgout(2:nj+1, 2:ni+1);
subplot(3,2,5); imagesc(imgout); axis('image');
    

%% Smooth with recursive filter
imgin = mosaic;
imgout = zeros(size(imgin));

disp('Recursive filtering');
tic;
for i = 1:ni
    imgout(:,i) = smooth_recursive(imgin(:,i), sigma);
end
imgin = imgout;
for j = 1:nj
    imgout(j,:) = smooth_recursive(imgin(j,:), sigma);
end
toc;
subplot(3,2,6); imagesc(imgout); axis('image');
