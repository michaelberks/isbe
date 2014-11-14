%Back from holidays and working at last...

% 19th August 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Testing Gabor filter and gradient response at for different line widths
% and lengths

% 1)
%%%%%%%%%%%%%%%

%Build some test patches
test_im = 120*ones(256);
for ii = 0:6
    centre = (ii+1)*32;
    test_im(33:224, centre-ii:centre+ii) = 200;
end
figure; image(test_im); colormap(gray(256)); axis image;
%%
%Gabor filter test image
GaborFilterArgs.ImageIn = test_im;
GaborFilterArgs.Tau = 4;
GaborFilterArgs.Len = 8;
GaborFilterArgs.Normalise = 0;

[gabor_mag gabor_ori] = gabor_filter(GaborFilterArgs);

figure; imagesc(gabor_mag); colormap(gray(256)); axis image;
figure; imagesc(gabor_mag); colormap(hsv(12)); axis image;
%%
%Compute gradient of test image using 2nd derivative of a Gaussian
sigma = 6;
width = 30; %

ssq = sigma^2;
[x,y] = meshgrid(-width:width,-width:width);
GaussDeriv = -x.*exp(-(x.*x+y.*y)/(2*ssq))/(pi*ssq);

%Calulate image gradient
grad_x = imfilter(test_im,GaussDeriv, 'conv','replicate');
grad_y = imfilter(test_im, GaussDeriv', 'conv','replicate');

gradient_ori = atan(grad_y ./ grad_x);
gradient_mag = grad_x.^2 + grad_y.^2;

%
figure; imagesc(gradient_mag); colormap(gray(256)); axis image;
figure; imagesc(gradient_mag); colormap(hsv(12)); axis image;