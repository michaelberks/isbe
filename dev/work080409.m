sw = sample_window(synthesised_image, 11, row, col);
figure; imagesc(sw); axis image;
figure; imagesc(reshape(sampled_window(1:121), 11, 11)); axis image;

[test_p test_chosen_cluster] = mb_gmm_sample('Model', conditioned_model, 'NumSamples', 7);
test_p'

%%
for ii = 1:10
    
    b1 = rand(1);
    b2 = rand(1);
    k = 10;
    a1 = exp(-k*b1);
    a2 = exp(-k*b2);
    
    display(['a1 = ', num2str(a1), ' and ' num2str(a2^(b1/b2)), ' a2 = ', num2str(a2)]);
end
%%
mean_ratios = zeros(15, 4,5);
sd_ratios = zeros(15, 4,5);

for pp = 1:15
    pyr_orig = u_load(['C:\isbe\dev\background\pyramid\normal_2\normal', zerostr(pp,3), '.bmp_pyramid.mat']);
    
    [rows cols] = size(pyr_orig{2,1});
    row_centre = round(rows / 2);
    col_centre = round(cols / 2);

    % Make a the biggest circular mask
    m = min([rows cols]);
    rad = floor((m - 128) / 2);
    [x y] = meshgrid(1:cols, 1:rows);
    filled_image2 = (x - col_centre).^2 + (y - row_centre).^2 < rad.^2;
    clear x y m rad row_centre col_centre;
    filled_image3 = imresize(filled_image2, size(pyr_orig{3,1}));
    filled_image4 = imresize(filled_image2, size(pyr_orig{4,1}));
    filled_image5 = imresize(filled_image2, size(pyr_orig{5,1}));
    filled_image6 = imresize(filled_image2, size(pyr_orig{6,1}));
    
    for ori = 1:5
        
        mean_ratios(pp, 1, ori) = mean(pyr_orig{2,ori}(filled_image2)) / mean(pyr_orig{3,ori}(filled_image3));
        mean_ratios(pp, 2, ori) = mean(pyr_orig{3,ori}(filled_image3)) / mean(pyr_orig{4,ori}(filled_image4));
        mean_ratios(pp, 3, ori) = mean(pyr_orig{4,ori}(filled_image4)) / mean(pyr_orig{5,ori}(filled_image5));
        mean_ratios(pp, 4, ori) = mean(pyr_orig{5,ori}(filled_image5)) / mean(pyr_orig{6,ori}(filled_image6));
        
        sd_ratios(pp, 1, ori) = std(pyr_orig{2,ori}(filled_image2)) / std(pyr_orig{3,ori}(filled_image3));
        sd_ratios(pp, 2, ori) = std(pyr_orig{3,ori}(filled_image3)) / std(pyr_orig{4,ori}(filled_image4));
        sd_ratios(pp, 3, ori) = std(pyr_orig{4,ori}(filled_image4)) / std(pyr_orig{5,ori}(filled_image5));
        sd_ratios(pp, 4, ori) = std(pyr_orig{5,ori}(filled_image5)) / std(pyr_orig{6,ori}(filled_image6));
    end
    
    clear pyr_orig filled_image*
end

%%
f1 = fftshift(fft2(i1));
f1_r = real(f1);
f1_i = imag(f1);
i2 = 10*i1;
%i2 = imresize(i1, 0.5);
f2 = fftshift(fft2(i2));
f2_r = real(f2);
f2_i = imag(f2);

f3 = 10*f1;
i3 = ifft2(ifftshift(f3))

%%
dims = size(i1);
ctr = ceil((dims+0.5)/2);

[xramp,yramp] = meshgrid( ([1:dims(2)]-ctr(2))./(dims(2)/2), ...
    ([1:dims(1)]-ctr(1))./(dims(1)/2) );
angle = atan2(yramp,xramp);
log_rad = sqrt(xramp.^2 + yramp.^2);
log_rad(ctr(1),ctr(2)) =  log_rad(ctr(1),ctr(2)-1);
log_rad  = log2(log_rad);

% Radial transition function (a raised cosine in log-frequency):
[Xrcos,Yrcos] = rcosFn(1,(-1/2),[0 1]);
Yrcos = sqrt(Yrcos);

YIrcos = sqrt(1.0 - Yrcos.^2);
low_mask = pointOp(log_rad, YIrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);
imdft = fftshift(fft2(i1));
low_dft =  imdft .* low_mask;
figure; imagesc(ifft2(ifftshift(low_dft))); axis image; colormap(gray(256)); colorbar;
figure; imagesc(low_mask); axis image; colormap(jet(256)); colorbar;

%%
for lev = 2:2;
    for ori = 1:3;

        temp_band = pyr{lev,ori};
        [rows cols] = size(temp_band);
        row_centre = round(rows / 2);
        col_centre = round(cols / 2);

        % Make a the biggest circular mask
        m = min([rows cols]);
        rad = floor((m - 2^(9-lev) ) / 2);
        [x y] = meshgrid(1:cols, 1:rows);
        filled_image = (x - col_centre).^2 + (y - row_centre).^2 < rad.^2;
        clear x y m rad row_centre col_centre;


        [rows cols] = find(filled_image);
        for ii = 1:length(rows);
            temp_band(rows(ii), cols(ii)) = pyr{lev+1,ori}(ceil(rows(ii)/2), ceil(cols(ii)/2));
        end

        figure; imagesc(temp_band); axis image; caxis([-5 5]);
    end
end

%%
[synthesised_image_2_3 cluster_image_2_3] = mb_gmm_tex_synthesis_42(...
    'PathToTextureGMM', 'C:\isbe\dev\background\models\normal_2\normal_2_model_2_3',...
      'TargetImage', temp_band,... 
      'FilledImage', ~filled_image);
%%
clear
im1 = double(imread(['C:\isbe\dev\background\images\normal_2\normal016.bmp']));
im1 = im1(1:256, 1:256);
im2 = im1(1:2:end, 1:2:end);

[x1 y1] = meshgrid(1:256, 1:256);
lowmask1 = (x1 - 128).^2 + (y1 - 128).^2 < 64.^2;

[x2 y2] = meshgrid(1:128, 1:128);
lowmask2 = (x2 - 64).^2 + (y2 - 64).^2 < 64.^2;

figure; imagesc(im1); axis image; colormap(gray(256));
figure; imagesc(im2); axis image; colormap(gray(256));
figure; imagesc(lowmask1); axis image; colormap(gray(256));
figure; imagesc(lowmask2); axis image; colormap(gray(256));
%%

imdft1 = fftshift(fft2(im1));
imdft2 = fftshift(fft2(im2));
imdft12 = imdft1(65:192, 65:192);
figure; imagesc(abs(imdft12)); axis image;
figure; imagesc(abs(imdft2)); axis image;

sm11 = real(ifft2(ifftshift(imdft1 .* lowmask1)));
sm12 = real(ifft2(ifftshift(imdft12 .* lowmask2)));
sm13 = real(ifft2(ifftshift(imdft12 .* lowmask2 / 4)));
sm22 = real(ifft2(ifftshift(imdft2 .* lowmask2)));

figure; imagesc(sm11); axis image; colormap(gray(256));
figure; imagesc(sm12); axis image; colormap(gray(256));
figure; imagesc(sm13); axis image; colormap(gray(256));
figure; imagesc(sm22); axis image; colormap(gray(256));

%%
mb_build_pyramids(...
    'ImageDir', 'C:\isbe\dev\background\images\normal_2\',...
    'OutputDir', 'C:\isbe\dev\background\pyramid\normal_2\');
mb_build_pyramids(...
    'ImageDir', 'C:\isbe\dev\background\images\mass_2\',...
    'OutputDir', 'C:\isbe\dev\background\pyramid\mass_2\');