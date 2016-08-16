%Task: find peaks/troughs in vertical profiles sampled from a noisy OCT
%image of the skin

%Load in the image - make grayscale and convert to doubles
oct_im = double(rgb2gray(imread('C:\isbe\nailfold\misc\OCT\001_1.bmp')));

%Median filter the image to reduce specular noise - can play around with
%different size smoothing windows. May make sense to have non-square
%windows, smoothing more vertically?
oct_im_mf7 = medfilt2(oct_im, [7 7]);
oct_im_mf15 = medfilt2(oct_im, [15 15]);

%Look at some raw and smoothed profiles
sample_cols = 128:128:512;
%%
%1 - sampling profiles and then median filtering in 1D
figure;
subplot(1,3,1); imgray(oct_im); title('Original image');
subplot(1,3,2); imgray(medfilt1(oct_im, 7)); title('Median filtered vertically - w = 7');
subplot(1,3,3); imgray(medfilt1(oct_im, 15)); title('Median filtered vertically - w = 15');

figure;
for i_c = 1:4
    profile_raw = oct_im(:, sample_cols(i_c));
    profile_sm7 = medfilt1(profile_raw, 7);
    profile_sm15 = medfilt1(profile_raw, 15);
    
    subplot(2,2,i_c); hold all;
    plot(profile_raw);
    plot(profile_sm7, 'linewidth', 2);
    plot(profile_sm15, 'linewidth', 2);
    xlabel(['Column ' num2str(sample_cols(i_c))]);
    ylabel('OCT signal');
    legend({'Original profile', 'Smoothed w = 7', 'Smoothed w = 15'});
    
    if i_c == 1
        title('Sampling profiles and then apply median filtering in 1D');
    end

end
%%
%2 - median filtering in 2D then sampling profiles
figure;
subplot(1,3,1); imgray(oct_im); title('Original image');
subplot(1,3,2); imgray(oct_im_mf7); title('Median filtered - w = [7 7]');
subplot(1,3,3); imgray(oct_im_mf15); title('Median filtered - w = [15 15]');

figure;
for i_c = 1:4
    profile_raw = oct_im(:, sample_cols(i_c));
    profile_sm7 = oct_im_mf7(:, sample_cols(i_c));
    profile_sm15 = oct_im_mf15(:, sample_cols(i_c));
    
    subplot(2,2,i_c); hold all;
    plot(profile_raw);
    plot(profile_sm7, 'linewidth', 2);
    plot(profile_sm15, 'linewidth', 2);
    xlabel(['Column ' num2str(sample_cols(i_c))]);
    ylabel('OCT signal');
    legend({'Original profile', 'Smoothed w = 7', 'Smoothed w = 15'});
    
    if i_c == 1
        title('Median filtering in 2D and then sampling profiles');
    end

end
%%
% 3 - Finding peaks and troughs, using Gaussian derivatives - try different
% values of sigma. Increasing sigma increases the smoothing prior to
% feature detection, so you detect coarser features, but potentially reduce
% the spatial accuracy of locating them
sigma = 8;
[g dg ddg] = gaussian_filters_1d(8);

figure;
for i_c = 1:4
    profile_sm7 = oct_im_mf7(:, sample_cols(i_c));
    profile_dg = conv(profile_sm7, dg, 'same');
    
    subplot(2,2,i_c); hold all;
    plot(profile_dg, 'linewidth', 2);
    xlabel(['Column ' num2str(sample_cols(i_c))]);
    ylabel('Derivative');
    legend({'1st deriv. (edges)'});
    
    if i_c == 1
        title('Detetcing features by applying Gaussian deriavtives to median smoothed profiles');
    end

end
figure;
for i_c = 1:4
    profile_sm7 = oct_im_mf7(:, sample_cols(i_c));
    profile_ddg = -sigma*conv(profile_sm7, ddg, 'same');
    
    subplot(2,2,i_c); hold all;
    plot(profile_ddg, 'linewidth', 2);
    xlabel(['Column ' num2str(sample_cols(i_c))]);
    ylabel('Derivative');
    legend({'2nd deriv. (peaks/troughs)'});
    
    if i_c == 1
        title('Detetcing features by applying Gaussian deriavtives to median smoothed profiles');
    end

end
figure;
for i_c = 1:4
    profile_sm7 = oct_im_mf7(:, sample_cols(i_c));
    profile_dg = conv(profile_sm7, dg, 'same');
    profile_ddg = -sigma*conv(profile_sm7, ddg, 'same');
    
    subplot(2,2,i_c); hold all;
    plot(profile_dg, 'linewidth', 2);
    plot(profile_ddg, 'linewidth', 2);
    plot(sqrt(profile_dg.^2 + profile_ddg.^2), 'linewidth', 2);
    xlabel(['Column ' num2str(sample_cols(i_c))]);
    ylabel('Derivative');
    legend({'1st deriv. (edges)', '2nd deriv. (peaks/troughs)', 'Combined magnitude'});
    
    if i_c == 1
        title('Detetcing features by applying Gaussian deriavtives to median smoothed profiles');
    end

end
figure;
for i_c = 1:4
    profile_sm7 = oct_im_mf7(:, sample_cols(i_c));
    profile_dg = conv(profile_sm7, dg, 'same');
    profile_ddg = -sigma*conv(profile_sm7, ddg, 'same');
    
    subplot(2,2,i_c); hold all;
    plot(atan2(profile_dg, profile_ddg), 'linewidth', 2);
    xlabel(['Column ' num2str(sample_cols(i_c))]);
    ylabel('Derivative');
    legend({'Combined phase'});
    
    if i_c == 1
        title('Detetcing features by applying Gaussian deriavtives to median smoothed profiles');
    end

end
%%
oct_im_dg = imfilter(oct_im_mf7, dg', 'replicate');
oct_im_ddg = imfilter(oct_im_mf7, ddg', 'replicate');

[~, min_dg] = min(oct_im_dg);
[~, max_dg] = max(oct_im_dg);
[~, min_ddg] = min(oct_im_ddg);
[~, max_ddg] = max(oct_im_ddg);

figure; 
subplot(1,2,1); imgray(oct_im_dg); title('Gaussian 1st Deriv.');
plot(1:512, min_dg, 'r');
plot(1:512, max_dg, 'g');

subplot(1,2,2); imgray(oct_im); title('Original OCT');
plot(1:512, min_dg, 'r');
plot(1:512, max_dg, 'g');

figure; 
subplot(1,2,1); imgray(oct_im_ddg); title('Gaussian 2nd Deriv.');
plot(1:512, min_ddg, 'r');
plot(1:512, max_ddg, 'g');
subplot(1,2,2); imgray(oct_im); title('Original OCT');
plot(1:512, min_ddg, 'r');
plot(1:512, max_ddg, 'g');
%%
oct_im_mag = oct_im_dg.^2 + oct_im_ddg.^2;
[~, max_mag] = max(oct_im_mag(1:200,:));

figure; 
subplot(1,2,1); imgray(oct_im_mag); title('Gaussian 1st/2nd combined magnitude.');
plot(1:512, max_mag, 'g');
subplot(1,2,2); imgray(oct_im); title('Original OCT');
plot(1:512, max_mag, 'g');
%%
[im_h im_w] = size(oct_im);
smallest_offset = min(max_mag);
largest_offset = max(max_mag);

extended_h = im_h + (largest_offset - smallest_offset);

oct_reg = zeros(extended_h, im_w);

r_idx = 1:im_h;
for i_c = 1:im_w
    offset = largest_offset - max_mag(i_c);
    oct_reg(offset+r_idx,i_c) = oct_im(:,i_c);
end

figure; 
subplot(1,2,1); imgray(oct_im); title('Original OCT');
plot(1:512, max_mag, 'g');
subplot(1,2,2); imgray(oct_reg); title('Registered OCT');









