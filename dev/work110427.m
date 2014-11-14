mam_names = get_mammo_info(dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\*.mat'));
k_dir = 'C:\isbe\asymmetry_project\data\k_stellate_maps\rf_thin\2004_screening_processed\abnormals\';
missing_idx = zeros(292,3);
for ii = 1:292
    missing_idx(ii,1) = ~exist([k_dir mam_names{ii} '_f1.mat'], 'file');
    missing_idx(ii,2) = ~exist([k_dir mam_names{ii} '_f2.mat'], 'file');
    missing_idx(ii,3) = ~exist([k_dir mam_names{ii} '_mask.mat'], 'file');
end
%%
o_map_g2 = load_uint8('Z:\data\orientation_maps\g2d\2004_screening_processed\mass_roi\002LCC_roi.mat');
o_map_rf = load_uint8('Z:\data\orientation_maps\rf\2004_screening_processed\mass_roi\002LCC_roi.mat');
o_map_rt = load_uint8('Z:\data\orientation_maps\rf_thin\2004_screening_processed\mass_roi\002LCC_roi.mat');

l_map_g2 = u_load('Z:\data\line_maps\g2d\2004_screening_processed\mass_roi\002LCC_roi.mat');

a_map_g2 = exp(i * 2 * o_map_g2);
a_map_rf = exp(i * 2 * angle(o_map_rf));
a_map_rt = exp(i * 2 * angle(o_map_rt));

figure; image(complex2rgb(a_map_g2)); axis image;
figure; image(complex2rgb(a_map_rf)); axis image;
figure; image(complex2rgb(a_map_rt)); axis image;
%
a_diff = angle(a_map_g2 .* conj(a_map_rt)) / 2;
figure; imagesc(abs(a_diff)); axis image; colorbar;
%
l_diff = a_diff;
l_diff(~l_map_g2) = 0;
figure; imagesc(abs(l_diff)); axis image; colorbar;
%
xy = repmat(-400:399, 800, 1);
%mask = xy.^2 + xy'.^2 < 200^2;
mask = true(800);

px_per_mm = 100/9;
sigma_range = [3.2 4.0 5.4 6.2 7.8];
r_max = 3*round(px_per_mm*sigma_range);

%Gaussian
[f1_g2 f2_g2 mask] = karssemeijer_radial_projection_multiscale(...
    l_map_g2, o_map_g2, 'mask', mask, 'spacing', 4, 'r_max', r_max(3));

%RF with mag dispersion lines
[f1_rf f2_rf] = karssemeijer_radial_projection_multiscale(...
    abs(o_map_rf), angle(o_map_rf), 'mask', mask, 'spacing', 4, 'r_max', r_max(3));
[f1_rt f2_rt] = karssemeijer_radial_projection_multiscale(...
    abs(o_map_rt), angle(o_map_rt), 'mask', mask, 'spacing', 4, 'r_max', r_max(3));
%%
%RF with mag dispersion lines, angles mod pi
[f1_rf_mp f2_rf_mp] = karssemeijer_radial_projection_multiscale(...
    abs(o_map_rf), mod(angle(o_map_rf),pi), 'mask', mask, 'spacing', 4, 'r_max', r_max(3));
[f1_rt_mp f2_rt_mp] = karssemeijer_radial_projection_multiscale(...
    abs(o_map_rt), mod(angle(o_map_rt),pi), 'mask', mask, 'spacing', 4, 'r_max', r_max(3));

%RF with g2d lines
[f1_rfl f2_rfl] = karssemeijer_radial_projection_multiscale(...
    l_map_g2, angle(o_map_rf), 'mask', mask, 'spacing', 4, 'r_max', r_max(3));
[f1_rtl f2_rtl] = karssemeijer_radial_projection_multiscale(...
    l_map_g2, angle(o_map_rt), 'mask', mask, 'spacing', 4, 'r_max', r_max(3));

%RF with g2d lines, angles mod pi
[f1_rfl_mp f2_rfl_mp] = karssemeijer_radial_projection_multiscale(...
    l_map_g2, mod(angle(o_map_rt),pi), 'mask', mask, 'spacing', 4, 'r_max', r_max(3));
[f1_rtl_mp f2_rtl_mp] = karssemeijer_radial_projection_multiscale(...
    l_map_g2, mod(angle(o_map_rt),pi), 'mask', mask, 'spacing', 4, 'r_max', r_max(3));

%%
clims1 = [...
    min([f1_g2(:); f1_rf(:); f1_rt(:); f1_rfl(:); f1_rtl(:); f1_rtl_mp(:); f1_rt_mp(:); f1_rf_mp(:); f1_rfl_mp(:)])...
    max([f1_g2(:); f1_rf(:); f1_rt(:); f1_rfl(:); f1_rtl(:); f1_rtl_mp(:); f1_rt_mp(:); f1_rf_mp(:); f1_rfl_mp(:)])];
clims2 = [...
    min([f2_g2(:); f2_rf(:); f2_rt(:); f2_rfl(:); f2_rtl(:); f2_rtl_mp(:); f2_rt_mp(:); f2_rf_mp(:); f2_rfl_mp(:)])...
    max([f2_g2(:); f2_rf(:); f2_rt(:); f2_rfl(:); f2_rtl(:); f2_rtl_mp(:); f2_rt_mp(:); f2_rf_mp(:); f2_rfl_mp(:)])];
%
%Gaussian
f1_map_g2 = zeros(size(mask)); f1_map_g2(mask) = f1_g2;
f2_map_g2 = zeros(size(mask)); f2_map_g2(mask) = f2_g2;

%RF with mag dispersion lines
f1_map_rf = zeros(size(mask)); f1_map_rf(mask) = f1_rf;
f2_map_rf = zeros(size(mask)); f2_map_rf(mask) = f2_rf;
f1_map_rt = zeros(size(mask)); f1_map_rt(mask) = f1_rt;
f2_map_rt = zeros(size(mask)); f2_map_rt(mask) = f2_rt;

%RF with mag dispersion lines, angles mod pi
f1_map_rf_mp = zeros(size(mask)); f1_map_rf_mp(mask) = f1_rf_mp;
f2_map_rf_mp = zeros(size(mask)); f2_map_rf_mp(mask) = f2_rf_mp;
f1_map_rt_mp = zeros(size(mask)); f1_map_rt_mp(mask) = f1_rt_mp;
f2_map_rt_mp = zeros(size(mask)); f2_map_rt_mp(mask) = f2_rt_mp;

%RF with g2d lines
f1_map_rfl = zeros(size(mask)); f1_map_rfl(mask) = f1_rfl;
f2_map_rfl = zeros(size(mask)); f2_map_rfl(mask) = f2_rfl;
f1_map_rtl = zeros(size(mask)); f1_map_rtl(mask) = f1_rtl;
f2_map_rtl = zeros(size(mask)); f2_map_rtl(mask) = f2_rtl;

%RF with g2d lines, angles mod pi
f1_map_rfl_mp = zeros(size(mask)); f1_map_rfl_mp(mask) = f1_rfl_mp;
f2_map_rfl_mp = zeros(size(mask)); f2_map_rfl_mp(mask) = f2_rfl_mp;
f1_map_rtl_mp = zeros(size(mask)); f1_map_rtl_mp(mask) = f1_rtl_mp;
f2_map_rtl_mp = zeros(size(mask)); f2_map_rtl_mp(mask) = f2_rtl_mp;

%
%Gaussian
figure; 
subplot(1,2,1); imagesc(f1_map_g2(1:4:799, 1:4:799)); axis image; caxis(clims1); colorbar;
subplot(1,2,2); imagesc(f2_map_g2(1:4:799, 1:4:799)); axis image; caxis(clims2); colorbar;
title('Gaussian 2nd deriv.');

%RF with mag dispersion lines
figure; 
subplot(1,2,1); imagesc(f1_map_rf(1:4:799, 1:4:799)); axis image; caxis(clims1); colorbar;
subplot(1,2,2); imagesc(f2_map_rf(1:4:799, 1:4:799)); axis image; caxis(clims2); colorbar;
title('RF, mag dis lines');
figure; 
subplot(1,2,1); imagesc(f1_map_rt(1:4:799, 1:4:799)); axis image; caxis(clims1); colorbar;
subplot(1,2,2); imagesc(f2_map_rt(1:4:799, 1:4:799)); axis image; caxis(clims2); colorbar;
title('RF thin, mag dis lines');

%RF with mag dispersion lines, angles mod pi
figure; 
subplot(1,2,1); imagesc(f1_map_rf_mp(1:4:799, 1:4:799)); axis image; caxis(clims1); colorbar;
subplot(1,2,2); imagesc(f2_map_rf_mp(1:4:799, 1:4:799)); axis image; caxis(clims2); colorbar;
title('RF, mag dis lines, angles mod pi');
figure; 
subplot(1,2,1); imagesc(f1_map_rt_mp(1:4:799, 1:4:799)); axis image; caxis(clims1); colorbar;
subplot(1,2,2); imagesc(f2_map_rt_mp(1:4:799, 1:4:799)); axis image; caxis(clims2); colorbar;
title('RF thin, mag dis lines, angles mod pi');

%RF with g2d lines
figure; 
subplot(1,2,1); imagesc(f1_map_rfl(1:4:799, 1:4:799)); axis image; caxis(clims1); colorbar;
subplot(1,2,2); imagesc(f2_map_rfl(1:4:799, 1:4:799)); axis image; caxis(clims2); colorbar;
title('RF, g2d lines');
figure; 
subplot(1,2,1); imagesc(f1_map_rtl(1:4:799, 1:4:799)); axis image; caxis(clims1); colorbar;
subplot(1,2,2); imagesc(f2_map_rtl(1:4:799, 1:4:799)); axis image; caxis(clims2); colorbar;
title('RF thin, g2d lines');

%RF with g2d lines, angles mod pi
figure; 
subplot(1,2,1); imagesc(f1_map_rfl_mp(1:4:799, 1:4:799)); axis image; caxis(clims1); colorbar;
subplot(1,2,2); imagesc(f2_map_rfl_mp(1:4:799, 1:4:799)); axis image; caxis(clims2); colorbar;
title('RF, g2d lines, angles mod pi');
figure; 
subplot(1,2,1); imagesc(f1_map_rtl_mp(1:4:799, 1:4:799)); axis image; caxis(clims1); colorbar;
subplot(1,2,2); imagesc(f2_map_rtl_mp(1:4:799, 1:4:799)); axis image; caxis(clims2); colorbar;
title('RF thin, g2d lines, angles mod pi');
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
px_per_mm = 100/9;
sigma_range = [3.2 4.0 5.4 6.2 7.8];
r_max = 3*round(px_per_mm*sigma_range);
spacing = 8;

o_map_rf = load_uint8(...
    'C:\isbe\asymmetry_project\data\orientation_maps\rf_thin\2004_screening_processed\abnormals\013LML_class.mat');
o_map_g2 = load_uint8(...
    'C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\abnormals\013LML_ori.mat');
l_map_g2 = load_uint8(...
    'C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\abnormals\013LML_lines.mat');
l_map_rf = abs(o_map_rf);

mam_mask = u_load(...
    'C:\isbe\asymmetry_project\data\masks\2004_screening_processed\abnormals\013LML_mask.mat');
pec_mask = u_load(...
    'C:\isbe\asymmetry_project\data\masks_pectoral\2004_screening_processed\abnormals\013LML_mask.mat');
mam_mask = mam_mask & ~pec_mask; clear pec_mask;

%Discard lines from outside the mask
l_map_rf(~mam_mask) = 0;
l_map_g2(~mam_mask) = 0;

%Compute map of pixels to discard (i.e. near the edge of the image
%where artefacts cause strong aligned orientations
[discard_rf] = discard_orientations(angle(o_map_rf), 'mask', mam_mask);  
l_map_rf(discard_rf) = 0;
[discard_g2] = discard_orientations(o_map_g2, 'mask', mam_mask);  
l_map_g2(discard_g2) = 0;

%Compute the Karssemeijer maps
[f1_rf f2_rf mask_out] = karssemeijer_radial_projection_multiscale(...
    l_map_rf, angle(o_map_rf),...
    'r_max', r_max(3),...
    'spacing', spacing,...
    'mask', mam_mask);
%%
%Compute the Karssemeijer maps
[f1_rm f2_rm] = karssemeijer_radial_projection_multiscale(...
    mam_mask & ~discard_rf, angle(o_map_rf),...
    'r_max', r_max(3),...
    'spacing', spacing,...
    'mask', mam_mask);

%Compute the Karssemeijer maps
[f1_gm f2_gm] = karssemeijer_radial_projection_multiscale(...
    mam_mask & ~discard_g2, o_map_g2,...
    'r_max', r_max(3),...
    'spacing', spacing,...
    'mask', mam_mask);
%%
%Compute the Karssemeijer maps
[f1_rt f2_rt] = karssemeijer_radial_projection_multiscale(...
    l_map_rf > 0.35, angle(o_map_rf),...
    'r_max', r_max(3),...
    'spacing', spacing,...
    'mask', mam_mask);

%Compute the Karssemeijer maps
[f1_rg f2_rg] = karssemeijer_radial_projection_multiscale(...
    l_map_g2, angle(o_map_rf),...
    'r_max', r_max(3),...
    'spacing', spacing,...
    'mask', mam_mask);
%%
%Compute the Karssemeijer maps
[f1_g2 f2_g2] = karssemeijer_radial_projection_multiscale(...
    l_map_g2, o_map_g2,...
    'r_max', r_max(3),...
    'spacing', spacing,...
    'mask', mam_mask);
[f1_gf f2_gf] = karssemeijer_radial_projection_multiscale(...
    l_map_rf > 0.35, o_map_g2,...
    'r_max', r_max(3),...
    'spacing', spacing,...
    'mask', mam_mask);
%%    
f1_map_rf = zeros(size(mask_out)); 
f1_map_rf(mask_out) = f1_rf; f1_map_rf = f1_map_rf(1:spacing:end,1:spacing:end);
f2_map_rf = zeros(size(mask_out)); 
f2_map_rf(mask_out) = f2_rf; f2_map_rf = f2_map_rf(1:spacing:end,1:spacing:end);
%
f1_map_rm = zeros(size(mask_out)); 
f1_map_rm(mask_out) = f1_rm; f1_map_rm = f1_map_rm(1:spacing:end,1:spacing:end);
f2_map_rm = zeros(size(mask_out)); 
f2_map_rm(mask_out) = f2_rm; f2_map_rm = f2_map_rm(1:spacing:end,1:spacing:end);
%
f1_map_gm = zeros(size(mask_out)); 
f1_map_gm(mask_out) = f1_gm; f1_map_gm = f1_map_gm(1:spacing:end,1:spacing:end);
f2_map_gm = zeros(size(mask_out)); 
f2_map_gm(mask_out) = f2_gm; f2_map_gm = f2_map_gm(1:spacing:end,1:spacing:end);
%%
f1_map_rt = zeros(size(mask_out)); 
f1_map_rt(mask_out) = f1_rt; f1_map_rt = f1_map_rt(1:spacing:end,1:spacing:end);
f2_map_rt = zeros(size(mask_out)); 
f2_map_rt(mask_out) = f2_rt; f2_map_rt = f2_map_rt(1:spacing:end,1:spacing:end);
%
f1_map_rg = zeros(size(mask_out)); 
f1_map_rg(mask_out) = f1_rg; f1_map_rg = f1_map_rg(1:spacing:end,1:spacing:end);
f2_map_rg = zeros(size(mask_out)); 
f2_map_rg(mask_out) = f2_rg; f2_map_rg = f2_map_rg(1:spacing:end,1:spacing:end);
%
%%
f1_map_g2 = zeros(size(mask_out)); 
f1_map_g2(mask_out) = f1_g2; f1_map_g2 = f1_map_g2(1:spacing:end,1:spacing:end);
f2_map_g2 = zeros(size(mask_out)); 
f2_map_g2(mask_out) = f2_g2; f2_map_g2 = f2_map_g2(1:spacing:end,1:spacing:end);
%
f1_map_gf = zeros(size(mask_out)); 
f1_map_gf(mask_out) = f1_gf; f1_map_gf = f1_map_gf(1:spacing:end,1:spacing:end);
f2_map_gf = zeros(size(mask_out)); 
f2_map_gf(mask_out) = f2_gf; f2_map_gf = f2_map_gf(1:spacing:end,1:spacing:end);

%
%%
clims1 = [...
    min([f1_g2(:); f1_rf(:); f1_rt(:); f1_rg(:); f1_gf(:); f1_rm(:); f1_gm(:)])...
    max([f1_g2(:); f1_rf(:); f1_rt(:); f1_rg(:); f1_gf(:); f1_rm(:); f1_gm(:)])];
clims2 = [...
    min([f2_g2(:); f2_rf(:); f2_rt(:); f2_rg(:); f2_gf(:); f2_rm(:); f2_gm(:)])...
    max([f2_g2(:); f2_rf(:); f2_rt(:); f2_rg(:); f2_gf(:); f2_rm(:); f2_gm(:)])];

figure; 
subplot(1,2,1); imagesc(f1_map_rf); axis image; caxis(clims1); colorbar('east');
title('RF line map (continuous)');
subplot(1,2,2); imagesc(f2_map_rf); axis image; caxis(clims2); colorbar('east');
title('RF orientation map');

figure; 
subplot(1,2,1); imagesc(f1_map_rt); axis image; caxis(clims1); colorbar('east');
title('RF line map (thresholded)');
subplot(1,2,2); imagesc(f2_map_rt); axis image; caxis(clims2); colorbar('east');
title('RF orientation map');

figure; 
subplot(1,2,1); imagesc(f1_map_rm); axis image; caxis(clims1); colorbar('east');
title('Uniform line map');
subplot(1,2,2); imagesc(f2_map_rm); axis image; caxis(clims2); colorbar('east');
title('RF orientation map');

figure; 
subplot(1,2,1); imagesc(f1_map_gm); axis image; caxis(clims1); colorbar('east');
title('Uniform line map');
subplot(1,2,2); imagesc(f2_map_gm); axis image; caxis(clims2); colorbar('east');
title('Gaussian orientation map');

figure; 
subplot(1,2,1); imagesc(f1_map_rg); axis image; caxis(clims1); colorbar('east');
title('Gaussian binary line map');
subplot(1,2,2); imagesc(f2_map_rg); axis image; caxis(clims2); colorbar('east');
title('RF orientation map');

figure; 
subplot(1,2,1); imagesc(f1_map_g2); axis image; caxis(clims1); colorbar('east');
title('Gaussian binary line map');
subplot(1,2,2); imagesc(f2_map_g2); axis image; caxis(clims2); colorbar('east');
title('Gaussian orientation map');

figure; 
subplot(1,2,1); imagesc(f1_map_gf); axis image; caxis(clims1); colorbar('east');
title('RF line map (continuous)');
subplot(1,2,2); imagesc(f2_map_gf); axis image; caxis(clims2); colorbar('east');
title('Gaussian orientation map');
%%
mam = u_load(...
    'C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals_half\013LML.mat');
figure; imagesc(mam); axis image; colormap(gray(256));
figure;
subplot(1,2,1); imagesc(l_map_rf); axis image; colormap(gray(256));
subplot(1,2,2); imagesc(l_map_g2); axis image; colormap(gray(256));

