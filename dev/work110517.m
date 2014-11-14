g2_tp = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_g2d_uni_max_scale_f1_abnormals_h.mat', 'tp');
rf_tp = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_rf_uni_max_scale_f1_abnormals_h.mat', 'tp');

g2_fp = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_g2d_uni_max_scale_f1_normals_h.mat', 'fp');
rf_fp = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_rf_uni_max_scale_f1_normals_h.mat', 'fp');
%%

figure; imagesc(g2_tp.tp);
figure; imagesc(rf_tp.tp);
figure; imagesc(g2_tp.tp - rf_tp.tp);

rf_sum = sum(rf_tp.tp > 0, 2);
g2_sum = sum(g2_tp.tp > 0, 2);

[comparison r_idx] = sortrows([g2_sum, rf_sum]);
figure; bar(1:145, comparison);

sum_diff = g2_sum - rf_sum;
[sorted_diff m_idx] = sort(sum_diff);
figure; bar(1:145, sorted_diff);

%Masses better under RF:
rf_wins = m_idx(sorted_diff < 0);
g2_wins = m_idx(sorted_diff > 0);
%%
abnormal_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals.mat');
%Lets look at each type:
for ii = 1:10%length(rf_wins)
    %roi = load_uint8(['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\mass_roi\'...
    %    abnormal_names{rf_wins(ii)} '_roi.mat']);
    %figure; imagesc(roi); axis image; colormap(gray(256));
    
    g2_ori = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\mass_roi\'...
        abnormal_names{rf_wins(ii)} '_roi.mat']);
    g2_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\mass_roi\'...
        abnormal_names{rf_wins(ii)} '_roi.mat']);
    rf_ori = angle(load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_prob\2004_screening_processed\mass_roi\'...
        abnormal_names{rf_wins(ii)} '_roi.mat']));
    rf_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\rf_prob\2004_screening_processed\mass_roi\'...
        abnormal_names{rf_wins(ii)} '_roi.mat']);
    
    figure; 
    subplot(1,2,1); image(complex2rgb(g2_line .* exp(2*i*g2_ori))); axis image;
    subplot(1,2,2); image(complex2rgb(rf_line .* exp(2*i*rf_ori))); axis image;
end
%%
for ii = 1:length(g2_wins)
    roi = load_uint8(['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\mass_roi\'...
       abnormal_names{g2_wins(ii)} '_roi.mat']);
    %figure; imagesc(roi); axis image; colormap(gray(256));
    g2_ori = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\mass_roi\'...
        abnormal_names{g2_wins(ii)} '_roi.mat']);
    g2_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\mass_roi\'...
        abnormal_names{g2_wins(ii)} '_roi.mat']);
    rf_ori = angle(load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_prob\2004_screening_processed\mass_roi\'...
        abnormal_names{g2_wins(ii)} '_roi.mat']));
    rf_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\rf_prob\2004_screening_processed\mass_roi\'...
        abnormal_names{g2_wins(ii)} '_roi.mat']);
    
    figure; 
    a1 = subplot(2,2,1:2); imagesc(roi); axis image; colormap(gray(256));
    a2 = subplot(2,2,3); image(complex2rgb(g2_line .* exp(2*i*g2_ori))); axis image;
    a3 = subplot(2,2,4); image(complex2rgb(rf_line .* exp(2*i*rf_ori))); axis image;
    linkaxes([a1 a2 a3]);
end
%%
for ii = 50:length(g2_wins)
    roi = load_uint8(['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\mass_roi\'...
       abnormal_names{g2_wins(ii)} '_roi.mat']);
    %figure; imagesc(roi); axis image; colormap(gray(256));
    g2_ori = mod(load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\mass_roi\'...
        abnormal_names{g2_wins(ii)} '_roi.mat']),pi);
    g2_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\mass_roi\'...
        abnormal_names{g2_wins(ii)} '_roi.mat']);
    rf_ori = mod(angle(load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_prob\2004_screening_processed\mass_roi\'...
        abnormal_names{g2_wins(ii)} '_roi.mat'])),pi);
    rf_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\rf_prob\2004_screening_processed\mass_roi\'...
        abnormal_names{g2_wins(ii)} '_roi.mat']);
    
%     g2_thin = bwmorph(g2_line, 'skel', 'inf');
%     [y_g2 x_g2] = find(g2_thin);
%     
%     rf_thin = bwmorph(rf_line > 0.3, 'skel', 'inf');
%     [y_rf x_rf] = find(rf_thin);
%     
%     u_g2 = 4*cos(g2_ori(g2_thin));
%     v_g2 = -4*sin(g2_ori(g2_thin));
%     
%     u_rf = 4*cos(rf_ori(g2_thin));
%     v_rf = -4*sin(rf_ori(g2_thin));
%     
%     figure; 
%     a1 = subplot(1,2,1);
%     imagesc(roi); axis image; colormap(gray(256)); hold on;
%     quiver_orientations(g2_ori, 'line_map', g2_line, 'spacing', 2);
%     a2 = subplot(1,2,2);
%     imagesc(roi); axis image; colormap(gray(256)); hold on;
%     quiver_orientations(rf_ori, 'line_map', rf_line>0.3, 'spacing', 2);
%     linkaxes([a1 a2]);
%     figure; 
%     imagesc(roi); axis image; colormap(gray(256)); hold on;
%     quiver(x_g2 - u_g2/2, y_g2 - v_g2/2, u_g2, v_g2, 0, 'r', 'ShowArrowHead', 'off');
%     quiver(x_g2 - u_rf/2, y_g2 - v_rf/2, u_rf, v_rf, 0, 'g', 'ShowArrowHead', 'off');
    
    figure; 
    a1 = subplot(1,2,1);
    imagesc(g2_ori*2); axis image; caxis([0 2*pi]); colormap(hsv(180));
    a2 = subplot(1,2,2);
    imagesc(rf_ori*2); axis image; caxis([0 2*pi]); colormap(hsv(180));
    linkaxes([a1 a2]);
end
%%
for ii = 1:10
    roi = load_uint8(['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\mass_roi\'...
       abnormal_names{rf_wins(ii)} '_roi.mat']);
    %figure; imagesc(roi); axis image; colormap(gray(256));
    g2_ori = mod(load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\mass_roi\'...
        abnormal_names{rf_wins(ii)} '_roi.mat']),pi);
    g2_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\mass_roi\'...
        abnormal_names{rf_wins(ii)} '_roi.mat']);
    rf_ori = mod(angle(load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_prob\2004_screening_processed\mass_roi\'...
        abnormal_names{rf_wins(ii)} '_roi.mat'])),pi);
    rf_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\rf_prob\2004_screening_processed\mass_roi\'...
        abnormal_names{rf_wins(ii)} '_roi.mat']);
    
    figure; 
    a1 = subplot(1,2,1);
    imagesc(g2_ori*2); axis image; caxis([0 2*pi]); colormap(hsv(180));
    a2 = subplot(1,2,2);
    imagesc(rf_ori*2); axis image; caxis([0 2*pi]); colormap(hsv(180));
    linkaxes([a1 a2]);
end
%%
px_per_mm = 100/9;
sigma_range = [3.2 4.0 5.4 6.2 7.8];
r_max = 3*round(px_per_mm*sigma_range);

mask = false(800);
mask(100:700, 100:700) = 1;
%%
%RF with mag dispersion lines
[f1_g2 f2_g2 mask_out] = karssemeijer_radial_projection_multiscale(...
    g2_line, g2_ori, 'mask', mask, 'spacing', 4, 'r_max', r_max(3));
%RF with mag dispersion lines
[f1_rf f2_rf] = karssemeijer_radial_projection_multiscale(...
    rf_line, rf_ori, 'mask', mask, 'spacing', 4, 'r_max', r_max(3));
%%
%RF
[f1_fg f2_fg] = karssemeijer_radial_projection_multiscale(...
    rf_line, g2_ori, 'mask', mask, 'spacing', 4, 'r_max', r_max(3));
%RF
[f1_gf f2_gf] = karssemeijer_radial_projection_multiscale(...
    g2_line, rf_ori, 'mask', mask, 'spacing', 4, 'r_max', r_max(3));
%%
%RF with mag dispersion lines
[f1_ga f2_ga] = karssemeijer_radial_projection_multiscale(...
    g2_line, g2_ori, 'mask', mask, 'spacing', 4, 'r_max', r_max(3), 'r_min', 90);
%RF with mag dispersion lines
[f1_ra f2_ra] = karssemeijer_radial_projection_multiscale(...
    rf_line, rf_ori, 'mask', mask, 'spacing', 4, 'r_max', r_max(3), 'r_min', 90);
%%
f1_map_rf = zeros(size(mask_out));
f1_map_rf(mask_out) = f1_rf; f1_map_rf = f1_map_rf(1:4:end, 1:4:end);

f2_map_rf = zeros(size(mask_out));
f2_map_rf(mask_out) = f2_rf; f2_map_rf = f2_map_rf(1:4:end, 1:4:end);

f1_map_g2 = zeros(size(mask_out));
f1_map_g2(mask_out) = f1_g2; f1_map_g2 = f1_map_g2(1:4:end, 1:4:end);

f2_map_g2 = zeros(size(mask_out));
f2_map_g2(mask_out) = f2_g2; f2_map_g2 = f2_map_g2(1:4:end, 1:4:end);
%%
f1_map_ra = zeros(size(mask_out));
f1_map_ra(mask_out) = f1_ra; f1_map_ra = f1_map_ra(1:4:end, 1:4:end);

f2_map_ra = zeros(size(mask_out));
f2_map_ra(mask_out) = f2_ra; f2_map_ra = f2_map_ra(1:4:end, 1:4:end);

f1_map_ga = zeros(size(mask_out));
f1_map_ga(mask_out) = f1_ga; f1_map_ga = f1_map_ga(1:4:end, 1:4:end);

f2_map_ga = zeros(size(mask_out));
f2_map_ga(mask_out) = f2_ga; f2_map_ga = f2_map_ga(1:4:end, 1:4:end);
%%
f1_map_fg = zeros(size(mask_out));
f1_map_fg(mask_out) = f1_fg; f1_map_fg = f1_map_fg(1:4:end, 1:4:end);

f2_map_fg = zeros(size(mask_out));
f2_map_fg(mask_out) = f2_fg; f2_map_fg = f2_map_fg(1:4:end, 1:4:end);

f1_map_gf = zeros(size(mask_out));
f1_map_gf(mask_out) = f1_gf; f1_map_gf = f1_map_gf(1:4:end, 1:4:end);

f2_map_gf = zeros(size(mask_out));
f2_map_gf(mask_out) = f2_gf; f2_map_gf = f2_map_gf(1:4:end, 1:4:end);
%%
figure;
subplot(1,2,1); imagesc(f1_map_rf); axis image;
subplot(1,2,2); imagesc(f2_map_rf); axis image;

figure;
subplot(1,2,1); imagesc(f1_map_g2); axis image;
subplot(1,2,2); imagesc(f2_map_g2); axis image;
%%
figure;
subplot(1,2,1); imagesc(f1_map_ra); axis image;
subplot(1,2,2); imagesc(f2_map_ra); axis image;

figure;
subplot(1,2,1); imagesc(f1_map_ga); axis image;
subplot(1,2,2); imagesc(f2_map_ga); axis image;
%%
figure;
subplot(1,2,1); imagesc(f1_map_fg); axis image;
subplot(1,2,2); imagesc(f2_map_fg); axis image;

figure;
subplot(1,2,1); imagesc(f1_map_gf); axis image;
subplot(1,2,2); imagesc(f2_map_gf); axis image;
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
g2_tp = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_g2d_max_scale_f2_abnormals_h.mat', 'tp');
rf_tp = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_rf_mix_max_scale_f2_abnormals_h.mat', 'tp');

figure; imagesc(g2_tp.tp);
figure; imagesc(rf_tp.tp);
figure; imagesc(g2_tp.tp - rf_tp.tp);

rf_sum = sum(rf_tp.tp > 0, 2);
g2_sum = sum(g2_tp.tp > 0, 2);

[comparison r_idx] = sortrows([g2_sum, rf_sum]);
figure; bar(1:146, comparison);

sum_diff = g2_sum - rf_sum;
[sorted_diff m_idx] = sort(sum_diff);
figure; bar(1:146, sorted_diff);

%Masses better under RF:
rf_wins = m_idx(sorted_diff < 0);
g2_wins = m_idx(sorted_diff > 0);
%%
abnormal_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals.mat');
%Lets look at each type:
for ii = 1:10%length(rf_wins)
    %roi = load_uint8(['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\mass_roi\'...
    %    abnormal_names{rf_wins(ii)} '_roi.mat']);
    %figure; imagesc(roi); axis image; colormap(gray(256));
    
    g2_ori = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\mass_roi\'...
        abnormal_names{rf_wins(ii)} '_roi.mat']);
    g2_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\mass_roi\'...
        abnormal_names{rf_wins(ii)} '_roi.mat']);
    rf_ori = angle(load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_prob\2004_screening_processed\mass_roi\'...
        abnormal_names{rf_wins(ii)} '_roi.mat']));
    rf_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\rf_prob\2004_screening_processed\mass_roi\'...
        abnormal_names{rf_wins(ii)} '_roi.mat']);
    
    figure; 
    subplot(1,2,1); image(complex2rgb(g2_line .* exp(2*i*g2_ori))); axis image;
    subplot(1,2,2); image(complex2rgb(rf_line .* exp(2*i*rf_ori))); axis image;
end
%%
for ii = 62:length(g2_wins)
    roi = load_uint8(['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\mass_roi\'...
       abnormal_names{g2_wins(ii)} '_roi.mat']);
    %figure; imagesc(roi); axis image; colormap(gray(256));
    g2_ori = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\mass_roi\'...
        abnormal_names{g2_wins(ii)} '_roi.mat']);
    g2_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\mass_roi\'...
        abnormal_names{g2_wins(ii)} '_roi.mat']);
    rf_ori = angle(load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_prob\2004_screening_processed\mass_roi\'...
        abnormal_names{g2_wins(ii)} '_roi.mat']));
    rf_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\rf_prob\2004_screening_processed\mass_roi\'...
        abnormal_names{g2_wins(ii)} '_roi.mat']);
    
    figure; 
    a1 = subplot(2,2,1:2); imagesc(roi); axis image; colormap(gray(256));
    a2 = subplot(2,2,3); image(complex2rgb(g2_line .* exp(2*i*g2_ori))); axis image;
    a3 = subplot(2,2,4); image(complex2rgb(rf_line .* exp(2*i*rf_ori))); axis image;
    linkaxes([a1 a2 a3]);
end
%%
figure;
subplot(2,1,1); hold on;
plot(49:-1:1, sum(g2_tp.tp(:,2:end) > 0), 'r');
plot(49:-1:1, sum(rf_tp.tp(:,2:end) > 0), 'b');

subplot(2,1,2); hold on;
plot(49:-1:1, sum(g2_fp.fp(:,2:end)), 'r');
plot(49:-1:1, sum(rf_fp.fp(:,2:end)), 'b');

