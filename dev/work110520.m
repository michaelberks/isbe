px_per_mm = 100/9;
sigma_range = [3.2 4.0 5.4 6.2 7.8];
r_max = 3*round(px_per_mm*sigma_range);
r_min = px_per_mm*4;
R = px_per_mm*2;

xy = repmat(-r_max(end):r_max(end), 2*r_max(end)+1, 1);

rij = (2 * R) ./ (pi*sqrt(xy.^2 + xy'.^2));
rij(r_max(end)+1, r_max(end)+1) = 0;

pN = zeros(5,1);
sqrt_Npq = zeros(5,1);

for ii = 1:5
    N_mask = (xy.^2 + xy'.^2 < r_max(ii)^2) & (xy.^2 + xy'.^2 > r_min^2);
    N = sum(N_mask(:));
    p = sum(rij(N_mask)) / N;
    
    pN(ii) = p* N;
    sqrt_Npq(ii) = sqrt(N*p*(1-p));
end

figure; 
subplot(1,2,1); plot(r_max, pN);
ylabel('pN');
xlabel('r_{max}');
subplot(1,2,2); plot(r_max, sqrt_Npq);
ylabel('$\sqrt{Np(1-p)}$', 'interpreter', 'latex');
xlabel('r_{max}');
%%
g2_uni = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_g2d_uni_max_scale_f1_abnormals_h.mat', 'tp');
s = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_g2d_uni_max_scale_f1_normals_h.mat', 'fp');
g2_uni.fp = s.fp; clear s;

rf_uni = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_rf_uni_max_scale_f1_abnormals_h.mat', 'tp');
s = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_rf_uni_max_scale_f1_normals_h.mat', 'fp');
rf_uni.fp = s.fp; clear s;

g2_pro = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_g2d_max_scale_f1_abnormals_h.mat', 'tp');
s = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_g2d_max_scale_f1_normals_h.mat', 'fp');
g2_pro.fp = s.fp; clear s;

rf_pro = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_rf_prob_max_scale_f1_abnormals_h.mat', 'tp');
s = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_rf_prob_max_scale_f1_normals_h.mat', 'fp');
rf_pro.fp = s.fp; clear s;

rf_thi = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_rf_thin_max_scale_f1_abnormals_h.mat', 'tp');
s = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_rf_thin_max_scale_f1_normals_h.mat', 'fp');
rf_thi.fp = s.fp; clear s;

rf_thr = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_rf_thresh_max_scale_f1_abnormals_h.mat', 'tp');
s = load('C:\isbe\asymmetry_project\experiments\k_maps\froc_rf_thresh_max_scale_f1_normals_h.mat', 'fp');
rf_thr.fp = s.fp; clear s;
%%
legend_text = {'g2d uniform', 'g2d line', 'rf uniform', 'rf prob', 'rf thin', 'rf thresh'};
figure;
subplot(2,1,1); hold all;
plot(49:-1:1, sum(g2_uni.tp(:,2:end) > 0));
plot(49:-1:1, sum(g2_pro.tp(:,2:end) > 0));
plot(49:-1:1, sum(rf_uni.tp(:,2:end) > 0));
plot(49:-1:1, sum(rf_pro.tp(:,2:end) > 0));
plot(49:-1:1, sum(rf_thi.tp(:,2:end) > 0));
plot(49:-1:1, sum(rf_thr.tp(:,2:end) > 0));
legend(legend_text, 'location', 'southeast')

subplot(2,1,2); hold all;
plot(49:-1:1, mean(g2_uni.fp(:,2:end)));
plot(49:-1:1, mean(g2_pro.fp(:,2:end)));
plot(49:-1:1, mean(rf_uni.fp(:,2:end)));
plot(49:-1:1, mean(rf_pro.fp(:,2:end)));
plot(49:-1:1, mean(rf_thi.fp(:,2:end)));
plot(49:-1:1, mean(rf_thr.fp(:,2:end)));
legend(legend_text, 'location', 'southeast')

%%
g2_uni.tp = [g2_uni.tp(1:65,:); g2_pro.tp(66,:); g2_uni.tp(66:end,:)]; 
rf_uni_sum = sum(rf_uni.tp > 0, 2);
g2_uni_sum = sum(g2_uni.tp > 0, 2);
rf_pro_sum = sum(rf_pro.tp > 0, 2);
g2_pro_sum = sum(g2_pro.tp > 0, 2);
rf_thi_sum = sum(rf_thi.tp > 0, 2);
rf_thr_sum = sum(rf_thr.tp > 0, 2);

rf_uni_diff = rf_uni_sum - rf_pro_sum;
rf_uni_wins = rf_uni_diff > 0;
rf_pro_wins = rf_uni_diff < 0;

g2_uni_diff = g2_uni_sum - g2_pro_sum;
g2_uni_wins = g2_uni_diff > 0;
g2_pro_wins = g2_uni_diff < 0;

gf_uni_diff = g2_uni_sum - rf_uni_sum;
gf_pro_diff = g2_pro_sum - rf_pro_sum;
gf_thr_diff = g2_pro_sum - rf_thr_sum;

%%
[mm m_idx] = max(rf_uni_diff);
abnormal_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals.mat');
display(abnormal_names{m_idx});
%%
map_type = {'g2d', 'rf_thin', 'rf_thresh', 'rf_prob', 'rf_mix', 'rf_uni', 'g2d_uni'};
ab_type = {'normals', 'abnormals'};
offset_r = 0;
offset_c = 0;

aa = 2;
    
%load in the lists of mammograms names
mam_names = u_load(['C:\isbe\asymmetry_project\data\mam_names\2004_screening_' ab_type{aa}]);
    
for mm = [6 7] %For each method rf, wrf, g2
    if mm > 5
        spacing = 4;
    else
        spacing = 2;
    end

    %create a directory for the combined maps
    mkdir(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
        '_template_scale\2004_screening_processed\' ab_type{aa} '\']);
    mkdir(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
        '_max_scale\2004_screening_processed\' ab_type{aa} '\']);

    for ii = m_idx;%1:length(mam_names) %for each mammogram

        if mm == 7 && ii == 66; continue; end
        %load in the f1 and f2 maps
        f1 = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
            '_max_scale\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_f1.mat']);
        f2 = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
            '_max_scale\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_f2.mat']);

        f1_scale = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
            '_max_scale\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_f1_scale.mat']);
        f2_scale = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
            '_max_scale\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_f2_scale.mat']);

        mask = u_load(['C:\isbe\asymmetry_project\data\k_stellate_maps\' map_type{mm}...
            '\2004_screening_processed\' ab_type{aa} '\' mam_names{ii} '_mask.mat']);
        mask = mask(offset_r + (1:spacing:end), offset_c + (1:spacing:end));

        f1_map = zeros(size(mask));
        f2_map = zeros(size(mask));
        f1_map(mask) = f1;
        f2_map(mask) = f2;   

        f1_scale_map = zeros(size(mask));
        f2_scale_map = zeros(size(mask));
        f1_scale_map(mask) = f1_scale;
        f2_scale_map(mask) = f2_scale; 

        %load in mass border
        mass_xy = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\'...
            mam_names{ii} '_meta.mat']);
        [r c] = size(f1_map);
        x = mass_xy(:,1) * c;
        y = mass_xy(:,2) * r;

        xc = round(mean(x));
        yc = round(mean(y));

        figure; 
        subplot(1,2,1); imagesc(f1_map); axis image; hold on;
        plot(x, y, 'k'); plot(xc, yc, 'kx');
        subplot(1,2,2); imagesc(f2_map); axis image; hold on;
        plot(x, y, 'k'); plot(xc, yc, 'kx');

        figure; 
        subplot(1,2,1); imagesc(f1_scale_map); axis image; colormap(lines(6));
        colorbar('location', 'west'); hold on;
        plot(x, y, 'k'); plot(xc, yc, 'kx');
        subplot(1,2,2); imagesc(f2_scale_map); axis image; colormap(lines(6)); 
        colorbar('location', 'west'); hold on;
        plot(x, y, 'k'); plot(xc, yc, 'kx');

    end
end
%%
%%
mam = load_uint8(['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\abnormals\'...
    abnormal_names{m_idx} '.mat']);
mass_xy = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\'...
     abnormal_names{m_idx} '_meta.mat']);

[r c] = size(mam);
x = mass_xy(:,1) * c;
y = mass_xy(:,2) * r;

xc = round(mean(x)); %= 4*369
yc = round(mean(y)); %= 4*268

roi = sample_window(mam, 800, yc, xc, 0);

g2_ori = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\abnormals\'...
    abnormal_names{m_idx} '_ori.mat']);
g2_ori = sample_window(g2_ori, 800, yc, xc, 0);

g2_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\abnormals\'...
    abnormal_names{m_idx} '_lines.mat']);
g2_line = sample_window(g2_line, 800, yc, xc, 0);

rf_ori = mod(angle(load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_prob\2004_screening_processed\abnormals\'...
    abnormal_names{m_idx} '_class.mat'])),pi);
rf_ori = sample_window(rf_ori, 800, yc, xc, 0);

rf_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\rf_prob\2004_screening_processed\abnormals\'...
    abnormal_names{m_idx} '_class.mat']);
rf_line = sample_window(rf_line, 800, yc, xc, 0);
pack;
%%
figure; imagesc(g2_line); axis image; colormap(gray(256));
figure; imagesc(rf_line > 0.5); axis image; colormap(gray(256));
figure; imagesc(2*g2_ori); axis image; caxis([0 2*pi]); colormap(hsv(256));
figure; imagesc(2*rf_ori); axis image; caxis([0 2*pi]); colormap(hsv(256));
%%
px_per_mm = 100/9;
sigma_range = [3.2 4.0 5.4 6.2 7.8];
r_max = 3*round(px_per_mm*sigma_range);
r_min = px_per_mm*4;
R = px_per_mm*2;

xx = repmat(-400:399, 800, 1);
yy = repmat((399:-1:-400)', 1, 800);

theta = mod(atan2(yy, xx), pi);
rij = sqrt(xx.^2 + yy.^2);
phi = abs(asin(R./rij));

g2_theta_diff = abs(theta - g2_ori);
g2_theta_mask = ((g2_theta_diff < phi) | (abs(g2_theta_diff-pi) < phi)) & (rij > r_min);

rf_theta_diff = abs(theta - rf_ori);
rf_theta_mask = ((rf_theta_diff < phi) | (abs(rf_theta_diff-pi) < phi)) & (rij > r_min);

figure; imagesc(g2_theta_mask); axis image; colormap(gray(256));
figure; imagesc(rf_theta_mask); axis image; colormap(gray(256));
%

rf_uni_r = double(rf_line > 0.5 & (rij > r_min));
rf_uni_g = double(rf_theta_mask);
rf_uni_rgb = cat(3, rf_uni_r, rf_uni_g, zeros(800,800));

g2_uni_r = double(g2_line & (rij > r_min));
g2_uni_g = double(g2_theta_mask);
g2_uni_rgb = cat(3, g2_uni_r, g2_uni_g, zeros(800,800));

figure; image(rf_uni_rgb); axis image
figure; image(g2_uni_rgb); axis image
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
px_per_mm = 100/9;
sigma_range = [3.2 4.0 5.4 6.2 7.8];
r_max = 3*round(px_per_mm*sigma_range);
r_min = px_per_mm*4;
R = px_per_mm*2;

xx = repmat(-400:399, 800, 1);
yy = repmat((399:-1:-400)', 1, 800);

theta = mod(atan2(yy, xx), pi);
rij = sqrt(xx.^2 + yy.^2);
phi = abs(asin(R./rij));

pij = 2*R./(pi*rij);

offset_r = 0;
offset_c = 0;
spacing = 2;

r_sums = zeros(5,1);
r_mask = false(800,800,5);
for rr = 1:length(r_max)
    r_mask(:,:,rr) = (rij < r_max(rr)) & (rij > r_min);
    r_sums(rr) = sum(sum(r_mask(:,:,rr)));
end

%load in the lists of mammograms names
abnormal_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals.mat');
%%
[dummy sort_idx] = sort(gf_thr_diff, 'descend');
%%
for kk = 136:146;
 
    ii = kk;% sort_idx(kk);
    
    %load in the f1 map and mask
    f1 = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\g2d_max_scale'...
        '\2004_screening_processed\abnormals\' abnormal_names{ii} '_f1.mat']);
    mask = u_load(['C:\isbe\asymmetry_project\data\k_stellate_maps\g2d_max_scale'...
        '\2004_screening_processed\abnormals\' abnormal_names{ii} '_mask.mat']);

    %Convert f1 vector into map and smooth
    mask = mask(offset_r + (1:spacing:end), offset_c + (1:spacing:end));
    f1_map = zeros(size(mask));
    f1_map(mask) = f1;
    f1_map = imfilter(f1_map, fspecial('gaussian', 21, 4));
    
    %load in mass border
    mass_xy = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\'...
        abnormal_names{ii} '_meta.mat']);
    [r c] = size(f1_map);
    x = mass_xy(:,1) * c;
    y = mass_xy(:,2) * r;

    %Workout location of maximum f1 score within mass border
    mass_mask = poly2mask(x,y,r,c);
    [yc xc] = find(mass_mask);
    [dummy mass_max_idx] = max(f1_map(mass_mask));
    xc = spacing*xc(mass_max_idx);
    yc = spacing*yc(mass_max_idx);
    
%     figure; imagesc(f1_map); axis image; hold on;
%     plot(x, y, 'k');
%     plot(xc / spacing, yc / spacing, 'kx');
    clear f1_map mask; pack;

    %Load in regions associated with mass
    roi = load_uint8(['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '.mat']);
    roi = sample_window(roi, 800, yc, xc, 0);

    g2_ori = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_ori.mat']);
    g2_ori = sample_window(g2_ori, 800, yc, xc, 0);

    g2_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_lines.mat']);
    g2_line = sample_window(g2_line, 800, yc, xc, 0);

    rf_ori = mod(angle(load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_prob\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_class.mat'])),pi);
    rf_ori = sample_window(rf_ori, 800, yc, xc, 0);

    rf_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\rf_prob\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_class.mat']);
    rf_line = sample_window(rf_line, 800, yc, xc, 0);
    
    rf_pm = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_thin\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_class.mat']);
    rf_pm = sample_window(rf_pm, 800, yc, xc, 0);

    rf_ori2 = mod(angle(rf_pm),pi);
    rf_line2 = abs(rf_pm);
    
    pack;
    
    %Look at target pixels
    g2_theta_diff = abs(theta - g2_ori);
    g2_theta_mask = ((g2_theta_diff < phi) | (abs(g2_theta_diff-pi) < phi)) & (rij > r_min);

    rf_theta_diff = abs(theta - rf_ori);
    rf_theta_mask = ((rf_theta_diff < phi) | (abs(rf_theta_diff-pi) < phi)) & (rij > r_min);
    
    rf_theta_diff2 = abs(theta - rf_ori2);
    rf_theta_mask2 = ((rf_theta_diff2 < phi) | (abs(rf_theta_diff2-pi) < phi)) & (rij > r_min);

    %
    rf_lm = rf_line > 0.5 & (rij > r_min);
    rf_uni_r = double(rf_lm);
    rf_uni_g = double((rf_theta_mask & ~rf_lm) | (~rf_theta_mask & rf_lm));
    rf_uni_b = double(~rf_theta_mask & rf_lm);
    rf_uni_rgb = cat(3, rf_uni_r, rf_uni_g, rf_uni_b);
    
    rf_lm2 = rf_line2 > 0.35 & (rij > r_min);
    rf_uni_r2 = double(rf_lm2);
    rf_uni_g2 = double((rf_theta_mask2 & ~rf_lm2) | (~rf_theta_mask2 & rf_lm2));
    rf_uni_b2 = double(~rf_theta_mask2 & rf_lm2);
    rf_uni_rgb2 = cat(3, rf_uni_r2, rf_uni_g2, rf_uni_b2);

    g2_lm = g2_line & (rij > r_min);
    g2_uni_r = double(g2_lm);
    g2_uni_g = double((g2_theta_mask & ~g2_lm) | (~g2_theta_mask & g2_lm));
    g2_uni_b = double(~g2_theta_mask & g2_lm);
    g2_uni_rgb = cat(3, g2_uni_r, g2_uni_g, g2_uni_b);

    alpha = linspace(0, 2*pi, 200);
    
    figure; 
    a1 = subplot(2,2,1); imagesc(roi); axis image; colormap(gray(256));
    
    a2 = subplot(2,2,2); image(rf_uni_rgb); axis image; hold on;
    rf_u_scores = zeros(1,5);
    rf_n_scores = zeros(1,5);
    rf_f_scores = zeros(1,5);
    for rr = 1:length(r_max)
        plot(r_max(rr)*cos(alpha)+399, r_max(rr)*sin(alpha)+399);
        rf_u_scores(rr) = 100*sum(rf_theta_mask(r_mask(:,:,rr))) / r_sums(rr);
        rf_n_scores(rr) = sum(rf_theta_mask(r_mask(:,:,rr)) & rf_lm(r_mask(:,:,rr)));
        rf_f_scores(rr) = 100*sum(rf_theta_mask(r_mask(:,:,rr)) & rf_lm(r_mask(:,:,rr)))/...
            sum(rf_lm(r_mask(:,:,rr)));
    end
    title(['RF Line score: ' num2str(rf_pro_sum(ii)) ', uniform score: ' num2str(rf_uni_sum(ii))]); 
    xlabel({['U scores: ' num2str(rf_u_scores)];...
        ['N scores: ' num2str(rf_n_scores)];...
        ['f scores: ' num2str(rf_f_scores)]});
    
    a3 = subplot(2,2,3); image(rf_uni_rgb2); axis image; hold on;
    rf_u_scores = zeros(1,5);
    rf_n_scores = zeros(1,5);
    rf_f_scores = zeros(1,5);
    for rr = 1:length(r_max)
        plot(r_max(rr)*cos(alpha)+399, r_max(rr)*sin(alpha)+399);
        rf_u_scores(rr) = 100*sum(rf_theta_mask2(r_mask(:,:,rr))) / r_sums(rr);
        rf_n_scores(rr) = sum(rf_theta_mask2(r_mask(:,:,rr)) & rf_lm2(r_mask(:,:,rr)));
        rf_f_scores(rr) = 100*sum(rf_theta_mask2(r_mask(:,:,rr)) & rf_lm2(r_mask(:,:,rr)))/...
            sum(rf_lm2(r_mask(:,:,rr)));
    end
    title(['RF Line score: ' num2str(rf_thr_sum(ii)) ', uniform score: ' num2str(rf_uni_sum(ii))]);
    xlabel({['U scores: ' num2str(rf_u_scores)];...
        ['N scores: ' num2str(rf_n_scores)];...
        ['f scores: ' num2str(rf_f_scores)]});
    
    
    a4 = subplot(2,2,4); image(g2_uni_rgb); axis image; hold on;
    g2_u_scores = zeros(1,5);
    g2_n_scores = zeros(1,5);
    g2_f_scores = zeros(1,5);
    for rr = 1:length(r_max)
        plot(r_max(rr)*cos(alpha)+399, r_max(rr)*sin(alpha)+399);
        g2_u_scores(rr) = 100*sum(g2_theta_mask(r_mask(:,:,rr))) / r_sums(rr);
        g2_n_scores(rr) = sum(g2_theta_mask(r_mask(:,:,rr)) & g2_lm(r_mask(:,:,rr)));
        g2_f_scores(rr) = 100*sum(g2_theta_mask(r_mask(:,:,rr)) & g2_lm(r_mask(:,:,rr)))/...
            sum(g2_lm(r_mask(:,:,rr)));
    end
    title(['G 2nd Deriv. Line score: ' num2str(g2_pro_sum(ii)) ', uniform score: ' num2str(g2_uni_sum(ii))]);
    xlabel({['U scores: ' num2str(g2_u_scores)];...
        ['N scores: ' num2str(g2_n_scores)];...
        ['f scores: ' num2str(g2_f_scores)]});
    
    linkaxes([a1 a2 a3 a4]);

end
%%
warning('off', 'load_uint8:missing_variables');

px_per_mm = 100/9;
sigma_range = [3.2 4.0 5.4 6.2 7.8];
r_max = 3*round(px_per_mm*sigma_range);
r_min = px_per_mm*4;
R = px_per_mm*2;

xx = repmat(-400:399, 800, 1);
yy = repmat((399:-1:-400)', 1, 800);

theta = mod(atan2(yy, xx), pi);
rij = sqrt(xx.^2 + yy.^2);
phi = abs(asin(R./rij));

pij = 2*R./(pi*rij);

offset_r = 0;
offset_c = 0;
spacing = 2;

r_sums = zeros(5,1);
p_sums = zeros(5,1);
r_mask = false(800,800,5);
for rr = 1:length(r_max)
    r_mask(:,:,rr) = (rij < r_max(rr)) & (rij > r_min);
    r_sums(rr) = sum(sum(r_mask(:,:,rr)));
    p_sums(rr) = sum(pij(r_mask(:,:,rr))) / r_sums(rr);
end

%load in the lists of mammograms names
abnormal_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals.mat');

rp_u_scores = zeros(146,5);
rp_n_scores = zeros(146,5);
rp_N_scores = zeros(146,5);
rp_p_scores = zeros(146,5);

rt_u_scores = zeros(146,5);
rt_n_scores = zeros(146,5);
rt_N_scores = zeros(146,5);
rt_p_scores = zeros(146,5);

g2_u_scores = zeros(146,5);
g2_n_scores = zeros(146,5);
g2_N_scores = zeros(146,5);
g2_p_scores = zeros(146,5);

for kk = 1:146;
 
    ii = kk;
    display(['processing mammo ' num2str(ii)]);
    
    %load in the f1 map and mask
    f1 = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\g2d_max_scale'...
        '\2004_screening_processed\abnormals\' abnormal_names{ii} '_f1.mat']);
    mask = u_load(['C:\isbe\asymmetry_project\data\k_stellate_maps\g2d_max_scale'...
        '\2004_screening_processed\abnormals\' abnormal_names{ii} '_mask.mat']);

    %Convert f1 vector into map and smooth
    mask = mask(offset_r + (1:spacing:end), offset_c + (1:spacing:end));
    f1_map = zeros(size(mask));
    f1_map(mask) = f1;
    f1_map = imfilter(f1_map, fspecial('gaussian', 21, 4));
    
    %load in mass border
    mass_xy = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\'...
        abnormal_names{ii} '_meta.mat']);
    [r c] = size(f1_map);
    x = mass_xy(:,1) * c;
    y = mass_xy(:,2) * r;

    %Workout location of maximum f1 score within mass border
    mass_mask = poly2mask(x,y,r,c);
    [yc xc] = find(mass_mask);
    [dummy mass_max_idx] = max(f1_map(mass_mask));
    xc = spacing*xc(mass_max_idx);
    yc = spacing*yc(mass_max_idx);
    
    clear f1_map mask; pack;

    %Load in regions associated with mass
    roi = load_uint8(['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '.mat']);
    roi = sample_window(roi, 800, yc, xc, 0);

    g2_ori = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_ori.mat']);
    g2_ori = sample_window(g2_ori, 800, yc, xc, 0);

    g2_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_lines.mat']);
    g2_line = sample_window(g2_line, 800, yc, xc, 0);

    rp_ori = mod(angle(load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_prob\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_class.mat'])),pi);
    rp_ori = sample_window(rp_ori, 800, yc, xc, 0);

    rp_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\rf_prob\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_class.mat']);
    rp_line = sample_window(rp_line, 800, yc, xc, 0);
    
    rf_pm = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_thin\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_class.mat']);
    rf_pm = sample_window(rf_pm, 800, yc, xc, 0);

    rt_ori = mod(angle(rf_pm),pi);
    rt_line = abs(rf_pm);
    
    pack;
    
    %Look at target pixels
    g2_theta_diff = abs(theta - g2_ori);
    g2_theta_mask = ((g2_theta_diff < phi) | (abs(g2_theta_diff-pi) < phi)) & (rij > r_min);

    rp_theta_diff = abs(theta - rp_ori);
    rp_theta_mask = ((rp_theta_diff < phi) | (abs(rp_theta_diff-pi) < phi)) & (rij > r_min);
    
    rt_theta_diff = abs(theta - rt_ori);
    rt_theta_mask = ((rt_theta_diff < phi) | (abs(rt_theta_diff-pi) < phi)) & (rij > r_min);

    %
    rp_lm = rp_line > 0.5 & (rij > r_min);  
    rt_lm = rt_line > 0.35 & (rij > r_min);
    g2_lm = g2_line & (rij > r_min);


    for rr = 1:length(r_max)

        rp_u_scores(ii,rr) = sum(rp_theta_mask(r_mask(:,:,rr)));
        rp_n_scores(ii,rr) = sum(rp_theta_mask(r_mask(:,:,rr)) & rp_lm(r_mask(:,:,rr)));
        rp_N_scores(ii,rr) = sum(rp_lm(r_mask(:,:,rr)));
        rp_p_scores(ii,rr) = sum(rp_lm(r_mask(:,:,rr)) .* pij(r_mask(:,:,rr)));

        rt_u_scores(ii,rr) = sum(rt_theta_mask(r_mask(:,:,rr)));
        rt_n_scores(ii,rr) = sum(rt_theta_mask(r_mask(:,:,rr)) & rt_lm(r_mask(:,:,rr)));
        rt_N_scores(ii,rr) = sum(rt_lm(r_mask(:,:,rr)));
        rt_p_scores(ii,rr) = sum(rt_lm(r_mask(:,:,rr)) .* pij(r_mask(:,:,rr)));


        g2_u_scores(ii,rr) = sum(g2_theta_mask(r_mask(:,:,rr)));
        g2_n_scores(ii,rr) = sum(g2_theta_mask(r_mask(:,:,rr)) & g2_lm(r_mask(:,:,rr)));
        g2_N_scores(ii,rr) = sum(g2_lm(r_mask(:,:,rr)));
        g2_p_scores(ii,rr) = sum(g2_lm(r_mask(:,:,rr)) .* pij(r_mask(:,:,rr)));
    end


end

rp_p_scores = rp_p_scores ./ rp_N_scores;
rt_p_scores = rt_p_scores ./ rt_N_scores;
g2_p_scores = g2_p_scores ./ g2_N_scores;
%%
rt_f1 = (rt_n_scores - rt_N_scores.*rt_p_scores) ./ sqrt(rt_N_scores.*rt_p_scores.*(1-rt_p_scores));
rp_f1 = (rp_n_scores - rp_N_scores.*rp_p_scores) ./ sqrt(rp_N_scores.*rp_p_scores.*(1-rp_p_scores));
g2_f1 = (g2_n_scores - g2_N_scores.*g2_p_scores) ./ sqrt(g2_N_scores.*g2_p_scores.*(1-g2_p_scores));

rt_fu = bsxfun(@rdivide, bsxfun(@minus, rt_u_scores, (r_sums.*p_sums)'), sqrt(r_sums.*p_sums.*(1-p_sums))');
rp_fu = bsxfun(@rdivide, bsxfun(@minus, rp_u_scores, (r_sums.*p_sums)'), sqrt(r_sums.*p_sums.*(1-p_sums))');
g2_fu = bsxfun(@rdivide, bsxfun(@minus, g2_u_scores, (r_sums.*p_sums)'), sqrt(r_sums.*p_sums.*(1-p_sums))');
%%
cdf_rt = load('C:\isbe\asymmetry_project\experiments\k_maps\cdf_rf_thresh_max_scale_f1.mat');
cdf_rp = load('C:\isbe\asymmetry_project\experiments\k_maps\cdf_rf_prob_max_scale_f1.mat');
cdf_g2 = load('C:\isbe\asymmetry_project\experiments\k_maps\cdf_g2d_max_scale_f1.mat');
cdf_gu = load('C:\isbe\asymmetry_project\experiments\k_maps\cdf_g2d_uni_max_scale_f1.mat');
cdf_ru = load('C:\isbe\asymmetry_project\experiments\k_maps\cdf_rf_uni_max_scale_f1.mat');

figure; hold on;
plot(cdf_rt.x, cdf_rt.P_x);
plot(cdf_rp.x, cdf_rp.P_x, 'r');
plot(cdf_g2.x, cdf_g2.P_x, 'g');
plot(cdf_gu.x, cdf_gu.P_x, 'm');
plot(cdf_ru.x, cdf_ru.P_x, 'k');

legend({'RF thresh', 'RF prob', 'Gauss.', 'Gauss. uniform', 'RF uniform'});
%%
aa_type = {'abnormals\', 'normals\'};
for ii = 1:10
    
    %load in the f1 map and mask
    f1 = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\g2d_max_scale'...
        '\2004_screening_processed\abnormals\' abnormal_names{ii} '_f1.mat']);
    mask = u_load(['C:\isbe\asymmetry_project\data\k_stellate_maps\g2d_max_scale'...
        '\2004_screening_processed\abnormals\' abnormal_names{ii} '_mask.mat']);

    %Convert f1 vector into map and smooth
    mask = mask(offset_r + (1:spacing:end), offset_c + (1:spacing:end));
    f1_map = zeros(size(mask));
    f1_map(mask) = f1;
    f1_map = imfilter(f1_map, fspecial('gaussian', 21, 4));
    
    %load in mass border
    mass_xy = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\'...
        abnormal_names{ii} '_meta.mat']);
    [r c] = size(f1_map);
    x = mass_xy(:,1) * c;
    y = mass_xy(:,2) * r;

    %Workout location of maximum f1 score within mass border
    mass_mask = poly2mask(x,y,r,c);
    [yc xc] = find(mass_mask);
    [dummy mass_max_idx] = max(f1_map(mass_mask));
    xc = spacing*xc(mass_max_idx);
    yc = spacing*yc(mass_max_idx);
    
%     figure; imagesc(f1_map); axis image; hold on;
%     plot(x, y, 'k');
%     plot(xc / spacing, yc / spacing, 'kx');
    clear f1_map mask; pack;
    
    mam_name{1} = abnormal_names{ii};
    if abnormal_names{ii}(4) == 'R'
        mam_name{2} = [mam_name{1}(1:3) 'L' mam_name{1}(5:6)];
    else
        mam_name{2} = [mam_name{1}(1:3) 'R' mam_name{1}(5:6)];
    end
    
    figure;
    for aa = 1:2

        if aa == 2
            xc = c - xc;
        end
        %Load in regions associated with mass
        roi = load_uint8(['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\' aa_type{aa}...
             mam_name{aa} '.mat']);
        roi = sample_window(roi, 800, yc, xc, 0);

        g2_ori = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\' aa_type{aa}...
            mam_name{aa} '_ori.mat']);
        g2_ori = sample_window(g2_ori, 800, yc, xc, 0);

        g2_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\' aa_type{aa}...
            mam_name{aa} '_lines.mat']);
        g2_line = sample_window(g2_line, 800, yc, xc, 0);

        rf_ori = mod(angle(load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_prob\2004_screening_processed\' aa_type{aa}...
            mam_name{aa} '_class.mat'])),pi);
        rf_ori = sample_window(rf_ori, 800, yc, xc, 0);

        rf_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\rf_prob\2004_screening_processed\' aa_type{aa}...
            mam_name{aa} '_class.mat']);
        rf_line = sample_window(rf_line, 800, yc, xc, 0);

        rf_pm = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_thin\2004_screening_processed\' aa_type{aa}...
            mam_name{aa} '_class.mat']);
        rf_pm = sample_window(rf_pm, 800, yc, xc, 0);

        rf_ori2 = mod(angle(rf_pm),pi);
        rf_line2 = abs(rf_pm);

        pack;

        %Look at target pixels
        g2_theta_diff = abs(theta - g2_ori);
        g2_theta_mask = ((g2_theta_diff < phi) | (abs(g2_theta_diff-pi) < phi)) & (rij > r_min);

        rf_theta_diff = abs(theta - rf_ori);
        rf_theta_mask = ((rf_theta_diff < phi) | (abs(rf_theta_diff-pi) < phi)) & (rij > r_min);

        rf_theta_diff2 = abs(theta - rf_ori2);
        rf_theta_mask2 = ((rf_theta_diff2 < phi) | (abs(rf_theta_diff2-pi) < phi)) & (rij > r_min);

        roi = 0.8*(roi - min(roi(:)))/(max(roi(:))-min(roi(:)));
        roi_rgb = cat(3, 0.2*g2_theta_mask + roi, 0.2*rf_theta_mask2 + roi, roi);

        subplot(1,2,aa); image(roi_rgb); axis image;
    end
    
    
end
