mammo_names = get_mammo_info(dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\*.mat'));

for jj = 1:10
    mammo = imresize(u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\o04_' mammo_names{jj} '.mat']), 0.5, 'bilinear');
    map = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening\abnormals\' mammo_names{jj} '_template.mat']);
    scales = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening\abnormals\scales\' mammo_names{jj} '_scales.mat']);

    [maxima_pos] = local_image_maxima(map, 80, [], 0.5);
    
    figure; 
    h1 = subplot(1,2,1); imagesc(map); axis image; colormap(jet(256)); hold on;
    plot(maxima_pos(:,1), maxima_pos(:,2), 'g*');
    h2 = subplot(1,2,2); imagesc(mammo); axis image; colormap(jet(256)); hold on;
    plot(maxima_pos(:,1), maxima_pos(:,2), 'g*');
    for ii = 1:size(maxima_pos, 1)
        x = maxima_pos(ii,1);
        y = maxima_pos(ii,2);
        r = double(scales(y,x));
        axes(h1); 
        plot_ellipse(r, r, [1 0], x, y, 'g:');
        plot_ellipse(3*r, 3*r, [1 0], x, y, 'g');
        axes(h2); plot_ellipse(r, r, [1 0], x, y, 'g');
        plot_ellipse(r, r, [1 0], x, y, 'g:');
        plot_ellipse(3*r, 3*r, [1 0], x, y, 'g');
    end
end
%%
ori_map = zeros(101);
line_map = zeros(101); line_map(51, 51) = 1;
xx = repmat(1:101, 101, 1);
yy = flipud(xx');
R = 5;
r_min = 10;
r_max = 30;
f_map = zeros(101);
%%
for xi = 1:101
    for yi = 1:101
        d = sqrt((xx - xi).^2 + (yy - yi).^2);
        
        v = atan((yy - yi) ./ (xx - xi));
        
        w = (abs(ori_map - v) < (R ./ d)) & line_map & (d > r_min) & (d < r_max);
        
        f_map(yi, xi) = sum(w(:));
    end
end
%%
%Baby, make me a spiral
bg = u_load('C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\bg00001.mat');
label = zeros(512);
label_orientation = zeros(512);
for ii = 0:24:168
    [bar, labeli] = create_sin_bar(4, 4, ii, 512, 512, 0, 256, 256);
    bg = bg + bar;
    label = label | labeli;
    label_orientation(labeli) = ii;
end
%
figure; imagesc(label); axis image;
[fi1a fi2a nia Nia prob_ia ni_plusa Kia] = karssemeijer_radial_projection(label, label_orientation, 20, 80, 15, 24, 10);
figure; 
subplot(1,2,1); imagesc(fi1a); axis image; colormap(jet(256)); colorbar;
subplot(1,2,2); imagesc(fi2a); axis image; colormap(jet(256)); colorbar;
%%    
bg = u_load('C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\bg00001.mat');
label = zeros(512);
label_orientation = zeros(512);
for ii = 0:12:84
    [bar, labeli] = create_sin_bar(4, 4, ii, 512, 512, 0, 256, 256);
    bg = bg + bar;
    label = label | labeli;
    label_orientation(labeli) = ii;
end
%
figure; imagesc(label); axis image;
[fi1b fi2b nib Nib prob_ib ni_plusb Kib] = karssemeijer_radial_projection(label, label_orientation, 20, 80, 15, 24, 10);
figure; 
subplot(1,2,1); imagesc(fi1b); axis image; colormap(jet(256)); colorbar;
subplot(1,2,2); imagesc(fi2b); axis image; colormap(jet(256)); colorbar;
%%
line_map = u_load('C:\isbe\asymmetry_project\data\synthetic_lines\spiral512\results\191905\spiral001_class.mat');
ori_map = u_load('C:\isbe\asymmetry_project\data\synthetic_lines\spiral512\results\191934\spiral001_class.mat');
[fi1 fi2 ni Ni prob_i ni_plus Ki] = karssemeijer_radial_projection(line_map, ori_map, 20, 80, 15, 24, 10);
%%
line_map = imresize(load_uint8('C:\isbe\asymmetry_project\data\line_maps\2004_screening_processed\abnormals\024RCC_class.mat'), 0.5, 'bilinear');
ori_map = imresize(load_uint8('C:\isbe\asymmetry_project\data\orientation_maps\2004_screening_processed\abnormals\024RCC_class.mat'), 0.5, 'bilinear');

px_per_mm = 50 / 9;
R = round(2*px_per_mm);
r_min = round(4*px_per_mm);
r_max = round(16*px_per_mm);
num_angles = [48 24];

[fi1_mb2 fi2_mb2 ni_mb2 Ni2 prob_i_mb2 ni_plus_mb2 Ki_mb2] = karssemeijer_radial_projection(line_map, ori_map, r_min, r_max, R, num_angles, 20);
%%
figure; imagesc(nik); axis image;
figure; imagesc(mik); axis image;
figure; imagesc(Nik); axis image;
%%
mammo = double(imresize(u_load('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\024RCC.mat'), 0.5, 'bilinear'));
mask = imresize(u_load('C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\024RCC_mask.mat'), 0.5);
new_mask = mask;
orig_area = sum(mask(:));
new_area = orig_area;
while new_area > 0.9*orig_area;
    new_mask = imerode(new_mask, strel('disk', 1));
    new_area = sum(new_mask(:));
end
%%
Te = mean(mammo(new_mask));
smooth_mammo = imfilter(mammo, fspecial('average', round(3*px_per_mm)));

figure; 
subplot(1,2,1); imagesc(smooth_mammo); axis image; colormap(gray(256));
subplot(1,2,2); imagesc(smooth_mammo < Te); axis image; colormap(gray(256));
%%
edge_idx = smooth_mammo < Te;
new_mammo = mammo;
new_mammo(edge_idx) = mammo(edge_idx) - smooth_mammo(edge_idx) + Te;
figure; imagesc(new_mammo); axis image; colormap(gray(256));
%%
m_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\*.mat');
for ii = 1:length(m_list);
    movefile(...
        ['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\' m_list(ii).name],...
        ['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\' m_list(ii).name(5:end)]);
    movefile(...
        ['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\' m_list(ii).name(1:end-4) '_mask.mat'],...
        ['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\' m_list(ii).name(5:end-4) '_mask.mat']);
end
%%
px_per_mm = 100 / 9;

% mkdir C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\abnormals\
correct_mammo_dropoff(...
    'C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals',...
    'C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals',...
    'C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\abnormals', round(3*px_per_mm));

% mkdir C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\normals\
correct_mammo_dropoff(...
    'C:\isbe\asymmetry_project\data\mammograms\2004_screening\normals',...
    'C:\isbe\asymmetry_project\data\masks\2004_screening\normals',...
    'C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\normals', round(3*px_per_mm));

%%
px_per_mm = 50 / 9;
mammo = double(imresize(u_load('C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\abnormals\024RCC.mat'), 0.5, 'bilinear'));
mask = imresize(u_load('C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\024RCC_mask.mat'), 0.5);

[line_strength, orientation, scale] = gaussian_2nd_derivative_line(mammo, [1 2 4 8]);
[grad_strength2, grad_orientation] = gaussian_1st_derivative_gradient(mammo, 1*px_per_mm);

discard_idx = ...
    ((abs(mb_mod(orientation - grad_orientation + pi/2, pi)) < pi/6) &...
    (grad_strength > 25)) | ...
    (line_strength > -0.1) | ...
    ~mask;

line_strength(discard_idx) = 0;
figure; imagesc(-line_strength); axis image;

[fi1 fi2 ni Ni prob_i ni_plus Ki] = ...
    karssemeijer_radial_projection(~discard_idx, mod(180*orientation/pi,180), r_min, r_max, R, num_angles, 10);
[fi1a fi2a nia Nia prob_ia ni_plusa Kia] =...
    karssemeijer_radial_projection2(~discard_idx, mod(180*orientation/pi,180), r_min, r_max, R, [48 24], 10);
%%
R = 5;
r_min = 10;
r_max = 100;
num_angles = [48 24];

xy = repmat(1:512, 512, 1);
xx = xy - 256;
yy = xy' - 256; clear xy;

dist_map = sqrt(xx.^2 + yy.^2);
gamma_map = atan(yy ./ xx);
gamma_map(256, 256) = 0;
line_map = ones(512);

S = (dist_map > r_min) & (dist_map < r_max);
Ni = sum(S(:));

np = zeros(100,1);
for ii = 1:100
    ori_map = pi*rand(512) -  pi/2;
    hit_map = abs(mb_mod(ori_map - gamma_map,pi)) < atan(R ./ dist_map) & S;
    np(ii) = sum(hit_map(:));
end
%%
compute_mass_template_map_batch('abnormals',...
    'radii_scales', 5*(6:3:18), 'template_dir', 'template_maps',...
    'mammo_dir', 'mammograms\2004_screening_processed', 'mask_dir', 'masks\2004_screening');
compute_mass_template_map_batch('normals',...
    'radii_scales', 5*(6:3:18), 'template_dir', 'template_maps',...
    'mammo_dir', 'mammograms\2004_screening_processed', 'mask_dir', 'masks\2004_screening');
%%
line_map_o = imresize(load_uint8('C:\isbe\asymmetry_project\data\line_maps\2004_screening\abnormals\024RCC_data.mat'), 0.5, 'bilinear');
ori_map_o = imresize(load_uint8('C:\isbe\asymmetry_project\data\orientation_maps\2004_screening\abnormals\024RCC_data.mat'), 0.5, 'bilinear');
line_map_p = imresize(load_uint8('C:\isbe\asymmetry_project\data\line_maps\2004_screening_processed\abnormals\024RCC_class.mat'), 0.5, 'bilinear');
ori_map_p = imresize(load_uint8('C:\isbe\asymmetry_project\data\orientation_maps\2004_screening_processed\abnormals\024RCC_class.mat'), 0.5, 'bilinear');
relevance_map = imresize(load_uint8('C:\isbe\asymmetry_project\data\relevance_maps\2004_screening\abnormals\024RCC_class.mat'), 0.5, 'bilinear');
%%
px_per_mm = 50 / 9;
R = round(2*px_per_mm);
r_min = round(4*px_per_mm);
r_max = round(16*px_per_mm);
num_angles = [48 24];

[f_i1_oo f_i2_oo n_i_oo N_i_oo p_i_oo] = karssemeijer_radial_projection(line_map_o, ori_map_o, r_min, r_max, R, num_angles, 20);
[f_i1_po f_i2_po n_i_po N_i_po p_i_po] = karssemeijer_radial_projection(line_map_p, ori_map_o, r_min, r_max, R, num_angles, 20);
[f_i1_op f_i2_op n_i_op N_i_op p_i_op] = karssemeijer_radial_projection(line_map_o, ori_map_p, r_min, r_max, R, num_angles, 20);
[f_i1_pp f_i2_pp n_i_pp N_i_pp p_i_pp] = karssemeijer_radial_projection(line_map_p, ori_map_p, r_min, r_max, R, num_angles, 20);

[f_i1_oor f_i2_oor n_i_oor N_i_oor p_i_oor] = karssemeijer_radial_projection(line_map_o .* relevance_map, ori_map_o, r_min, r_max, R, num_angles, 20);
%%
line_map_o = imresize(load_uint8('C:\isbe\asymmetry_project\data\line_maps\2004_screening\abnormals\024RCC_data.mat'), 0.5, 'bilinear');
ori_map_o = imresize(load_uint8('C:\isbe\asymmetry_project\data\orientation_maps\2004_screening\abnormals\024RCC_data.mat'), 0.5, 'bilinear');

[f_i1_a f_i2_a n_i_a N_i_a p_i_a] = karssemeijer_radial_projection(line_map_o, ori_map_o, r_min, r_max, R, num_angles, 20);
[f_i1_b f_i2_b n_i_b N_i_b p_i_b] = karssemeijer_radial_projection(line_map_o > 0.5, ori_map_o, r_min, r_max, R, num_angles, 20);