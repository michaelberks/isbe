clear;

load C:\isbe\dev\masses\m_files.mat
load(['C:\isbe\dev\masses\', m_files(31).name]);
load C:\isbe\dev\temp\model031
mass_model = combine_model(mass_model);

mean_shape  = mass_model.mean_shape;
P_shape     = mass_model.P_shape;
mean_scale  = mass_model.mean_scale;
P_scale     = mass_model.P_scale;
mean_tex    = mass_model.mean_tex;
P_tex       = mass_model.P_tex;
mean_c      = mass_model.mean_c;
P_c         = mass_model.P_c;
W_shape     = mass_model.W_shape;
W_tex       = mass_model.W_tex;
W_scale     = mass_model.W_scale;

mean_shape_pl = mass_model.mean_shape_pl;
size_shape_vec = length(mean_shape) / 2;
k_shape     = size(P_shape, 2);
k_tex       = size(P_tex, 2);

scale_factor = mass_model.scale_factor;
scaled_mean = mean_shape / scale_factor;

mean_off_r  = mass_model.mean_off_r;
mean_off_c  = mass_model.mean_off_c;

idx = round(linspace(1, length(mass.mass_outline(:,1)), size_shape_vec+1));
idx = idx(1:size_shape_vec); %ensures first point is not equal to the last point!
x_shape = mass.mass_outline(idx,:);
[dd Z t] = procrustes(reshape(scaled_mean, size_shape_vec, 2), x_shape);

x_pro = [Z(:,1)', Z(:,2)'];
x_scale = t.b;

%Define points to be interpolated by TPS - as row vectors
i_x = mean_shape_pl(:,1)';
i_y = mean_shape_pl(:,2)';

%Define source points for TPS - as row vectors
s_x = mean_shape(1:size_shape_vec) + mean_off_c;
s_y = mean_shape(size_shape_vec+1:end) + mean_off_r;

tps_L_inv = tps_weights(s_x, s_y);

%Define displacement to target points
z_x = x_shape(:,1)' + mass.mass_centroid(1);
z_y = x_shape(:,2)' + mass.mass_centroid(2);

%Compute displacement of interpolated points
f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);

figure;
imagesc(mass.subtract_ROI); colormap(gray(256)); axis image; hold on;
plot(f_x, f_y, 'r.');

%Get texture vector
x_tex = interp2(mass.subtract_ROI, f_x, f_y);

x_tex(isnan(x_tex)) = 0;

b_shape = P_shape' * (x_pro-scaled_mean)';
b_tex = P_tex' * (x_tex-mean_tex)';
b_scale = P_scale' * (x_scale-mean_scale)';
x_c = [W_shape*b_shape; W_tex*b_tex; W_scale*b_scale]';
b_c = P_c' * (x_c - mean_c)';

Q_shape = P_c(1:k_shape,:); 
Q_tex = P_c(k_shape+1:k_shape + k_tex,:);
Q_scale = P_c(end, :);

new_shape = scaled_mean + (P_shape*Q_shape*b_c)' / W_shape;
new_tex = mean_tex + (P_tex*Q_tex*b_c)' / W_tex;
new_scale = mean_scale + (P_scale*Q_scale*b_c)' / W_scale;

%new_shape = reshape(new_shape, size_shape_vec, 2);
%new_shape = new_shape*inv(t.T) / new_scale;

%new_shape(:,1) = new_shape(:,1) + mass.mass_centroid(1);
%new_shape(:,2) = new_shape(:,2) + mass.mass_centroid(2);

b_shape_c = Q_shape*b_c / W_shape;
b_tex_c = Q_tex*b_c / W_tex;
b_scale_c = Q_scale*b_c / W_scale;
x_shape_new = P_shape*b_shape + scaled_mean';
x_scale_new = P_scale*b_scale + mean_scale';
x_tex_new = P_tex*b_tex + mean_tex';

figure('WindowStyle', 'docked'); plot(x_pro(1:500), x_pro(501:1000), 'b');
hold on;
plot(x_shape_new(1:500), x_shape_new(501:1000), 'r:');
plot(new_shape(1:500), new_shape(501:1000), 'g:');

mean_row = ceil(max(mean_shape(size_shape_vec+1:2*size_shape_vec)))...
    + mean_off_r + 100;
mean_col = ceil(max(mean_shape(1:size_shape_vec)))...
    + mean_off_c + 100;

figure('WindowStyle', 'docked');
temp1 = zeros(mean_row, mean_col);
temp1(sub2ind(size(temp1), mean_shape_pl(:,1), mean_shape_pl(:,2))) = x_tex;
temp1(temp1 < 0) = 0;
subplot(2,2,1);
imagesc(temp1); axis image; colormap(gray(256)); hold on;
clear temp1;

temp2 = zeros(mean_row, mean_col);
temp2(sub2ind(size(temp2), mean_shape_pl(:,1), mean_shape_pl(:,2))) = x_tex_new;
temp2(temp2 < 0) = 0;
subplot(2,2,2);
imagesc(temp2); axis image; colormap(gray(256)); hold on;
clear temp1;

temp3 = zeros(mean_row, mean_col);
temp3(sub2ind(size(temp3), mean_shape_pl(:,1), mean_shape_pl(:,2))) = new_tex;
temp3(temp3 < 0) = 0;
subplot(2,2,3);
imagesc(temp3); axis image; colormap(gray(256)); hold on;
clear temp1;
%%

N = length(m_files);

for ii = 58:59
    file_out = ['C:\isbe\dev\masses\left_out_models\model', zerostr(ii, 3)];
    one_out_files = m_files;
    one_out_files(ii) = [];

    mass_model = generate_mass_AM(...
        one_out_files, file_out, 500, 'C:\isbe\dev\masses\', 20000);
end
%%
idx2 = setdiff(1:179, idx);
d_files = m_files(idx2);
calculate_weights('mass_model', d_files, 'masses\', 'weights\RMS_random1');