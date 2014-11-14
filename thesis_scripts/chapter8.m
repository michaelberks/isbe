%--------------------------------------------------------------------------
% Chapter 8: Combining regions (this may actually be chapter 9, but hey,
% who's counting?)
%--------------------------------------------------------------------------
mkdir C:\isbe\thesis\figures\9
%%
%--------------------------------------------------------------------------
% Show plot of first 2 combined mass model modes, with contours of Gaussian
% plot
%--------------------------------------------------------------------------
load C:\isbe\dev\mass_model\models\model_w500_50K.mat
[xx yy] = meshgrid(...
    linspace(-3*sqrt(mass_model.L_com(1)), 3*sqrt(mass_model.L_com(1)), 50),...
    linspace(-3*sqrt(mass_model.L_com(2)), 3*sqrt(mass_model.L_com(2)), 50));
gxy = diag(exp(-.5*[xx(:) yy(:)]*[1/mass_model.L_com(1) 0; 0 1/mass_model.L_com(2)]*[xx(:) yy(:)]'));
figure; contour(xx,yy,reshape(gxy, 50,50)); hold on;
plot(mass_model.B_com(1,:), mass_model.B_com(2,:), 'rx');

text(-2*sqrt(mass_model.L_com(1)),0, '-2 \sigma_1');
text(2*sqrt(mass_model.L_com(1)),0, '+2 \sigma_1');
text(0,-2*sqrt(mass_model.L_com(2)), '-2 \sigma_2');
text(0,2*sqrt(mass_model.L_com(2)), '+2 \sigma_2');
xlabel('Model parameters along 1st mode');
ylabel('Model parameters along 2nd mode');
%%
%--------------------------------------------------------------------------
% Conditionally sampling a mass
%--------------------------------------------------------------------------
my_files = dir('C:\isbe\dev\masses1024x1024\*.mat');
for ii = 173
    try
        figure;
        mass = u_load(['C:\isbe\dev\masses1024x1024\', my_files(ii).name]);
        mass_sub = mass.mass_ROI-mass.background_ROI;
        mass_sub(mass_sub < 0) = 0;
        
        mass_bg = mass.mass_ROI(201:824, 201:824);
        mass_bgr = rot90(mass_sub, 1) + mass.background_ROI;
        mass_bgr = mass_bgr(201:824, 201:824);
         
        %subplot(2,2,1); imagesc(mass.mass_ROI); axis image; colormap(gray(256)); hold on;
        %plot(mass.mass_outline(:,1), mass.mass_outline(:,2), 'r--');
        %subplot(2,2,2); imagesc(mass.background_ROI); axis image; colormap(gray(256)); hold on;
        %plot(mass.mass_outline(:,1), mass.mass_outline(:,2), 'r--');
        %subplot(2,2,3); imagesc(mass_sub); axis image; colormap(gray(256)); hold on;
        %plot(mass.mass_outline(:,1), mass.mass_outline(:,2), 'r--');
        %subplot(2,2,4); imagesc(rot90(mass_sub)+mass.background_ROI); axis image; colormap(gray(256)); hold on;
        %plot(mass.mass_outline(:,2), 1024-mass.mass_outline(:,1), 'r--');
        subplot(1,2,1); imagesc(mass_bg); axis image; colormap(gray(256)); hold on;
        plot(mass.mass_outline(:,1) - 200, mass.mass_outline(:,2) - 200, 'r--');
        subplot(1,2,2); imagesc(mass_bgr); axis image; colormap(gray(256)); hold on;
        plot(mass.mass_outline(:,2) - 200, 1024-mass.mass_outline(:,1) - 200, 'r--');
        clear mass
    end
end
%%
% 135RML
anno = u_load('C:\isbe\dev\annotations\an04_135RML.mat');
figure; imagesc(anno.mass_ROI); axis image; colormap(gray(256)); hold on;
plot(anno.mass_outline(:,1), anno.mass_outline(:,2), 'r--');
%%
for ii = 1:5
    plot(anno.mass_spicules(ii).outline(:,1), anno.mass_spicules(ii).outline(:,2));
end
%%
mass = u_load('C:\isbe\dev\masses1024x1024\m04_135RML.mat');
mass_sub = mass.mass_ROI-mass.background_ROI;
mass_sub(mass_sub < 0) = 0;

mass_bg = mass.mass_ROI;
mass_bgr = rot90(mass_sub, 1) + mass.background_ROI;

%%
f1 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [100 100 512 512],...
    'PaperPositionMode','auto');
axes('units', 'pixels', 'position', [0 0 512 512]);

imagesc(mass.mass_ROI); colormap(gray(256)); axis image; hold on; caxis([140 max(mass_bgr(:))]);
plot(mass.mass_outline(:,1), mass.mass_outline(:,2), 'r--', 'LineWidth', 1.5);
plot(x1, y1, 'r--', 'LineWidth', 1.5);
plot(x2, y2, 'r--', 'LineWidth', 1.5);
plot(x3, y3, 'r--', 'LineWidth', 1.5);
plot(x4, y4, 'r--', 'LineWidth', 1.5);
plot(x5, y5, 'r--', 'LineWidth', 1.5);
plot(x6, y6, 'r--', 'LineWidth', 1.5);
plot(x7, y7, 'r--', 'LineWidth', 1.5);
plot(x8, y8, 'r--', 'LineWidth', 1.5);

set(gca, 'xticklabel', [], 'yticklabel', []);

print('-dtiff', '-noui', '-painters', f1, '-r300', 'C:\isbe\thesis\figures\9\mass_unrotated.tif');
%%
f2 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [100 100 512 512],...
    'PaperPositionMode','auto');
axes('units', 'pixels', 'position', [0 0 512 512]);

imagesc(mass_bgr); colormap(gray(256)); axis image; hold on; caxis([140 max(mass_bgr(:))]);
plot(mass.mass_outline(:,2), 1024-mass.mass_outline(:,1), 'r--', 'LineWidth', 1.5);
plot(x1, y1, 'r--', 'LineWidth', 1.5);
plot(x2, y2, 'r--', 'LineWidth', 1.5);
plot(x3, y3, 'r--', 'LineWidth', 1.5);
plot(x4, y4, 'r--', 'LineWidth', 1.5);
plot(x5, y5, 'r--', 'LineWidth', 1.5);
plot(x6, y6, 'r--', 'LineWidth', 1.5);
plot(x7, y7, 'r--', 'LineWidth', 1.5);
plot(x8, y8, 'r--', 'LineWidth', 1.5);

set(gca, 'xticklabel', [], 'yticklabel', []);

print('-dtiff', '-noui', '-painters', f2, '-r300', 'C:\isbe\thesis\figures\9\mass_rotated.tif');
%%
save C:\isbe\thesis\figures\9\rotate_mass_data mass x* y*
%%
load C:\isbe\dev\mass_model\models\model_w500_50K.mat
mass_idx = 23;
sample_args.TargetRegion = zeros(1024);
sample_args.TargetCentre = [512 512];
sample_args.ConditionMode = mass_idx;
sample_args.Rotation = mass_model.rotations(:,:,mass_idx);
sample_args.Origin = mass_model.origins(mass_idx);
sample_args.MassModel = mass_model;
sample_args.NumModes = length(mass_model.L_com);
%
for jj = 0
    sample_args.NumModes = length(mass_model.L_com);
    [mass] = sample_new_mass_in_region(sample_args);
    
    figure('Name', ['Number of consitioning modes = ', num2str(jj)]); 
    subplot(2,3,1); imagesc(mass.subtract_ROI); axis image; colormap(gray(256));
    sample_args.NumModes = jj;
    for ii = 2:6
        [mass] = sample_new_mass_in_region(sample_args);
        subplot(2,3,ii); imagesc(mass.subtract_ROI); axis image; colormap(gray(256));
    end
end
%%
load C:\isbe\dev\mass_model\models\model_w500_50K.mat
mass_idx = 23;
sample_args.TargetRegion = zeros(1024);
sample_args.TargetCentre = [512 512];
sample_args.ConditionMode = mass_idx;
sample_args.Rotation = mass_model.rotations(:,:,mass_idx);
sample_args.Origin = mass_model.origins(mass_idx);
sample_args.MassModel = mass_model;

sample_args.NumModes = length(mass_model.L_com);
[mass] = sample_new_mass_in_region(sample_args);
mass_region = mass.subtract_ROI(65:end-64, 65:end-64);
write_im_from_colormap(mass_region, 'C:\isbe\thesis\figures\9\mass23_f.bmp', gray(256), [0 50]);
for m = [1 2 3 4 5 10 0]
    sample_args.NumModes = m;
    for ii = 1:4
        mass_name = ['C:\isbe\thesis\figures\9\mass23_f', zerostr(m,2), 's', num2str(ii), '.bmp'];
        [mass] = sample_new_mass_in_region(sample_args);
        mass_region = mass.subtract_ROI(65:end-64, 65:end-64); 
        write_im_from_colormap(mass_region, mass_name, gray(256), [0 50]);
    end
end

%%
sample_args.NumModes = length(mass_model.L_com);
[mass] = sample_new_mass_in_region(sample_args);

figure('Name', ['Number of consitioning modes = ', num2str(0)]); 
subplot(2,3,1); imagesc(mass.subtract_ROI); axis image; colormap(gray(256));
sample_args.ConditionMode = [];
for ii = 2:6
    [mass] = sample_new_mass_in_region(sample_args);
    subplot(2,3,ii); imagesc(mass.subtract_ROI); axis image; colormap(gray(256));
end
%%
for ii = 1:6;
    subplot(2,3,ii); colorbar off; caxis([0 45]);
end
%%
load C:\isbe\dev\mass_model\models\model_w500_50K.mat
%
P_shape = mass_model.P_shape;
L_shape = mass_model.L_shape;
mean_shape = mass_model.mean_shape;
mean_shape_pl = mass_model.mean_shape_pl;

P_tex = mass_model.P_tex;
L_tex = mass_model.L_tex;
mean_tex = mass_model.mean_tex;

P_scale = mass_model.P_scale;
mean_scale = mass_model.mean_scale;

P_com = mass_model.P_com;
L_com = mass_model.L_com;
mean_com = mass_model.mean_com;

W_shape     = mass_model.W_shape;
W_tex       = mass_model.W_tex;
W_scale     = mass_model.W_scale;

k_shape     = length(L_shape);
k_tex       = length(L_tex);

Q_shape = P_com(1:k_shape,:); 
Q_tex = P_com(k_shape+1:k_shape + k_tex,:);
Q_scale = P_com(end, :);

scaled_mean_shape = mass_model.scaled_mean_shape;

%Define source points for TPS - as row vectors
s_x = scaled_mean_shape(1:end/2);% + mean_centre(1);
s_y = scaled_mean_shape(end/2+1:end);% + mean_centre(2);

%Define points to be interpolated by TPS - as row vectors
i_x = mean_shape_pl(:,1)';
i_y = mean_shape_pl(:,2)';

height = 500;
width = 500;
off_r = height / 2;
off_c = width / 2;
%
for ii = 1:10%length(L_com)
    c_vec = zeros(length(L_com),1);
    c_vec(ii) = 2*sqrt(L_com(ii));
    
    c_shape = mean_shape + (P_shape*inv(W_shape)*Q_shape*c_vec)';
    c_tex = mean_tex + (P_tex*inv(W_tex)*Q_tex*c_vec)';
    c_scale = mean_scale + (P_scale*Q_scale*c_vec)' / W_scale;    
    c_shape = c_shape*c_scale;
    
    display(['Scale = ', num2str(c_scale)]);
    
    c_shape(end/2+1:end) = c_shape(end/2+1:end) + off_r;
    c_shape(1:end/2) = c_shape(1:end/2) + off_c;

    c_bw = roipoly(height, width, c_shape(1:end/2),...
        c_shape(end/2+1:end));

    rp = regionprops(bwlabel(c_bw, 4), 'PixelList');
    c_shape_pl = rp.PixelList; clear rp;

    %Define displacement to target points
    z_x = c_shape(1:end/2);
    z_y = c_shape(end/2+1:end);

    %Compute displacement of interpolated points        
    %         f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
    %         f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
    T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
            'transform', 'spline');
    [pts] = geom_transformpoints([i_x; i_y], T);
    f_x = pts(1,:);
    f_y = pts(2,:);

    %display('completed warping');
    c_shape_tex = griddata(f_x, f_y, c_tex,...
        c_shape_pl(:,1), c_shape_pl(:,2));

    clear new_tex;

    c_ROIp = zeros(height, width);
    c_ROIp(sub2ind([height width], c_shape_pl(:,2), c_shape_pl(:,1))) = c_shape_tex;
    
    c_vec(ii) = -2*sqrt(L_com(ii));
    
    c_shape = mean_shape + (P_shape*inv(W_shape)*Q_shape*c_vec)';
    c_tex = mean_tex + (P_tex*inv(W_tex)*Q_tex*c_vec)';
    c_scale = mean_scale + (P_scale*Q_scale*c_vec)' / W_scale;    
    c_shape = c_shape*c_scale;
    
    display(['Scale = ', num2str(c_scale)]);

    c_shape(end/2+1:end) = c_shape(end/2+1:end) + off_r;
    c_shape(1:end/2) = c_shape(1:end/2) + off_c;

    c_bw = roipoly(height, width, c_shape(1:end/2),...
        c_shape(end/2+1:end));

    rp = regionprops(bwlabel(c_bw, 4), 'PixelList');
    c_shape_pl = rp.PixelList; clear rp;

    %Define displacement to target points
    z_x = c_shape(1:end/2);
    z_y = c_shape(end/2+1:end);

    %Compute displacement of interpolated points        
    %         f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
    %         f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
    T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
            'transform', 'spline');
    [pts] = geom_transformpoints([i_x; i_y], T);
    f_x = pts(1,:);
    f_y = pts(2,:);

    %display('completed warping');
    c_shape_tex = griddata(f_x, f_y, c_tex,...
        c_shape_pl(:,1), c_shape_pl(:,2));

    clear new_tex;

    c_ROIn = zeros(height, width);
    c_ROIn(sub2ind([height width], c_shape_pl(:,2), c_shape_pl(:,1))) = c_shape_tex;
    
    figure; 
    subplot(1,2,1); imagesc(c_ROIn); axis image; colormap(gray(256)); caxis([0 50]);
    subplot(1,2,2); imagesc(c_ROIp); axis image; colormap(gray(256)); caxis([0 50]);
end
%%
for ii = 1:length(L_com)
    c_vec = zeros(length(L_com),1);
    c_vec(1:ii) = 2*sqrt(L_com(1:ii));
    
    c_shape = mean_shape + (P_shape*inv(W_shape)*Q_shape*c_vec)';
    c_tex = mean_tex + (P_tex*inv(W_tex)*Q_tex*c_vec)';
    c_scale = mean_scale + (P_scale*Q_scale*c_vec)' / W_scale;    
    c_shape = c_shape*c_scale;
    
    display(['Scale = ', num2str(c_scale)]);
    
    c_shape(end/2+1:end) = c_shape(end/2+1:end) + off_r;
    c_shape(1:end/2) = c_shape(1:end/2) + off_c;

    c_bw = roipoly(height, width, c_shape(1:end/2),...
        c_shape(end/2+1:end));

    rp = regionprops(bwlabel(c_bw, 4), 'PixelList');
    c_shape_pl = rp.PixelList; clear rp;

    %Define displacement to target points
    z_x = c_shape(1:end/2);
    z_y = c_shape(end/2+1:end);

    %Compute displacement of interpolated points        
    %         f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
    %         f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
    T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
            'transform', 'spline');
    [pts] = geom_transformpoints([i_x; i_y], T);
    f_x = pts(1,:);
    f_y = pts(2,:);

    %display('completed warping');
    c_shape_tex = griddata(f_x, f_y, c_tex,...
        c_shape_pl(:,1), c_shape_pl(:,2));

    clear new_tex;

    c_ROIp = zeros(height, width);
    c_ROIp(sub2ind([height width], c_shape_pl(:,2), c_shape_pl(:,1))) = c_shape_tex;
    
    c_vec(1:ii) = -2*sqrt(L_com(1:ii));
    
    c_shape = mean_shape + (P_shape*inv(W_shape)*Q_shape*c_vec)';
    c_tex = mean_tex + (P_tex*inv(W_tex)*Q_tex*c_vec)';
    c_scale = mean_scale + (P_scale*Q_scale*c_vec)' / W_scale;    
    c_shape = c_shape*c_scale;
    
    display(['Scale = ', num2str(c_scale)]);

    c_shape(end/2+1:end) = c_shape(end/2+1:end) + off_r;
    c_shape(1:end/2) = c_shape(1:end/2) + off_c;

    c_bw = roipoly(height, width, c_shape(1:end/2),...
        c_shape(end/2+1:end));

    rp = regionprops(bwlabel(c_bw, 4), 'PixelList');
    c_shape_pl = rp.PixelList; clear rp;

    %Define displacement to target points
    z_x = c_shape(1:end/2);
    z_y = c_shape(end/2+1:end);

    %Compute displacement of interpolated points        
    %         f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
    %         f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
    T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
            'transform', 'spline');
    [pts] = geom_transformpoints([i_x; i_y], T);
    f_x = pts(1,:);
    f_y = pts(2,:);

    %display('completed warping');
    c_shape_tex = griddata(f_x, f_y, c_tex,...
        c_shape_pl(:,1), c_shape_pl(:,2));

    clear new_tex;

    c_ROIn = zeros(height, width);
    c_ROIn(sub2ind([height width], c_shape_pl(:,2), c_shape_pl(:,1))) = c_shape_tex;
    
    figure; 
    subplot(1,2,1); imagesc(c_ROIn); axis image; colormap(gray(256)); caxis([0 50]);
    subplot(1,2,2); imagesc(c_ROIp); axis image; colormap(gray(256)); caxis([0 50]);
end
%%
%--------------------------------------------------------------------------
sample_args.NumModes = length(mass_model.L_com);
[mass] = sample_new_mass_in_region(sample_args);

figure('Name', ['Number of consitioning modes = ', num2str(0)]); 
subplot(2,3,1); imagesc(mass.subtract_ROI); axis image; colormap(gray(256));
sample_args.ConditionMode = [];
for ii = 2:6
    [mass] = sample_new_mass_in_region(sample_args);
    subplot(2,3,ii); imagesc(mass.subtract_ROI); axis image; colormap(gray(256));
end
%%
for ii = 1:6;
    subplot(2,3,ii); colorbar off; caxis([0 45]);
end
%%
%--------------------------------------------------------------------------
load C:\isbe\dev\mass_model\models\model_w500_50K.mat

D = length(mass_model.L_com);
%%
scale_change = zeros(D,1);
Q_scale = mass_model.P_com(end, :);
for ii = 1:D
    b_com = zeros(D,1);
    b_com(ii) = 2*sqrt(mass_model.L_com(ii));
    
    b_scale = Q_scale * b_com / mass_model.W_scale;
    
    scale_change(ii) = b_scale / sqrt(mass_model.L_scale);
end
%%
shape_change = zeros(D,5);
for mode = 1:5
    Q_shape = mass_model.P_com(mode, :);
    for ii = 1:D
        b_com = zeros(D,1);
        b_com(ii) = 2*mass_model.L_com(ii);

        b_shape = Q_shape * b_com / mass_model.W_shape;

        shape_change(ii,mode) = b_shape / sqrt(mass_model.L_shape(mode));
    end
end
%%
tex_change = zeros(D,5);
for mode = 1:5
    Q_tex = mass_model.P_com(length(mass_model.L_shape)+mode, :);
    for ii = 1:D
        b_com = zeros(D,1);
        b_com(ii) = 2*mass_model.L_com(ii);

        b_tex = Q_tex * b_com / mass_model.W_tex;

        tex_change(ii,mode) = b_tex / sqrt(mass_model.L_tex(mode));
    end
end
%%
mr_list = dir('D:\isbe\dev\background\modified_regions\*.mat');
for ii = 61:80
    load(['D:\isbe\dev\background\modified_regions\', mr_list(ii).name]);
    figure;
    subplot(1,2,1); imagesc(target_region); axis image; colormap(gray(256));
    subplot(1,2,2); imagesc(modified_region); axis image; colormap(gray(256));
end
%%
norm_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');
load C:\isbe\dev\mass_model\models\model_w500_50K.mat
mass_dir = 'C:\isbe\dev\masses1024x1024';
load C:\isbe\dev\files\bg_files.mat
%%
combo = [6 51 6 15 16 31 40 51 8 15 51 6 16 16; 2 2 3 3 3 3 3 3 4 4 4 5 4 5]';
for ii = 13:size(combo,1)
    mm = combo(ii,1);
    nn = combo(ii,2);
    target_region = double(imread(['C:\isbe\dev\background\images\normal1024\', norm_list(nn).name]));
    [mass] = ...
        mb_synthesise_mass(target_region, [512 512], mass_model,...
        mass_dir, model_id.mass_files, bg_good(mm), []);
    figure;
    subplot(1,2,1); imagesc(target_region); axis image; colormap(gray(256));
    subplot(1,2,2); imagesc(mass.mass_ROI); axis image; colormap(gray(256));
    
end
%%
norm_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');
load C:\isbe\dev\mass_model\models\model_w500_50K.mat
mass_dir = 'C:\isbe\dev\masses1024x1024';
for nn = 1:length(norm_list)
    for mm = [14 26 27 77]
        try
            target_region = double(imread(['C:\isbe\dev\background\images\normal1024\', norm_list(nn).name]));
            [mass] = ...
                mb_synthesise_mass(target_region, [512 512], mass_model,...
                mass_dir, model_id.mass_files, mm, []);
            mass_name = ['C:\isbe\thesis\figures\9\norm', zerostr(nn,3), 'mass', zerostr(mm,3), '.bmp'];
            write_im_from_colormap(mass.mass_ROI, mass_name, colormap(gray(256)));
        catch
            display('Skipped');
        end
    end
%     figure;
%     subplot(1,2,1); imagesc(target_region); axis image; colormap(gray(256));
%     subplot(1,2,2); imagesc(mass.mass_ROI); axis image; colormap(gray(256));
    
end
%%
load C:\isbe\dev\files\bg_files.mat

for ii = 1:length(bg_files)
    mass = u_load(['C:\isbe\dev\masses1024x1024\', bg_files(ii).name]);
    figure; imagesc(mass.background_ROI); axis image; colormap(gray(256));
    clear mass;
end
%%
norm003 = double(imread(['C:\isbe\dev\background\images\normal1024\', norm_list(3).name]));
norm005 = double(imread(['C:\isbe\dev\background\images\normal1024\', norm_list(5).name]));
norm028 = double(imread(['C:\isbe\dev\background\images\normal1024\', norm_list(28).name]));

write_im_from_colormap(norm003,'C:\isbe\thesis\figures\9\norm003.bmp', colormap(gray(256)));
write_im_from_colormap(norm005,'C:\isbe\thesis\figures\9\norm005.bmp', colormap(gray(256)));
write_im_from_colormap(norm028,'C:\isbe\thesis\figures\9\norm028.bmp', colormap(gray(256)));

norm003m = imread('C:\isbe\thesis\figures\9\norm003mass026.bmp');
norm005m = imread('C:\isbe\thesis\figures\9\norm005mass014.bmp');
norm028m = imread('C:\isbe\thesis\figures\9\norm028mass077.bmp');

figure; 
subplot(1,2,1); imagesc(norm003); axis image; colormap(gray(256)); %caxis([min(norm003m(:)) max(norm003m(:))]);
subplot(1,2,2); imagesc(norm003m); axis image; colormap(gray(256)); %caxis([min(norm003m(:)) max(norm003m(:))]);

figure; 
subplot(1,2,1); imagesc(norm005); axis image; colormap(gray(256)); %caxis([min(norm005m(:)) max(norm005m(:))]);
subplot(1,2,2); imagesc(norm005m); axis image; colormap(gray(256)); %caxis([min(norm005m(:)) max(norm005m(:))]);

figure; 
subplot(1,2,1); imagesc(norm028); axis image; colormap(gray(256)); %caxis([min(norm015m(:)) max(norm015m(:))]);
subplot(1,2,2); imagesc(norm028m); axis image; colormap(gray(256)); %caxis([min(norm015m(:)) max(norm015m(:))]);
