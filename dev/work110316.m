bg = u_load('C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\bg00001.mat');


for ii = 18:18:180
    
    [test_im label] = create_sin_bar(3, 4, ii, 512, 512, 0.5, 256, 256);
    
    figure;
    subplot(2,2,1); imagesc(test_im); axis image; colormap(gray(256));
    subplot(2,2,2); imagesc(test_im+bg); axis image; colormap(gray(256));

    [d d orientation_image] = monogenic_phase_cong(test_im, 2, 8, 2, 0.65);
    
    orientation_image = label .* complex(cosd(orientation_image), sind(orientation_image));
    subplot(2,2,3); image(complex2rgb(orientation_image(231+(1:50), 231+(1:50)).^2)); axis image;
    
    [d d orientation_image] = monogenic_phase_cong(test_im+bg, 2, 8, 2, 0.65);
    orientation_image = label .* complex(cosd(orientation_image), sind(orientation_image));
    subplot(2,2,4); image(complex2rgb(orientation_image(231+(1:50), 231+(1:50)).^2)); axis image;
    
    
end
%%
for ii = 1:100
    
    s = load(['C:\isbe\asymmetry_project\data\synthetic_lines\real512\labels\label' zerostr(ii,3) '.mat']);
    mask = s.label == 1;
    save(['C:\isbe\asymmetry_project\data\synthetic_lines\real512\masks\mask' zerostr(ii,3) '.mat'], 'mask');
end
%%
mass_xy = u_load('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\024RCC_meta.mat');

for ii = [60 90 120 150 180]
    f1_map_g2 = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\g2d\2004_screening_processed\abnormals\'...
        '024RML_f1_' zerostr(ii,3) '.mat']);
    f1_map_dt = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\old_rf\2004_screening\abnormals\'...
        '024RML_f1_' zerostr(ii,3) '.mat']);
    
    f2_map_g2 = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\g2d\2004_screening_processed\abnormals\'...
        '024RML_f2_' zerostr(ii,3) '.mat']);
    f2_map_dt = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\old_rf\2004_screening\abnormals\'...
        '024RML_f2_' zerostr(ii,3) '.mat']);

    [r c] = size(f1_map_g2);
    x = mass_xy(:,1) * c;
    y = mass_xy(:,2) * r;
    xc = round(mean(x));
    yc = round(mean(y));

    f1_roi_g2 = sample_window(f1_map_g2, 400, yc, xc,0);
    f1_roi_dt = sample_window(f1_map_dt, 400, yc, xc,0);
    f2_roi_g2 = sample_window(f2_map_g2, 400, yc, xc,0);
    f2_roi_dt = sample_window(f2_map_dt, 400, yc, xc,0);
    
    figure; 
    subplot(2,2,1); imagesc(f1_roi_g2); axis image;
    subplot(2,2,2); imagesc(f1_roi_dt); axis image;
    subplot(2,2,3); imagesc(f2_roi_g2); axis image;
    subplot(2,2,4); imagesc(f2_roi_dt); axis image;
end
%%
mass_xy = u_load('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\024RCC_meta.mat');
ori_map = load_uint8('C:\isbe\asymmetry_project\data\orientation_maps\rf\2004_screening_processed\abnormals\024RML_class.mat');
[r c] = size(ori_map);
x = mass_xy(:,1) * c;
y = mass_xy(:,2) * r;
xc = round(mean(x));
yc = round(mean(y));
ori_map = sample_window(ori_map, 800, yc, xc,0);
pack;
r_max = [120 180 240 300 360];
%%
profile on;
[f_i1 f_i2] = karssemeijer_radial_projection_multiscale(...
    abs(ori_map), angle(ori_map), 10, r_max, 5, 24, 10, 2);
profsave(profile('info'),'C:\isbe\asymmetry_project\data\misc\profiles\k_full2');