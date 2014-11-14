load C:\isbe\asymmetry_project\data\misc\ori_maps.mat ori_map_1_5
ori_map = angle(ori_map_1_5);
line_map = abs(ori_map_1_5) > 0.25;

r_max = [60 90 120 150 180];
%
tic;
[f_i1a f_i2a] = karssemeijer_radial_projection_pix2(line_map, ori_map, 10, r_max, 5, 24, 50, 16);
toc;

tic;
for scales = 1:5
    [f_i1b f_i2b] = karssemeijer_radial_projection_pix(line_map, ori_map, 10, r_max(scales), 5, 24, 50, 16);
    
    diff1 = f_i1b(1:16:end,1:16:end) - f_i1a(1:16:end,1:16:end,scales);
    diff2 = f_i2b(1:16:end,1:16:end) - f_i2a(1:16:end,1:16:end,scales);
    
    display(['Max diff1 = ' num2str(max(diff1(:))) ', mean diff1 = ' num2str(mean(diff1(:)))]);
    display(['Max diff2 = ' num2str(max(diff2(:))) ', mean diff2 = ' num2str(mean(diff2(:)))]);
end
toc;
%%
profile on;
[f_i1b f_i2b] = karssemeijer_radial_projection_pix2(line_map, ori_map, 10, r_max, 5, 24, 50, 4);
profile viewer;
%%
load C:\isbe\asymmetry_project\data\misc\ori_maps.mat ori_map_1_5
ori_map_1_5 = imresize(ori_map_1_5, 0.5);
ori_map = angle(ori_map_1_5);
line_map = abs(ori_map_1_5) > 0.25;
r_max = [60 90 120 150 180];
tic;
[f_i1 f_i2] = karssemeijer_radial_projection_multiscale(line_map, ori_map, 10, r_max, 5, 24, 10, 1);    
toc;
%%
ori_map_3_4 = imresize(ori_map_3_4, 0.5);
ori_map = angle(ori_map_3_4);
line_map = abs(ori_map_3_4) > 0.5;
r_max = [60 90 120 150 180];
tic;
[f_i1a f_i2a] = karssemeijer_radial_projection_multiscale(line_map, ori_map, 10, r_max, 5, 24, 10, 1);    
toc;
%%
ori_map_3_4_old = imresize(ori_map_3_4_old, 0.5);
ori_map = angle(ori_map_3_4);
line_map = abs(ori_map_3_4);
tic;
[f_i1b f_i2b] = karssemeijer_radial_projection_multiscale(line_map, ori_map, 10, r_max, 5, 24, 10, 1);    
toc;
%%
ori_map_g2d = load_uint8('C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\mass_roi\024RCC_roi.mat');
ori_map_g2d = pi*imresize(ori_map_g2d, 0.5)/180;
line_map_g2d = load_uint8('C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\mass_roi\024RCC_roi.mat');
line_map_g2d = imresize(line_map_g2d, 0.5);
tic;
[f_i1g f_i2g] = karssemeijer_radial_projection_multiscale(line_map_g2d, ori_map_g2d, 10, r_max, 5, 24, 10, 1);    
toc;