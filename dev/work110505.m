roi = load_uint8('C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\mass_roi\002LCC_roi.mat');
%
%Get orientation map of image
rf_ori = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\313384\random_forest.mat');
rf_ori.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';

args.image_in = roi;
args.forest_type = 'orientation';
args.sampling_args.do_max = 0;
args.sampling_args.rotate = 0;
args.sampling_args.use_nag = 0;
args.sampling_args.feature_shape = 'rect';
args.sampling_args.feature_type = 'conj';
args.sampling_args.win_size = 3;
args.sampling_args.num_levels = 3;
args.forest = rf_ori;
ori_map = classify_image(args);
%
%Get line map of image
rf_line = u_load('C:\isbe\asymmetry_project\data\line_detection_rfs\313943\random_forest.mat');
rf_line.tree_root = 'C:\isbe\asymmetry_project\data\line_detection_rfs\';

args.image_in = roi;
args.forest_type = 'isbe';
args.sampling_args.do_max = 0;
args.sampling_args.rotate = 0;
args.sampling_args.use_nag = 0;
args.sampling_args.feature_shape = 'rect';
args.sampling_args.feature_type = 'conj';
args.sampling_args.win_size = 3;
args.sampling_args.num_levels = 3;
args.forest = rf_line;
line_map = classify_image(args);
%
figure; 
subplot(1,2,1); imagesc(line_map); axis image; colormap(gray(256));
subplot(1,2,2); image(complex2rgb(ori_map.^2)); axis image;

%
px_per_mm = 100/9;
sigma_range = [3.2 4.0 5.4 6.2 7.8];
r_max = 3*round(px_per_mm*sigma_range);

mask = false(800);
mask(100:700, 100:700) = 1;
%%
%RF with mag dispersion lines
[f1_rfl f2_rfl mask_out] = karssemeijer_radial_projection_multiscale(...
    line_map, angle(ori_map), 'mask', mask, 'spacing', 4, 'r_max', r_max(3));
%RF with mag dispersion lines
[f1_rfu f2_rfu] = karssemeijer_radial_projection_multiscale(...
    true(800), angle(ori_map), 'mask', mask, 'spacing', 4, 'r_max', r_max(3));

f1_map_l = zeros(size(mask_out));
f1_map_l(mask_out) = f1_rfl; f1_map_l = f1_map_l(1:4:end, 1:4:end);

f2_map_l = zeros(size(mask_out));
f2_map_l(mask_out) = f2_rfl; f2_map_l = f2_map_l(1:4:end, 1:4:end);

f1_map_u = zeros(size(mask_out));
f1_map_u(mask_out) = f1_rfu; f1_map_u = f1_map_u(1:4:end, 1:4:end);

f2_map_u = zeros(size(mask_out));
f2_map_u(mask_out) = f2_rfu; f2_map_u = f2_map_u(1:4:end, 1:4:end);

figure;
subplot(1,2,1); imagesc(f1_map_l); axis image;
subplot(1,2,2); imagesc(f2_map_l); axis image;

figure;
subplot(1,2,1); imagesc(f1_map_u); axis image;
subplot(1,2,2); imagesc(f2_map_u); axis image;
