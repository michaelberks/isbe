roi = load_uint8('C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\mass_roi\024RCC_roi.mat');

linreg = u_load('Z:\data\line_linear_regression_rfs\310934\random_forest.mat');
sargs = u_load('Z:\data\line_linear_regression_rfs\310934\sampling_args.mat');
sampling_args.do_max = sargs.do_max;
sampling_args.rotate = sargs.rotate;
sampling_args.use_nag = sargs.use_nag;
sampling_args.feature_shape = sargs.feature_shape;
sampling_args.feature_type = sargs.feature_type;
sampling_args.num_levels = sargs.num_levels;
sampling_args.win_size = sargs.win_size;

ori_map_l = classify_image(...
    'image_in', roi,... % the mandatory arguments
    'sampling_args',sampling_args,...
    'forest', linreg, ...
    'decomp_type', 'dt',...
    'forest_type', 'linear_regression',...
    'use_probs', 0,...
    'mask', [],...
    'num_trees', [], ...
    'max_size', 128);

load C:\isbe\asymmetry_project\data\misc\ori_maps.mat ori_map_3_4_old
line_map_l = abs(ori_map_3_4_old) > 0.5;
ori_map_l = angle(ori_map_l)/2;
%
r_max = [60 90 120 150 180];
[f_i1l f_i2l] = karssemeijer_radial_projection_multiscale(line_map_l, ori_map_l, 10, r_max, 5, 24, 10, 1);

load C:\isbe\asymmetry_project\data\misc\k_maps.mat f_i*
save C:\isbe\asymmetry_project\data\misc\k_maps.mat f_i*
%%
[goat cow] = karssemeijer_radial_projection_multiscale(abs(ori_map_3_4_old) > 0.5, angle(ori_map_3_4_old), 10, r_max, 5, 24, 10, 1);