%%
roi = u_load('C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\mass_roi\002LCC_roi.mat');

% rf_3_4_4 = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\313303\random_forest.mat');
% rf_3_3_4 = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\313300\random_forest.mat');
% rf_3_4_2 = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\313371\random_forest.mat');
rf_3_3_2 = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\313384\random_forest.mat');

% rf_3_4_4.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
% rf_3_3_4.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
% rf_3_4_2.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
rf_3_3_2.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';

args.image_in = roi;
args.forest_type = 'orientation';

args.sampling_args.do_max = 0;
args.sampling_args.rotate = 0;
args.sampling_args.use_nag = 0;
args.sampling_args.feature_shape = 'rect';
args.sampling_args.feature_type = 'conj';
args.sampling_args.win_size = 3;
%
% args.forest = rf_3_4_4;
% args.sampling_args.num_levels = 4;
% ori_map_3_4_4 = classify_image(args);
% %
% args.forest = rf_3_3_4;
% args.sampling_args.num_levels = 3;
% ori_map_3_3_4 = classify_image(args);
%
args.forest = rf_3_4_2;
args.sampling_args.num_levels = 4;
ori_map_3_4_2 = classify_image(args);
%%
args.forest = rf_3_3_2;
args.sampling_args.num_levels = 3;
ori_map_3_3_2 = classify_image(args);

figure; image(complex2rgb(ori_map_3_4_4.^2)); axis image;
figure; image(complex2rgb(ori_map_3_3_4.^2)); axis image;
figure; image(complex2rgb(ori_map_3_4_2.^2)); axis image;
figure; image(complex2rgb(ori_map_3_3_2.^2)); axis image;
%%
bg_idxs = [];
widths = [];
contrasts = [];
orientations = [];
squashes = [];
radii = [];

for jj = 1:20
    for kk = 1:10
        load(['Z:\data\line_orientation_rfs\313303\line_parameters\' zerostr(kk,2) '\parameters' zerostr(jj,3) '.mat']);
        for ii = 1:length(parameters)
            bg_idxs(end+1,1)=parameters(ii).bg_idx; %#ok
            widths(end+1,1)=parameters(ii).width; %#ok
            contrasts(end+1,1)=parameters(ii).contrast; %#ok
            orientations(end+1,1)=parameters(ii).orientation; %#ok
            squashes(end+1,1)=parameters(ii).squash; %#ok
            radii(end+1,1)=parameters(ii).radius; %#ok
        end
    end
end