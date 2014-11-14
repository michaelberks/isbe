o2 = load_uint8('C:\isbe\asymmetry_project\data\orientation_maps\rf_thin\2004_screening_processed\abnormals\002LCC_class.mat');
m2 = u_load('C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\002LCC_mask.mat');

%figure; imagesc(m2 & (angle(o2) < -pi/2 | angle(o2) > pi/2)); axis image
cm = m2 & (angle(o2) < -pi/2 | angle(o2) > pi/2);
clear o2 m2;

rf = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\313384\random_forest.mat');
rf.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
args.image_in = u_load('C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\abnormals\002LCC.mat');
args.forest_type = 'orientation';
args.sampling_args.do_max = 0;
args.sampling_args.rotate = 0;
args.sampling_args.use_nag = 0;
args.sampling_args.feature_shape = 'rect';
args.sampling_args.feature_type = 'conj';
args.sampling_args.win_size = 3;
args.sampling_args.num_levels = 3;
args.max_size = 5e3;
args.mask = cm;
args.forest = rf;
clear cm;
pack;
classify_image(args);