clc; clear;

% set sampling arguments
sampling_args.levels = 1:4;
sampling_args.feature_shape = 'rect';
sampling_args.feature_type = 'conj';
sampling_args.do_max = 0;
sampling_args.rotate = 0;
sampling_args.win_size = 1;
size_sample_vector = length(sampling_args.levels)*compute_dt_feature_size(sampling_args);

% randomly select background - will load a matrix 'bg'
args.bg_dir = 'U:\projects\mammography\data\synthetic_backgrounds/real512/train/';
bg_list = dir([args.bg_dir, '*.mat']);
args.num_bgs = length(bg_list);
bg_idx = 241;
bg_filename = [args.bg_dir bg_list(bg_idx).name];
bg = u_load(bg_filename);

% generate the image
args.width_range = [8 8];
args.orientation_range = [45 45];
args.contrast_range = [16 16];
args.squash_range = 0.5*[1 1];
args.decay_rate = 4;
args.line_type = 'sin';


%Compute dual-tree transform of image
[image_in,label,label_centre,label_orientation,params] = ...
	generate_line_image(bg,args);
dt = compute_dual_tree(image_in,4,0);
dt_coeffs_pos = sample_dt_data(dt,256,256,sampling_args)

dt = compute_dual_tree(255-image_in,4,0);
dt_coeffs_neg = sample_dt_data(dt,256,256,sampling_args)

args.polarity = 'neg';
[image_in,label,label_centre,label_orientation,params] = ...
	generate_line_image(bg,args);
dt = compute_dual_tree(image_in,4,0);
dt_coeffs_neg2 = sample_dt_data(dt,256,256,sampling_args)

figure(2); clf; hold on;
	plot(dt_coeffs_pos,'b-');
	plot(dt_coeffs_neg,'r-');
	plot(dt_coeffs_neg2,'r--');
	axis([1 48 0 20]);

