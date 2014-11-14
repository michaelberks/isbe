% generate a synthetic image with a single line
bg_root = 'U:\projects\mammography\data\synthetic_backgrounds\smooth512\train\';
bg_dir = dir([bg_root,'bg*.mat']);
bg_idx = 1; %ceil(rand*length(bg_dir));
bg_filename = [bg_root,bg_dir(bg_idx).name];
load(bg_filename);

% generate the image
args.contrast_range = [4 4];
args.width_range = 16*[1 1];
args.orientation_range = [120 120];
args.line_type = 'sin';
args.decay_rate = 4;
img = generate_line_image(bg,args);

imagesc(img);

%%	set up arguments for sampling and dt-cwt
args.num_levels = 5;
args.use_nag = 0;

sampling_args.feature_shape = args.line_type;
sampling_args.do_max = 0;
sampling_args.rotate = 0;
sampling_args.win_size = 3;
sampling_args.feature_type = 'conj';


%%	compute dt-cwt coefficients in a 3x3 window around the centre pixel (both
%	in 'conj' representation and 'all' rep.
dt = compute_dual_tree(img, args.num_levels, args.use_nag);
samp1_conj = sample_dt_data(dt, 256, 256, sampling_args);


%%	rotate the image by 180 degrees
img = img(end:-1:1,end:-1:1);

%%	compute dt-cwt coefficients in a 3x3 window around the centre pixel (both
%	in 'conj' representation and 'all' rep.
dt = compute_dual_tree(img, args.num_levels, args.use_nag);
samp2_conj = sample_dt_data(dt, 257, 257, sampling_args);

[samp1_conj; samp2_conj]'
