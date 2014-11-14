args.image_dir = [asymmetryroot,'data/retinograms/DRIVE_clean/training/images/'];
args.vessel_mask_dir = [asymmetryroot,'data/retinograms/DRIVE_clean/training/vessel_masks/'];
args.fov_mask_dir = [asymmetryroot,'data/retinograms/DRIVE_clean/training/fov_masks/'];
args.ori_dir = [asymmetryroot,'data/retinograms/DRIVE_clean/training/orientations/'];
args.width_dir = [asymmetryroot,'data/retinograms/DRIVE_clean/training/width_maps/'];
args.pts_dir = [asymmetryroot,'data/retinograms/DRIVE_clean/training/vxl_rfs/sample_data/'];
args.num_pts = 2000;
args.repeat = 1;

rand('twister', 1000 * args.repeat);

image_list = dir([args.image_dir '*.mat']);
fov_list = dir([args.fov_mask_dir '*.mat']);
vessel_list = dir([args.vessel_mask_dir '*.mat']);
ori_list = dir([args.ori_dir '*.mat']);
width_list = dir([args.width_dir '*.mat']);

num_images = length(image_list);
curr_sample = 1;
%loop through each image sampling data
for i_image = 1:num_images

    display(['Sampling from image ' num2str(i_image)]);

    %Work out the number of samples to take
    num_samples_image = ...
        sample_from_binomial((args.num_pts + 1 - curr_sample), 1/(num_images+1-i_image), 1);

    if ~num_samples_image
        continue;
    end

    %Load in image and masks
    vessel_mask = load_uint8([args.vessel_mask_dir vessel_list(i_image).name]);
    fov_mask = u_load([args.fov_mask_dir fov_list(i_image).name]);
    ori_map = load_uint8([args.ori_dir ori_list(i_image).name]);
    width_map = load_uint8([args.width_dir width_list(i_image).name]);
    ret = rgb2gray(u_load([args.image_dir image_list(i_image).name]));
    centre_mask = bwmorph(vessel_mask, 'thin', 'inf');

    %Check we have enough samples in data
    total_v_pts = sum(vessel_mask(:) & fov_mask(:));
    total_b_pts = sum(~vessel_mask(:) & fov_mask(:));

    %Get random sample of vessel pixels
    v_idx = find(vessel_mask & fov_mask);
    r_idx = randperm(total_v_pts);
    v_idx = v_idx(r_idx(1:num_samples_image));
    [v_rows v_cols] = ind2sub(size(ret), v_idx);

    %Get random sample of background pixels
    b_idx = find(~vessel_mask & fov_mask);
    r_idx = randperm(total_b_pts);
    b_idx = b_idx(r_idx(1:num_samples_image));
    [b_rows b_cols] = ind2sub(size(ret), b_idx);

    true_oris = ori_map(v_idx);
    true_widths = width_map(v_idx);
    true_centres = centre_mask(v_idx);

    %Open a txt file to save points to
    image_path = [args.image_dir 'png/training_' zerostr(i_image, 2) '.png'];
    fid1 = fopen([args.pts_dir '/vessel_xy_' zerostr(i_image,2) '.pts'], 'wt');
    fprintf(fid1, '%s \n', 'version: 1');
    fprintf(fid1, '%s %d \n', 'n_points:', 2*num_samples_image);
    fprintf(fid1, '%s \n', '{'); 
    for i_pt = 1:num_samples_image
        fprintf(fid1,'%.2f %.2f \n', v_cols(i_pt), v_rows(i_pt));
    end
    for i_pt = 1:num_samples_image
        fprintf(fid1,'%.2f %.2f \n', b_cols(i_pt), b_rows(i_pt));
    end
    fprintf(fid1, '%s \n', '}');
    fprintf(fid1, 'nailfold: %s \n', image_path);
    fclose(fid1);
    
    %Open a txt file to save orientations to
    image_path = [args.image_dir 'png/training_' zerostr(i_image, 2) '.png'];
    fid1 = fopen([args.pts_dir 'vessel_oris_' zerostr(i_image,2) '.pts'], 'wt');
    fprintf(fid1, '%s \n', 'version: 1');
    fprintf(fid1, '%s %d \n', 'n_points:', num_samples_image);
    fprintf(fid1, '%s \n', '{'); 
    for i_pt = 1:num_samples_image
        fprintf(fid1,'%.2f %.2f  %.2f %.2f \n', v_cols(i_pt), v_rows(i_pt), real(true_oris(i_pt)), imag(true_oris(i_pt)));
    end
    fprintf(fid1, '%s \n', '}');
    fprintf(fid1, 'nailfold: %s \n', image_path);
    fclose(fid1);
    
    %Open a txt file to save widths to
    image_path = [args.image_dir 'png/training_' zerostr(i_image, 2) '.png'];
    fid1 = fopen([args.pts_dir 'vessel_widths_' zerostr(i_image,2) '.pts'], 'wt');
    fprintf(fid1, '%s \n', 'version: 1');
    fprintf(fid1, '%s %d \n', 'n_points:', num_samples_image);
    fprintf(fid1, '%s \n', '{'); 
    for i_pt = 1:num_samples_image
        fprintf(fid1,'%.2f %.2f  %.2f \n', v_cols(i_pt), v_rows(i_pt), true_widths(i_pt));
    end
    fprintf(fid1, '%s \n', '}');
    fprintf(fid1, 'nailfold: %s \n', image_path);
    fclose(fid1);
    
    %Open a txt file to save centre indices to
    image_path = [args.image_dir 'png/training_' zerostr(i_image, 2) '.png'];
    fid1 = fopen([args.pts_dir 'vessel_centres_' zerostr(i_image,2) '.pts'], 'wt');
    fprintf(fid1, '%s \n', 'version: 1');
    fprintf(fid1, '%s %d \n', 'n_points:', num_samples_image);
    fprintf(fid1, '%s \n', '{'); 
    for i_pt = 1:num_samples_image
        fprintf(fid1,'%.2f %.2f  %.2f \n', v_cols(i_pt), v_rows(i_pt), true_centres(i_pt));
    end
    fprintf(fid1, '%s \n', '}');
    fprintf(fid1, 'nailfold: %s \n', image_path);
    fclose(fid1);

    %Save a PNG copy of the image
    imwrite(ret, image_path);

    %Update count of samples
    curr_sample = curr_sample + num_samples_image;

end
%%
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES=3 LEVELS="[1:16]" NUM_ANGLES=6 EXP_NAME="martha" qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES=3 LEVELS="[1:16]" NUM_ANGLES=6 EXP_NAME="martha" qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh
MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES=3 LEVELS="[1 2 4 8 16]" NUM_ANGLES=6 EXP_NAME="orig" qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES=3 LEVELS="[1 2 4 8 16]" NUM_ANGLES=6 EXP_NAME="orig" qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh

MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES=3 LEVELS="[1:16]" EXP_NAME="martha" qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/compare_DRIVE_g2d_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES=3 LEVELS="[1:16]" EXP_NAME="martha" qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/compare_DRIVE_g2d_csf.sh

MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="[1]" LEVELS="[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]" NUM_ANGLES=6 EXP_NAME="martha" qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/comparing_DRIVE_gabor_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="[1]" LEVELS="[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]" NUM_ANGLES=6 EXP_NAME="martha" qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/comparing_DRIVE_gabor_csf.sh


MAKE_DATA=0 DO_DETECTION=1 DO_WIDTH=0 DO_ORIENTATION=0 NUM_SAMPLES=100000 WIN_SIZES="[3]" LEVELS="[1:16]" EXP_NAME="george" qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/
%%
comparing_DRIVE_gabor_csf(1, ...
    'exp_name',     unixenv('EXP_NAME', []), ...
    'image_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/images/'], ...
    'vessel_mask_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/vessel_masks/'], ...
    'fov_mask_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/fov_masks/'], ...
    'ori_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/orientations/'], ...
    'width_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/width_maps/'], ...
    'selected_images', [],...
    'num_pts',      unixenv('NUM_SAMPLES', 2000), ...
    'num_trees',    unixenv('NUM_TREES', 2), ...
    'make_data',    unixenv('MAKE_DATA',1), ...
    'num_angles',   unixenv('NUM_ANGLES', 6),...
    'levels',       unixenv('LEVELS', [1 2 4]),...
    'do_orientation', unixenv('DO_ORIENTATION',0), ...
    'do_detection', unixenv('DO_DETECTION',0), ...
    'do_width', unixenv('DO_WIDTH',0), ...
    'win_sizes', unixenv('WIN_SIZES',[1 3]) ...
);
%
comparing_DRIVE_experiment_csf(3, 1,  ... % non-strict mode
    'image_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/images/'], ...
    'vessel_mask_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/vessel_masks/'], ...
    'fov_mask_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/fov_masks/'], ...
    'ori_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/orientations/'], ...
    'width_dir',      [asymmetryroot,'data/retinograms/DRIVE/training/width_maps/'], ...
    'selected_images', [],...
    'num_pts',      unixenv('NUM_SAMPLES', 2000), ...
    'num_trees',    unixenv('NUM_TREES', 2), ...
    'make_data',    unixenv('MAKE_DATA',1), ...
    'do_orientation', unixenv('DO_ORIENTATION',0), ...
    'do_detection', unixenv('DO_DETECTION',0), ...
    'do_width', unixenv('DO_WIDTH',0), ...
    'do_tests', unixenv('DO_TESTS',1), ...
    'win_sizes', unixenv('WIN_SIZES',[1 3]) ...
);
%%
clear;
r_all = u_load('C:\isbe\asymmetry_project\experiments\DRIVE\comparing_representations\test\1\responses_gabor.mat');
r_1 = u_load('C:\isbe\asymmetry_project\experiments\DRIVE\comparing_representations\test\1\1\responses_gabor.mat');
r_2 = u_load('C:\isbe\asymmetry_project\experiments\DRIVE\comparing_representations\test\1\2\responses_gabor.mat');
r_4 = u_load('C:\isbe\asymmetry_project\experiments\DRIVE\comparing_representations\test\1\4\responses_gabor.mat');

d_args.decomp_type = {'gabor'};
d_args.num_angles = 18;
d_args.sigma_range = 1;	
d_args.do_max = 0;
d_args.rotate = 0;
d_args.feature_type = 'complex';
d_args.win_size = 3;
d_args.normalise = 0;
d_args.pca = [];

n_args = [];
n_args.bands = 1:3:16;
r1 = convert_decomp_form(r_1, d_args, n_args);
r2 = convert_decomp_form(r_2, d_args, n_args);
r4 = convert_decomp_form(r_4, d_args, n_args);

a_args = d_args;
a_args.sigma_range = [1 2 4];
a_args.num_angles = 6;

n_args = [];
n_args.levels = 1;
r1a = convert_decomp_form(r_all, a_args, n_args);
n_args = [];
n_args.levels = 2;
r2a = convert_decomp_form(r_all, a_args, n_args);
n_args = [];
n_args.levels = 3;
r4a = convert_decomp_form(r_all, a_args, n_args);

max(max(abs(r1a(:) - r1(:))))
max(max(abs(r2a(:) - r2(:))))
max(max(abs(r4a(:) - r4(:))))
%%
d_args{1}.decomp_type = {'g1da'};
d_args{1}.sigma_range = 2;
d_args{1}.num_angles = 1;
d_args{1}.do_max = 0;
d_args{1}.rotate = 0;

d_args{2}.decomp_type = {'g2da'};
d_args{2}.sigma_range = 2;
d_args{2}.num_angles = 1;
d_args{2}.do_max = 0;
d_args{2}.rotate = 0;

d_args{3}.decomp_type = {'gabor'};
d_args{3}.num_angles = 1;
d_args{3}.sigma_range = 2;	
d_args{3}.do_max = 0;
d_args{3}.rotate = 0;
d_args{3}.feature_type = 'complex';

%set common parameters for all decomp types and then compute vector sizes
for ii = 1:3
    d_args{ii}.win_size = 1;
    d_args{ii}.normalise = 0;
    d_args{ii}.pca = [];
end
num_angles = 24;

lrs1 = zeros(31,1);
lrs2 = zeros(31,1);
lrsg = zeros(31,1);

ers1 = zeros(31,1);
ers2 = zeros(31,1);
ersg = zeros(31,1);

for i_angle = -15:15
    line = create_sin_bar(1, 1, -90+i_angle, 16, 16, 0, 8, 8);
    edge = create_sin_step(1, 1, -90+i_angle, 16, 16, 8, 8);
    
    lr1 = compute_filter_responses(line, d_args{1});
    lr2 = compute_filter_responses(line, d_args{2});
    lrg = compute_filter_responses(line, d_args{3});

    er1 = compute_filter_responses(edge, d_args{1});
    er2 = compute_filter_responses(edge, d_args{2});
    erg = compute_filter_responses(edge, d_args{3});
    
    lrs1(i_angle + 16) = lr1.g1d(8,8);
    lrs2(i_angle + 16) = lr2.g2d(8,8);
    lrsg(i_angle + 16) = lrg.gabor(8,8);

    ers1(i_angle + 16) = er1.g1d(8,8);
    ers2(i_angle + 16) = er2.g2d(8,8);
    ersg(i_angle + 16) = erg.gabor(8,8);
    
end
figure; 
subplot(1,2,1); hold all;
plot(-15:15, lrs1);
plot(-15:15, lrs2);
plot(-15:15, sqrt(lrs1.^2 + lrs2.^2));
subplot(1,2,2); hold all;
plot(-15:15, ers1);
plot(-15:15, ers2);
plot(-15:15, sqrt(ers1.^2 + ers2.^2));
figure; 
subplot(1,2,1); hold all;
plot(-15:15, real(lrsg));
plot(-15:15, imag(lrsg));
plot(-15:15, abs(lrsg));
subplot(1,2,2); hold all;
plot(-15:15, real(ersg));
plot(-15:15, imag(ersg));
plot(-15:15, abs(ersg));

figure; 
subplot(2,2,1); imgray(abs(lrg.gabor));
subplot(2,2,2); imgray(abs(erg.gabor));
subplot(2,2,3); imgray(sqrt(lrs1.^2 + lrs2.^2));
subplot(2,2,4); imgray(sqrt(ers1.^2 + ers2.^2));

figure; 
subplot(2,2,1); imgray(angle(lrg.gabor)); colormap(hsv(180)); caxis([-pi pi]);
subplot(2,2,2); imgray(angle(erg.gabor)); colormap(hsv(180)); caxis([-pi pi]);
subplot(2,2,1); imgray(angle(lrg.gabor)); colormap(hsv(180)); caxis([-pi pi]);
subplot(2,2,2); imgray(angle(erg.gabor)); colormap(hsv(180)); caxis([-pi pi]);
subplot(2,2,3); imgray(angle(lrg.gabor)); colormap(hsv(180)); caxis([-pi pi]);
subplot(2,2,4); imgray(angle(erg.gabor)); colormap(hsv(180)); caxis([-pi pi]);
%%
d_args{1}.decomp_type = {'g1da'};
d_args{1}.sigma_range = 2;
d_args{1}.num_angles = 6;
d_args{1}.do_max = 0;
d_args{1}.rotate = 0;

d_args{2}.decomp_type = {'g2da'};
d_args{2}.sigma_range = 2;
d_args{2}.num_angles = 6;
d_args{2}.do_max = 0;
d_args{2}.rotate = 0;

d_args{3}.decomp_type = {'gabor'};
d_args{3}.num_angles = 6;
d_args{3}.sigma_range = 2;	
d_args{3}.do_max = 0;
d_args{3}.rotate = 0;
d_args{3}.feature_type = 'complex';

im = zeros(15);
im(8,8) = 1;
r1 = compute_filter_responses(im, d_args{1});
r2 = compute_filter_responses(im, d_args{2});
rg = compute_filter_responses(im, d_args{3});

figure;
for ii = 1:6
    subplot(2,3,ii); imgray(r1.g1d(:,:,1,ii));
end
figure;
for ii = 1:6
    subplot(2,3,ii); imgray(r2.g2d(:,:,1,ii));
end
figure;
for ii = 1:6
    subplot(2,3,ii); imgray(sqrt(r1.g1d(:,:,1,ii).^2 + r2.g2d(:,:,1,ii).^2));
end
figure;
for ii = 1:6
    subplot(2,3,ii); imgray(abs(rg.gabor(:,:,ii)));
end
%%
for sigma = [1 2 4 8]
    responses = zeros(180,4);
    steered_responses = zeros(180,1);
    for i_angle = 1:180
        edge = create_sin_step(sigma, 1, i_angle, 64, 64, 32, 32);
        [g2dh_responses] = compute_hilbert_2nd_derivatives_sep(edge, sigma);
        responses(i_angle,:) = squeeze(g2dh_responses(32,32,1,:))';
        steered_responses(i_angle,:) =...
            steer_hilbert_2nd_derivatives(g2dh_responses(32,32,1,:), pi * i_angle / 180);   
    end
    figure; hold all;
    plot(responses);
    plot(steered_responses);
end
%%
steered_responses = zeros(16,4);
for i_sigma = 1:16
    edge = create_sin_step(i_sigma, 1, 0, 16, 16, 8, 8);
    [g2dh_responses] = compute_hilbert_2nd_derivatives_sep(edge, [1 2 4 8]);
    steered_responses(i_sigma,:) =...
        squeeze( steer_hilbert_2nd_derivatives(g2dh_responses(8,8,:,:), pi * i_angle / 180) );   
end
figure; hold all;
plot(steered_responses);
%%
responses = zeros(180,4);
for i_angle = 1:180
    theta = pi * i_angle / 180;
    ccc = cos(theta)^3;
    sss = sin(theta)^3;
    ccs3 = 3*cos(theta)^2*sin(theta);
    ssc3 = 3*cos(theta)*sin(theta)^2;
    
    responses(i_angle,1) = ccc;
    responses(i_angle,2) = ccs3;
    responses(i_angle,3) = ssc3;
    responses(i_angle,4) = sss;
    
    
end
figure; hold all;
plot(responses);
%%
im = zeros(32); im(16,16) = 1;
[g2d_responses] = compute_gaussian_2nd_derivatives(im, [1 2 4 8]);
figure;
for ii = 1:3
    for jj = 1:4
        band = g2d_responses(:,:,jj,ii);
        subplot(3,4,4*(ii-1)+jj); imgray(band);
        title(['Max (scaled): ' num2str(max(band(:)) * 2^(2*jj))]);
    end
end
%%
im = zeros(32); im(16,16) = 1;
[h2d_responses] = compute_hilbert_2nd_derivatives_sep(im, [1 2 4 8]);
figure;
for ii = 1:3
    for jj = 1:4
        band = h2d_responses(:,:,jj,ii);
        subplot(3,4,4*(ii-1)+jj); imgray(band);
        title(['Max (scaled): ' num2str(max(band(:)) * 2^(2*jj))]);
        %title(['Max (scaled): ' num2str(max(band(:)))]);
    end
end
%%
im = zeros(32); im(16,16) = 1;
[h2d_responses] = compute_hilbert_2nd_derivatives_sep(im, 2);
[g2d_responses] = compute_gaussian_2nd_derivatives(im, 2);
h2d_responses_a = steer_hilbert_2nd_derivatives(h2d_responses, [], 6);
g2d_responses_a = steer_gaussian_2nd_derivatives(g2d_responses, [], 6);

abs_responses = sqrt(h2d_responses_a.^2 + g2d_responses_a.^2);
figure;
for i_angle = 1:6
    subplot(2,3,i_angle); imgray(h2d_responses_a(:,:,:,i_angle));
end
figure;
for i_angle = 1:6
    subplot(2,3,i_angle); imgray(g2d_responses_a(:,:,:,i_angle));
end
figure;
for i_angle = 1:6
    subplot(2,3,i_angle); imgray(abs_responses(:,:,:,i_angle));
end
%
figure; hold all;
plot(h2d_responses_a(8,:,:,1));
plot(g2d_responses_a(8,:,:,1));
plot(abs_responses(8,:,:,1));
%%
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="[3]" LEVELS="[1:16]" NUM_ANGLES=6 EXP_NAME="scales" qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh
MAKE_DATA=0 DO_DETECTION=0 DO_WIDTH=0 DO_ORIENTATION=1 NUM_SAMPLES=100000 WIN_SIZES="[3]" LEVELS="[1 2 4 8 16]" NUM_ANGLES=18 EXP_NAME="angles" qsub -V -t 1 -l twoday matlab_code/trunk/hydra/cuc/compare_DRIVE_gabor_csf.sh
%%
g2d_responses_a = zeros(180,6);
gabor_responses_a = zeros(180,6);
dt_responses_a = zeros(180,6);
for i_ori = 1:180
    
    im = create_gauss_bar(1,1,i_ori,32,32,16,16);
    
    [g2d_responses] = compute_gaussian_2nd_derivatives(im, 2);
    g2d_responses_a(i_ori,:) = ...
        squeeze( steer_gaussian_2nd_derivatives(g2d_responses(16,16,1,:), 6) )';
    
    [gabor_responses] = compute_gabor_responses(im, 2, 6);
    gabor_responses_a(i_ori,:) = ...
        squeeze( real(gabor_responses(16,16,1,:)) )';
    
    dt = compute_dual_tree(im, 2, 0);
    dt_responses_a(i_ori,:) = ...
        squeeze( abs(sample_dt_data(dt, 16, 16, 'levels', 2, 'win_size', 1)) )';
end
figure; plot(g2d_responses_a);
figure; plot(gabor_responses_a);
figure; plot(dt_responses_a);
%%
g2d_responses_a = zeros(32,4);
gabor_responses_a = zeros(32,4);
dt_responses_a = zeros(32,4);
for i_width = 1:32
    
    im = create_gauss_bar(i_width/2,1,67,32,32,16,16);
    
    [g2d_responses] = compute_gaussian_2nd_derivatives(im, [1 2 4 8]);
    g2d_responses_a(i_width,:) = ...
        squeeze( abs( steer_gaussian_2nd_derivatives(g2d_responses(16,16,:,:), 1)) )';

    [gabor_responses] = compute_gabor_responses(im, [1 2 4 8], 1);
    gabor_responses_a(i_width,:) = ...
        squeeze( real(gabor_responses(16,16,:,1)) )';
    
    dt = compute_dual_tree(im, 4, 0);
    temp = squeeze( abs(sample_dt_data(dt, 16, 16, 'levels', 1:4, 'win_size', 1)) );
    dt_responses_a(i_width,:) = temp(3,:);
        
end
figure; plot(g2d_responses_a);
figure; plot(gabor_responses_a);
figure; plot(dt_responses_a);
%%
oris = 4:4:180;
widths = (1:32)/2;

g2d_responses_a = zeros(45, 32, 6, 5);
gabor_responses_a = zeros(45, 32, 6, 5);
dt_responses_a = zeros(45, 32, 6, 5);
for i_ori = 1:45
    for i_width = 1:32
        
        im = create_gauss_bar(widths(i_width)/2,1,oris(i_ori),64,64,32,32);
        
        [g2d_responses] = compute_gaussian_2nd_derivatives(im, [1 2 4 8 16]);
        g2d_responses_a(i_ori,i_width,:,:) = ...
            permute( abs( steer_gaussian_2nd_derivatives(g2d_responses(32,32,:,:), 6)), [1 2 4 3]);

        [gabor_responses] = compute_gabor_responses(im, [1 2 4 8 16], 6);
        gabor_responses_a(i_ori,i_width,:,:) = ...
            permute( abs(gabor_responses(32,32,:,:)), [1 2 4 3]);
        
        dt = compute_dual_tree(im, 5, 0);
        dt_responses_a(i_ori,i_width,:,:) = ...
            abs(sample_dt_data(dt, 32, 32, 'levels', 1:5, 'win_size', 1));
    end
end
colors = 'rgbymc';
%%
figure; a1 = gca; hold all;
figure; a2 = gca; hold all;
figure; a3 = gca; hold all;

for i_scale = 1:5
    for i_band = 1:6    
        axes(a1); surface(widths, oris, g2d_responses_a(:,:,i_band,i_scale), 'edgecolor','none','facecolor', 'interp');%'edgecolor', colors(i_band), 
        axes(a2); surface(widths, oris, gabor_responses_a(:,:,i_band,i_scale), 'edgecolor','none','facecolor', 'interp');
        axes(a3); surface(widths, oris, dt_responses_a(:,:,i_band,i_scale), 'edgecolor','none','facecolor', 'interp');
    end
    colors = circshift(colors,[0 1]);
end
%%
save C:\isbe\asymmetry_project\experiments\filter_responses_surfaces.mat *_responses_a
%%
oris = 4:4:180;
widths = (1:32)/2;

mono_responses_a = zeros(45, 32, 5);
for i_ori = 1:45
    for i_width = 1:32
        
        im = create_gauss_bar(widths(i_width)/2,1,oris(i_ori),64,64,32,32);
        
        [local_amp] = ...
                monogenic(im, 5, 2, 2, 0.65, 1);
        
        mono_responses_a(i_ori,i_width,:) = local_amp(32,32,2:end);

    end
end
figure;
for i_scale = 1:5
    surface(widths, oris, mono_responses_a(:,:,i_scale), 'edgecolor','none','facecolor', 'interp');
end   
%%
build_predictor(...
    'predictor_name',       unixenv('PREDICTOR_NAME','predictor'), ...
    'task_id',				unixenv('SGE_TASK_ID',1), ...
    'job_id',				'mike', ...
    ... % Folders
    'image_root',           [asymmetryroot,'data/retinograms/STARE/training'],...
    'model_root',           [asymmetryroot,'data/models/vessel'], ...
	... % Output parameters
    'output_type',          unixenv('OUTPUT_TYPE', 'detection'), ...
	... % Sampling parameters
    'num_samples',			unixenv('NUM_SAMPLES',1000), ...
    'max_n_images',         unixenv('MAX_N_IMAGES',[]), ...
    'bg_ratio',				unixenv('BG_RATIO',1), ...
    'replace_sample',       unixenv('REPLACE_SAMPLE', false), ...
    'sampling_method',      unixenv('SAMPLING_METHOD','generate_training_data'), ...
    'shrink_fov',           unixenv('SHRINK_FOV',false), ...
    'image_type', 			unixenv('IMAGE_TYPE','real'), ...
        ... % Image sampling parameters
        'image_dir',            unixenv('IMAGE_DIR', 'images'),...
        'fov_mask_dir',         unixenv('FOV_MASK_DIR', 'fov_masks'),...
        'fg_mask_dir',          unixenv('FG_MASK_DIR', 'vessel_masks'),...
        'prediction_dir',       unixenv('PREDICTION_DIR', 'predictions'),...
        'probability_dir',      unixenv('PROBABILITY_DIR', ''),...
        'make_resampling_maps',          unixenv('MAKE_RESAMPLING_MAPS', 0.5), ...
        ... % Saved training data parameters
        'save_training_data',	unixenv('SAVE_TRAINING_DATA',false), ...
        'training_data_dir', 	unixenv('TRAINING_DATA_DIR','saved_training_data'), ...
        'training_data',		unixenv('TRAINING_DATA',''), ... % ???
        'training_labels',		unixenv('TRAINING_LABELS',''), ... % ???
	... % Image feature/decomposition parameters
    'num_levels', 			unixenv('NUM_LEVELS',1), ...
    'rgb_channel',          unixenv('RGB_CHANNEL','rgb'), ...
    'normalise', 			unixenv('NORMALISE',0), ...
    'win_size',				unixenv('WIN_SIZE',1), ...
    'pca_filename',         unixenv('PCA_FILENAME',[]), ...
    'do_max',				unixenv('DO_MAX',false), ...
    'rotate',				unixenv('ROTATE',false), ...
    'decomp_type', 			unixenv('DECOMP_TYPE','gabor'), ...
        ... % DTCWT parameters
        'feature_shape', 		unixenv('FEATURE_SHAPE','rect'), ...
        'feature_type',			unixenv('FEATURE_TYPE','conj'), ...
        ... % Gaussian derivative parameters
        'sigma_range', 			unixenv('SIGMA_RANGE',[1,2,4,8]), ...
    ... % Predictor parameters
    'prediction_type',		unixenv('PREDICTION_TYPE','rf_classification'), ...
        ... % Tree/Forest parameters
        'n_trees',				unixenv('NUM_TREES',2), ...
        'split_criterion_c',    unixenv('SPLIT_CRITERION_C','gdi'),...
        'split_criterion_r',    unixenv('SPLIT_CRITERION_R','ssq'),...
        'var_criterion_c',		unixenv('VAR_CRITERION_C','mabs'),...
        'var_criterion_r',		unixenv('VAR_CRITERION_R','ssq'),...
        'split_min',			unixenv('SPLIT_MIN',10), ...
        'end_cut_min',			unixenv('END_CUT_MIN',1), ...
        'do_ubound',			unixenv('DO_UBOUND',1), ...
        'do_circular',			unixenv('DO_CIRCULAR',[]), ...
        'w_prior',				unixenv('W_PRIOR',0), ...
        'impure_thresh',		unixenv('IMPURE_THRESH',1e-4), ...
        'minimise_size',		unixenv('MINIMIZE_TREE',0), ...
        'd',                    unixenv('d',[]), ...
        'quiet',                unixenv('QUIET', 0),...
    ... % Miscellaneous parameters
    'overwrite',			unixenv('OVERWRITE',false), ...
    'use_nag',				unixenv('USE_NAG',false),...
    'rand_seed',			unixenv('RAND_SEED',[]) ...
);
%%
for ii = 1:20
    load(['C:\isbe\asymmetry_project\data\retinograms\STARE\training\fov_masks\' zerostr(ii,2) '_training_f_mask.mat'])
%     f_mask_c = imclose(f_mask, strel('disk', 20));
%     f_mask_o = imopen(f_mask, strel('disk', 20));
% %     figure; 
% %     subplot(1,2,1); imgray(f_mask);
% %     subplot(1,2,2); imgray(f_mask_o);
%     f_mask = f_mask_o;
    f_mask([1 end],:) = 0;
    f_mask(:,[1 end]) = 0;
    save(['C:\isbe\asymmetry_project\data\retinograms\STARE\training\fov_masks\' zerostr(ii,2) '_training_f_mask.mat'], 'f_mask');

end
