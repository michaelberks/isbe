frame_dir = 'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\002wellcome\2015_02_27\L4_11_00_04\sequence_data\';
frame = imread([frame_dir, 'stationary_mosaic02.png']);
frame = double(frame(33:end-32, 33:end-32));
frame = 255*(frame - min(frame(:))) / (max(frame(:)) - min(frame(:)));
%%
% Get struct "o" containing SR data, and "gtruth" ground truth image.
K = 100;
gamma = 0.4;
sigma_noise = 5/255;
zm = 2;
transform = 'homography';
[o_gt,gtruth] = synthdata_demo([], K, transform, gamma, sigma_noise, zm); 
gtruth2 = imresize(gtruth, 0.5);

figure;
for i = 1:min(K,6)
    subplot(2,3,i); imgray(o_gt(i).im);
    title(['low-res ' num2str(i)]);
end

[ni_hi, nj_hi] = size(gtruth);
[~, theta, lambda] = vectorise_lo_res_data(o_gt);
theta_offset = mean(theta, 2);
theta_scaling = 0.35 ./ std(theta, 1, 2);
lambda_offset = mean(lambda, 2);
lambda_scaling = 0.35 ./ std(lambda, 1, 2);
%%
corr_points_x = zeros(K, 5);
corr_points_y = zeros(K, 5);
figure;
for i_k = 1:K
    clf;
    imgray(o_gt(i_k).im);
    title(i_k);
    x = [];
    y = [];
    while length(x) ~= 7
        [x, y, ~] = impixel();
    end
    corr_points_x(i_k, :) = x(:)';
    corr_points_y(i_k, :) = y(:)';
end
figure; plot(corr_points_x, corr_points_y, '.');
%%
figure; imgray(gtruth);
x_gt = [];
y_gt = [];
while length(x_gt) ~= 5
    [x_gt, y_gt, ~] = impixel();
end
%%
xy_gt = [x_gt y_gt];
o_man = o_gt;
theta_man = zeros(8, K);
for i_k = 1:K

    xy_lr = [corr_points_x(i_k,:)' corr_points_y(i_k,:)'];

    tform = fitgeotrans(xy_lr,xy_gt,'projective');
    [x_lrt y_lrt] = transformPointsForward(tform, xy_lr(:,1), xy_lr(:,2));
    
    o_man(i_k).H = tform.T';
    theta_man(:, i_k) = tform.T(1:8);
    
    figure;
    subplot(1,2,1); imgray(o_gt(i_k).im);
    plot(xy_lr(:,1), xy_lr(:,2), 'rx');
    
    subplot(1,2,2); imgray(gtruth);
    plot(xy_gt(:,1), xy_gt(:,2), 'rx');
    plot(x_lrt, y_lrt, 'go');
    
end
%%
%%
K = 20;
[Y, theta, lambda, gamma, ni_lo, nj_lo] = vectorise_lo_res_data(o_gt(1:K));
[avg_im_gt] = getAvim(ni_hi, nj_hi, ni_lo, nj_lo, Y, theta, gamma);
[W_gt] = makeW(ni_hi, nj_hi, ni_lo, nj_lo, theta, gamma);
%%
num_frames = K;
[num_rows, num_cols] = size(o_gt(1).im);
frames = zeros(num_rows, num_cols, num_frames);
gt_trans = zeros(3,3,num_frames);
rev_trans = zeros(3,3,num_frames);

cx = (1+num_cols)/2;
cy = (1+num_rows)/2;
offset_corrections = zeros(num_frames,2);
for i_frame = 1:num_frames
    gt_trans(:,:,i_frame) = o_gt(i_frame).H / 2;
    gt_trans(3,3,i_frame) = 1;
    offset_corrections(i_frame,:) = [cx cy] - [cx cy]*gt_trans(1:2,1:2,i_frame);
    gt_trans(1:2,3,i_frame) = gt_trans(1:2,3,i_frame) + offset_corrections(i_frame,:)';
    rev_trans(:,:,i_frame) = inv(gt_trans(:,:,i_frame));
    frames(:,:,i_frame) = o_gt(i_frame).im;
end

[gt_avg] = create_mosaic(frames, gt_trans);
figure; imgray(gt_avg);   
%%
KM = sum(ni_lo.*nj_lo);
lambda1 = ones(KM,1);
lambda2 = ones(KM,1);
curr_pt = 0;
for i_k = 1:K
    n_pts = v_lo(i_k)*h_lo(i_k);
    idx = curr_pt + (1:n_pts);
    curr_pt = curr_pt + n_pts;
    lambda1(idx) = lambda(1,i_k);
    lambda2(idx) = lambda(2,i_k);
end
%
alpha = 0.08;
nu = 0.04;

huber_im_gt = superres_huber(W_gt, Y, lambda1, lambda2, avg_im_gt, alpha, nu);
ml_im_gt = superres_ml(W_gt, Y, lambda1, lambda2, avg_im_gt);%

figure; 
subplot(1,3,1); imgray(avg_im_gt);
subplot(1,3,2); imgray(huber_im_gt);
subplot(1,3,3); imgray(ml_im_gt);
%%   
[init_reg, match_counts] = ...
        register_tiles_features(frames, ...
                            'ref_type', 'mosaic',...
                            'region', 'centre',...
                            'theta_range', -5:0.5:5, ...
                            'offset_lim', 50, ...
                            'mosaic', gtruth2,...
                            'sigma', 2,...
                            'tile_masks', [],...
                            'debug', 1);
                        
[init_avg, ~, t2m_reg ] = create_mosaic(frames, init_reg);
figure; imgray(init_avg);
%%
o_reg = o_gt;
for i_frame = 1:num_frames
    o_reg(i_frame).H = init_reg(:,:,i_frame);
    o_reg(i_frame).H(1:2,3) = o_reg(i_frame).H(1:2,3) - ([cx cy] - [cx cy]*init_reg(1:2,1:2,i_frame))';
    o_reg(i_frame).H(1:2,1:3) = o_reg(i_frame).H(1:2,1:3)*2;
    o_reg(i_frame).g = 0.4;
    
    o_reg(i_frame).l1 = 1;
    o_reg(i_frame).l2 = 0;
    
    %o2(i_frame).H(1:2,3) = o(i_frame).H(1:2,3);
end
%%
[Y, theta_reg, lambda_reg, gamma_reg, ni_lo, nj_lo] = vectorise_lo_res_data(o_reg(1:K));
[avg_im_reg, mask] = getAvim(ni_hi, nj_hi, ni_lo, nj_lo, Y, theta_reg, gamma_reg);
[W_reg] = makeW(ni_hi, nj_hi, ni_lo, nj_lo, theta_reg, gamma_reg);
%%
KM = sum(ni_lo.*nj_lo);
lambda1 = ones(KM,1);
lambda2 = ones(KM,1);
curr_pt = 0;
for i_k = 1:K
    n_pts = v_lo(i_k)*h_lo(i_k);
    idx = curr_pt + (1:n_pts);
    curr_pt = curr_pt + n_pts;
    lambda1(idx) = 1;
    lambda2(idx) = 0;
end
%
alpha = 0.08;
nu = 0.04;

huber_im_reg = superres_huber(W_reg, Y, lambda1, lambda2, avg_im_reg, alpha, nu, 3);
ml_im_reg = superres_ml(W_reg, Y, lambda1, lambda2, avg_im_reg, 3);%

figure; 
subplot(1,3,1); imgray(avg_im_reg);
subplot(1,3,2); imgray(huber_im_reg); caxis([-.5 .5]);
subplot(1,3,3); imgray(ml_im_reg); caxis([-.5 .5]);
%%
[W_reg] = makeW(ni_hi, nj_hi, ni_lo, nj_lo, theta_man, gamma_reg);
[avg_im_man] = getAvim(ni_hi, nj_hi, ni_lo, nj_lo, Y, theta_man, gamma_reg);    
huber_im_man = superres_huber(W_reg, Y, lambda1, lambda2, avg_im_man, alpha, nu, 3);
ml_im_man = superres_ml(W_reg, Y, lambda1, lambda2, avg_im_man, 3);%

figure; 
subplot(1,3,1); imgray(avg_im_man);
subplot(1,3,2); imgray(huber_im_man); caxis([-.5 .5]);
subplot(1,3,3); imgray(ml_im_man); caxis([-.5 .5]);
%%
[hi_res_img, lambda_out, theta_out, alpha_out, beta_out] = ...
    super_res_simultaneous(Y, theta, lambda_reg, gamma_reg, ni_lo, nj_lo, ni_hi, nj_hi,...
    'theta_offset', theta_offset,...
    'theta_scaling', theta_scaling, ...
    'lambda_offset', lambda_offset,...
    'lambda_scaling', lambda_scaling, ...
    'n_itr_all', 5,...
	'n_itr_alpha', 2,...
    'n_itr_x', 4);

figure; imgray(hi_res_img);
%%
theta_sigma = 0.35 ./ theta_scaling;

for i_noise = [1 2 5 10 20 100 1000]
    theta_noise = bsxfun(@times, randn(size(theta))/i_noise, theta_sigma);
    theta_degraded = theta + theta_noise;
    [avg_im_degraded] = getAvim(ni_hi, nj_hi, ni_lo, nj_lo, Y, theta_degraded, gamma_reg);
    figure; 
    subplot(1,2,1); imgray(avg_im_gt);
    subplot(1,2,2); imgray(avg_im_degraded);
    title(['noise = 1 / ' num2str(i_noise) ' s.d.']);
end
%%
theta_sigma = 0.35 ./ theta_scaling;

for i_noise = [1 2 5 10 20 100 1000]
    theta_noise = bsxfun(@times, randn(size(theta))/i_noise, theta_sigma);
    theta_degraded = theta + theta_noise;
    [avg_im_degraded] = getAvim(ni_hi, nj_hi, ni_lo, nj_lo, Y, theta_degraded, gamma_reg);
    [W_reg] = makeW(ni_hi, nj_hi, ni_lo, nj_lo, theta_degraded, gamma_reg);
    huber_im_degraded = superres_huber(W_reg, Y, lambda1, lambda2, avg_im_degraded, alpha, nu, 4);
    figure; 
    subplot(1,2,1); imgray(huber_im_gt);
    subplot(1,2,2); imgray(huber_im_degraded);
    title(['noise = 1 / ' num2str(i_noise) ' s.d.']);
end
%%
theta_noise = bsxfun(@times, randn(size(theta))/10, theta_sigma);
theta_degraded = theta + theta_noise;
[hi_res_img_10, lambda_out, theta_out, alpha_out, beta_out] = ...
    super_res_simultaneous(Y, theta_degraded, lambda_reg, gamma_reg, ni_lo, nj_lo, ni_hi, nj_hi,...
    'theta_offset', theta_offset,...
    'theta_scaling', theta_scaling, ...
    'theta_gt', theta,...
    'lambda_offset', lambda_offset,...
    'lambda_scaling', lambda_scaling, ...
    'n_itr_all', 5,...
	'n_itr_alpha', 2,...
    'n_itr_x', 4);
%%
[W_reg] = makeW(ni_hi, nj_hi, ni_lo, nj_lo, theta_degraded, gamma_reg);
huber_im_degraded = superres_huber(W_reg, Y, lambda1, lambda2, avg_im_degraded, alpha, nu, 4);

figure; 
subplot(1,2,1); imgray(huber_im_degraded);
subplot(1,2,2); imgray(hi_res_img_10);
%%
M = ni_lo(1)*nj_lo(1);
n_samples = ceil(0.01 * M);
[potential_cols, potential_rows] = meshgrid(5:nj_lo(1)-4, 5:ni_lo(1)-4);
potential_idx = sub2ind([ni_lo(1) nj_lo(1)], potential_rows(:), potential_cols(:));
r_samples = randperm(length(potential_idx), n_samples);
val_px = false(M,1);
val_px(potential_idx(r_samples)) = 1;
    
profile on;
W_kd = update_Wk(ni_hi, nj_hi, ni_lo(1), nj_lo(1), theta_reg(:,1), gamma(1), val_px);
profile viewer;
%%
angles = zeros(num_frames,1);
frames = zeros(size(gtruth,1), size(gtruth,2), num_frames);
frames(:,:,1) = gtruth;
for i_frame = 2:num_frames
    angles(i_frame) = round(10*rand-5);
    frames(:,:,i_frame) = double(imrotate(gtruth, angles(i_frame), 'crop'));
end
%%
[init_reg, match_counts] = ...
        register_tiles_features(frames, ...
                            'ref_type', 'previous',...
                            'region', 'centre',...
                            'theta_range', -5:5, ...
                            'offset_lim', 20, ...
                            'mosaic', [],...
                            'sigma', 6,...
                            'tile_masks', [],...
                            'debug', 1);
                        
[init_avg, ~, t2m_reg ] = create_mosaic(frames, init_reg);
        
figure; imgray(init_avg);

pred_angles = asind(squeeze(init_reg(2,1,:)));
[round(angles) pred_angles]
%%
frames_dir = 'C:\isbe\nailfold\camera_capture\003cic_2016\2016_12_06\L4_11_44_42\';
num_rows = 480;
num_cols = 640;
num_frames = 20;
start_row = 250;
start_col = 300;
num_rows_lo = 128;
num_cols_lo = 128;

frames_orig = zeros(num_rows,num_cols,3,num_frames, 'uint8');
frames = zeros(num_rows,640,num_frames);
frames_eq = zeros(num_rows,640,num_frames);
tile_masks = false(num_rows,640,num_frames);
frames_lo = zeros(num_rows_lo, num_cols_lo,num_frames);

for i_f = 1:num_frames
    frame = imread([frames_dir 'frame' zerostr(i_f,5) '.bmp']);
    frame_gray = rgb2gray(frame);
    
    glare_mask = all(frame > 150,3);
    glare_mask = imdilate(glare_mask, strel('disk', 20));
    frame_eq_vals = histeq(frame_gray(~glare_mask), 128);
    frame_eq = frame_gray;
    frame_eq(~glare_mask) = frame_eq_vals;

    frames_orig(:,:,:,i_f) = frame;
    frames(:,:,i_f) = double(double(frame_gray));
    frames_eq(:,:,i_f) = double(double(frame_eq)/255 - 0.5);
    tile_masks(:,:,i_f) = ~glare_mask;
    
    frames_lo(:,:,i_f) = frames_eq(start_row+(1:num_rows_lo), start_col+(1:num_cols_lo),i_f);
end
%%
[init_reg, match_counts] = register_tiles_features(frames, ...
    'ref_type', 'previous',...
    'region', 'centre',...
    'theta_range', -5:5, ...
    'offset_lim', 40, ...
    'mosaic', [],...
    'sigma', [1 2],...
    'tile_masks', tile_masks,...
    'debug', 1);
[init_avg] = create_mosaic(frames_eq, init_reg);

figure; imgray(init_avg);
%%
[init_reg, match_counts] = register_tiles_features(frames_lo, ...
    'ref_type', 'previous',...
    'region', 'centre',...
    'theta_range', -5:5, ...
    'offset_lim', 40, ...
    'mosaic', [],...
    'sigma', [1 2],...
    'tile_masks', [],...
    'debug', 1);
[init_avg, ~, t2m_reg ] = create_mosaic(frames_lo, init_reg);

figure; imgray(init_avg);
%%
zoom = 2;
cx = (1+num_cols_lo)/2;
cy = (1+num_rows_lo)/2;
num_rows_hi = zoom*num_rows_lo;
num_cols_hi = zoom*num_cols_lo;

for i_frame = 1:num_frames
    %Set fields 'H', 'la', 'lb', 'g', 'im'
    o_gt(i_frame).H = init_reg(:,:,i_frame);
    o_gt(i_frame).H(1:2,3) = o_gt(i_frame).H(1:2,3) - ([cx cy] - [cx cy]*init_reg(1:2,1:2,i_frame))';
    o_gt(i_frame).H(1:2,1:3) = o_gt(i_frame).H(1:2,1:3)*2;
    
    o_gt(i_frame).la = 1;
    o_gt(i_frame).b = 0;
    o_gt(i_frame).g = 0.4;
    
    o_gt(i_frame).im = frames_lo(:,:,i_frame);
    
end

[avg_im_reg,msk] = getAvim(num_rows_hi, num_cols_hi, o_gt);

[W_reg, Y, La, Lb] = makeW(num_rows_hi, num_cols_hi, o_gt);
opts = zeros(1,18); % This is the "options" vector for the Netlab "scg" routine.
opts(1) = 1; % verbose
opts(2:3) = 1e-3; % convergeance criteria
opts(14) = 50; % number of iterations before automatic termination.
alp = 0.08;
nu = 0.04;

huber_im_reg = superres_huber(W_reg,Y,La,Lb,avg_im_reg,alp,nu,opts);

figure; 
subplot(1,2,1); imgray(init_avg);
subplot(1,2,2); imgray(avg_im_reg);

figure; 
subplot(1,2,1); imgray(avg_im_reg);
subplot(1,2,2); imgray(huber_im_reg); caxis([-.5 .5]);
%%
frames_dir = 'C:\isbe\nailfold\camera_capture\003cic_2016\2016_12_06\L4_11_48_26\';

num_frames = 20;
start_row = 1000;
start_col = 750;
num_rows_lo = 512;
num_cols_lo = 1024;

frames_lo = zeros(num_rows_lo, num_cols_lo,num_frames);

for i_f = 1:num_frames
    frame = imread([frames_dir 'frame' zerostr(i_f,5) '.bmp']);
    frame_gray = rgb2gray(frame);
    
    glare_mask = all(frame > 150,3);
    glare_mask = imdilate(glare_mask, strel('disk', 20));
    frame_eq_vals = histeq(frame_gray(~glare_mask), 128);
    frame_eq = frame_gray;
    frame_eq(~glare_mask) = frame_eq_vals;

    frames_lo(:,:,i_f) = double(frame_eq(start_row+(1:num_rows_lo), start_col+(1:num_cols_lo)))/255 - 0.5;
end
%%
[init_reg, match_counts] = register_tiles_features(frames_lo, ...
    'ref_type', 'previous',...
    'region', 'centre',...
    'theta_range', linspace(-5, 5, 41), ...
    'offset_lim', 40, ...
    'mosaic', [],...
    'sigma', 8,...
    'tile_masks', [],...
    'debug', 1);
[init_avg, ~, t2m_reg ] = create_mosaic(frames_lo, init_reg);

figure; imgray(init_avg);
%%
zoom = 2;
cx = (1+num_cols_lo)/2;
cy = (1+num_rows_lo)/2;

sampled_rows_lo = 400;
sampled_cols_lo = 400;
num_rows_hi = zoom*sampled_rows_lo;
num_cols_hi = zoom*sampled_cols_lo;

for i_frame = 1:num_frames
    %Set fields 'H', 'la', 'lb', 'g', 'im'
    o_gt(i_frame).H = init_reg(:,:,i_frame);
    o_gt(i_frame).H(1:2,3) = o_gt(i_frame).H(1:2,3) - ([cx cy] - [cx cy]*init_reg(1:2,1:2,i_frame))';
    o_gt(i_frame).H(1:2,1:3) = o_gt(i_frame).H(1:2,1:3)*2;
    
    o_gt(i_frame).la = 1;
    o_gt(i_frame).lb = 0;
    o_gt(i_frame).g = 0.4;
    
    o_gt(i_frame).im = frames_lo(1:sampled_rows_lo,1:sampled_cols_lo,i_frame);
    
end

[avg_im_reg,msk] = getAvim(num_rows_hi, num_cols_hi, o_gt);
figure; imgray(avg_im_reg);
%%
[W_reg, Y, La, Lb] = makeW(num_rows_hi, num_cols_hi, o_gt);
%%
opts = zeros(1,18); % This is the "options" vector for the Netlab "scg" routine.
opts(1) = 1; % verbose
opts(2:3) = 1e-3; % convergeance criteria
opts(14) = 50; % number of iterations before automatic termination.
alp = 0.08;
nu = 0.04;

huber_im_reg = superres_huber(W_reg,Y,La,Lb,avg_im_reg,alp,nu,opts);

figure; 
subplot(1,2,1); imgray(init_avg);
subplot(1,2,2); imgray(avg_im_reg);

figure; 
subplot(1,2,1); imgray(avg_im_reg);
subplot(1,2,2); imgray(huber_im_reg); caxis([-.5 .5]);

%%
figure; 
a1 = subplot(1,2,2); imgray(o_gt(1).im);
apex_pts = zeros(4,2,num_frames);
for i_fr = 1:num_frames
    subplot(1,2,1); imgray(o_gt(i_fr).im);
    title('Select 4 points on the image');
    xi = [];
    while length(xi) ~= 4
        [xi,yi,~] = impixel();
    end 
    apex_pts(:,:,i_fr) = [xi yi];
    
    subplot(1,2,2); imgray(o_gt(i_fr).im);
    plot(a1, xi, yi, 'rx');
end
%%
frames_dir = 'C:\isbe\nailfold\camera_capture\cic_2016\003cic_2016\2016_12_06\L4_11_44_42\';
K = 20;
frames = zeros(240, 320, K);
for i_k = 1:K
    frame = imread([frames_dir 'frame' zerostr(i_k,5) '.bmp']);
    frame = double(rot90(rgb2gray(frame),2));
     
     frames(:,:,i_k) = frame(1:240,120+(1:320));
end
%%
n_pts = 10;
corr_points_x = zeros(n_pts, K);
corr_points_y = zeros(n_pts, K);
figure;
for i_k = 1:K
    clf;
    if i_k > 1
        subplot(1,2,1); imgray(frames(:,:,i_k-1));
        plot(corr_points_x(:,i_k-1), corr_points_y(:,i_k-1), 'g*');
    end
    subplot(1,2,2);  imgray(frames(:,:,i_k));
    title(i_k);
    if i_k > 1
        plot(corr_points_x(:,i_k-1), corr_points_y(:,i_k-1), 'g*');
    end
    
    x = [];
    y = [];
    while length(x) ~= n_pts
        [x, y, ~] = impixel();
    end
    corr_points_x(:,i_k) = x(:);
    corr_points_y(:,i_k) = y(:);
end
figure; axis image xy; plot(corr_points_x', corr_points_y', '.');
%%
xy_gt = 2*[corr_points_x(:,1) corr_points_y(:,1)];
theta_man = zeros(8, K);
gamma = 0.4*ones(1,K);

gtruth = imresize(frames(:,:,1), 2);

ni_hi = 480;
nj_hi = 640;

ni_lo = 240*ones(K,1);
nj_lo = 320*ones(K,2);

for i_k = 1:K

    xy_lr = [corr_points_x(:,i_k) corr_points_y(:,i_k)];

    tform = fitgeotrans(xy_lr,xy_gt,'projective');
    [x_lrt y_lrt] = transformPointsForward(tform, xy_lr(:,1), xy_lr(:,2));
    
    o_man(i_k).H = tform.T';
    theta_man(:, i_k) = tform.T(1:8);
    
    figure;
    subplot(1,2,1); imgray(frames(:,:,i_k));
    plot(xy_lr(:,1), xy_lr(:,2), 'rx');
    
    subplot(1,2,2); imgray(gtruth);
    plot(xy_gt(:,1), xy_gt(:,2), 'rx');
    plot(x_lrt, y_lrt, 'go');
    
end
    
[avg_im_gt] = getAvim(ni_hi, nj_hi, ni_lo, nj_lo, frames(:), theta_man, gamma);
[W_gt] = makeW(ni_hi, nj_hi, ni_lo, nj_lo, theta, gamma);
%%
frame_dir = 'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\002wellcome\2015_02_27\L4_11_00_04\sequence_data\';
frame = imread([frame_dir, 'stationary_mosaic02.png']);
frame = double(frame(80+(1:240), 200+(1:320)));
frame = 255*(frame - min(frame(:))) / (max(frame(:)) - min(frame(:)));
%%
% Get struct "o" containing SR data, and "gtruth" ground truth image.
K = 20;
gamma = 0.4;
sigma_noise = 5/255;
zm = 4;
transform = 'homography';
[o_gt,gtruth] = synthdata_demo(frame, K, transform, gamma, sigma_noise, zm); 
gtruth2 = imresize(gtruth, 0.5);

figure;
for i = 1:min(K,6)
    subplot(2,3,i); imgray(o_gt(i).im);
    title(['low-res ' num2str(i)]);
end

[ni_hi, nj_hi] = size(gtruth);
[Y, theta_gt, lambda_gt, gamma, ni_lo, nj_lo] = vectorise_lo_res_data(o_gt);
%theta_offset = mean(theta, 2);
%theta_scaling = 0.35 ./ std(theta, 1, 2);
%lambda_offset = mean(lambda, 2);
%lambda_scaling = 0.35 ./ std(lambda, 1, 2);
%%
n_pts = 4;
corr_points_x = zeros(n_pts, K);
corr_points_y = zeros(n_pts, K);
figure;
for i_k = 1:K
    clf;
    imgray(o_gt(i_k).im);
    title(i_k);
    
    
    x = [];
    y = [];
    while length(x) ~= n_pts
        [x, y, ~] = impixel();
    end
    corr_points_x(:,i_k) = x(:);
    corr_points_y(:,i_k) = y(:);
end
figure; axis image xy; plot(corr_points_x', corr_points_y', '.'); 
%%
figure; imgray(gtruth);
x_gt = [];
y_gt = [];
while length(x_gt) ~= n_pts
    [x_gt, y_gt, ~] = impixel();
end
%%
xy_gt = [x_gt y_gt];
o_man = o_gt;
theta_man = zeros(8, K);
for i_k = 1:K

    xy_lr = [corr_points_x(:,i_k) corr_points_y(:,i_k)];

    tform = fitgeotrans(xy_lr,xy_gt,'projective');
    [x_lrt y_lrt] = transformPointsForward(tform, xy_lr(:,1), xy_lr(:,2));
    
    o_man(i_k).H = tform.T';
    theta_man(:, i_k) = tform.T(1:8);
    
    figure;
    subplot(1,2,1); imgray(o_gt(i_k).im);
    plot(xy_lr(:,1), xy_lr(:,2), 'rx');
    
    subplot(1,2,2); imgray(gtruth);
    plot(xy_gt(:,1), xy_gt(:,2), 'rx');
    plot(x_lrt, y_lrt, 'go');
    
end
%%
[Y, theta, lambda, gamma, ni_lo, nj_lo] = vectorise_lo_res_data(o_gt(1:K));
[avg_im_gt] = getAvim(ni_hi, nj_hi, ni_lo, nj_lo, Y, theta, gamma);
[W_gt] = makeW(ni_hi, nj_hi, ni_lo, nj_lo, theta, gamma);


KM = sum(ni_lo.*nj_lo);
lambda1 = ones(KM,1);
lambda2 = ones(KM,1);
curr_pt = 0;
for i_k = 1:K
    n_pts = ni_lo(i_k)*nj_lo(i_k);
    idx = curr_pt + (1:n_pts);
    curr_pt = curr_pt + n_pts;
    lambda1(idx) = lambda(1,i_k);
    lambda2(idx) = lambda(2,i_k);
end
%
alpha = 0.08;
nu = 0.04;

huber_im_gt = superres_huber(W_gt, Y, lambda1, lambda2, avg_im_gt, alpha, nu);
ml_im_gt = superres_ml(W_gt, Y, lambda1, lambda2, avg_im_gt);%

figure; 
subplot(1,3,1); imgray(avg_im_gt);
subplot(1,3,2); imgray(huber_im_gt);
subplot(1,3,3); imgray(ml_im_gt);
%%
[W_reg] = makeW(ni_hi, nj_hi, ni_lo, nj_lo, theta_man, gamma);
[avg_im_man] = getAvim(ni_hi, nj_hi, ni_lo, nj_lo, Y, theta_man, gamma);    
huber_im_man = superres_huber(W_reg, Y, lambda1, lambda2, avg_im_man, alpha, nu, 3);
ml_im_man = superres_ml(W_reg, Y, lambda1, lambda2, avg_im_man, 3);%

figure; 
subplot(1,3,1); imgray(avg_im_man);
subplot(1,3,2); imgray(huber_im_man); caxis([-.5 .5]);
subplot(1,3,3); imgray(ml_im_man); caxis([-.5 .5]);
%%
frames_h5 = h5read('C:\isbe\nailfold\camera_capture\general_imaging\test001\2017_03_22\L4_10_25_32\frames.h5', '/frames/data');
frames_r = permute(frames_h5(3:3:end,:,:), [2 1 4 3]);
frames_g = permute(frames_h5(2:3:end,:,:), [2 1 4 3]);
frames_b = permute(frames_h5(1:3:end,:,:), [2 1 4 3]);
frames_rgb = cat(3, frames_r, frames_g, frames_b);
%%
for i_f = (1:20) + 120
    figure;
    imgray(rgb2gray(frames_rgb(:,:,:,i_f)));
end
%%
K = 9;
frames = zeros(480, 640, K);

n_pts = 10;
corr_points_x = zeros(n_pts, K);
corr_points_y = zeros(n_pts, K);
figure;

for i_k = 1:K
    frames(:,:,i_k) = rot90(rgb2gray(frames_rgb(:,:,:,i_k + 124)),2);

    clf;
    imgray(frames(:,:,i_k));
    if (i_k > 1)
        plot(x, y, 'bo');
    end
    
    title(i_k);
    
    x = [];
    y = [];
    while length(x) ~= n_pts
        [x, y, ~] = impixel();
    end
    corr_points_x(:,i_k) = x(:);
    corr_points_y(:,i_k) = y(:);
end
figure; axis image xy; plot(corr_points_x', corr_points_y', '.');
%%
zm = 1.5;
xy_gt = zm * [mean(corr_points_x,2) mean(corr_points_y,2)];

theta_man = zeros(8, K);
for i_k = 1:K

    xy_lr = [corr_points_x(:,i_k) corr_points_y(:,i_k)];

    tform = fitgeotrans(xy_lr,xy_gt,'projective');
    [x_lrt y_lrt] = transformPointsForward(tform, xy_lr(:,1), xy_lr(:,2));
    
    theta_man(:, i_k) = tform.T(1:8);
   
end
%%
[ni_lo nj_lo, K] = size(frames);
ni_hi = zm*ni_lo;
nj_hi = zm*nj_lo;
gamma = (0.4^2)*ones(1, K);
[avg_im_man] = getAvim(ni_hi, nj_hi, ni_lo, nj_lo, double(frames(:)), theta_man, gamma);
figure; imgray(avg_im_man);

imwrite(uint8(avg_im_man), 'C:\isbe\nailfold\camera_capture\general_imaging\test001\2017_03_21\L4_15_21_19\sr_mosaic.png');
%%
frames_small = frames(1:200, 1:300, :);
[ni_lo nj_lo, K] = size(frames_small);
ni_lo = ni_lo*ones(K,1);
nj_lo = nj_lo*ones(K,1);
ni_hi = zm*ni_lo;
nj_hi = zm*nj_lo;
gamma = (0.4^2)*ones(1, K);

[avg_im_man] = getAvim(ni_hi, nj_hi, ni_lo, nj_lo, double(frames_small(:))/255 - 0.5, theta_man, gamma);

%%
KM = sum(ni_lo.*nj_lo);
lambda1 = ones(KM,1);
lambda2 = zeros(KM,1);
%%
alpha = 0.08;
nu = 0.04;
[W_man] = makeW(ni_hi, nj_hi, ni_lo, nj_lo, theta_man, gamma);
huber_im_man = superres_huber(W_man, double(frames_small(:))/255 - 0.5, lambda1, lambda2, avg_im_man, alpha, nu, 3);
ml_im_man = superres_ml(W_man, double(frames_small(:))/255 - 0.5, lambda1, lambda2, avg_im_man, 3);%

figure; 
subplot(1,3,1); imgray(avg_im_man);
subplot(1,3,2); imgray(huber_im_man); caxis([-.5 .5]);
subplot(1,3,3); imgray(ml_im_man); caxis([-.5 .5]);
%%
frames_dir = 'C:\isbe\nailfold\camera_capture\general_imaging\test001\2017_03_21\L4_15_21_19\';
caps_sr = load([frames_dir 'detected_capillaries_cxx/sr_mosaicapex_candidates.txt']);
caps_lo = load([frames_dir 'detected_capillaries_cxx/raw_frameapex_candidates.txt']);

avg_im_man = imread([frames_dir 'sr_mosaic.png']);
raw_frame = imread([frames_dir 'raw_frame.png']);
%%
distal_sr = caps_sr(:,6) > 0.6;
distal_lo = caps_lo(:,6) > 0.6;

figure;
subplot(1,2,1);
imgray(avg_im_man);
for i_ap = 1:50
    text(caps_sr(i_ap,1)/0.622, caps_sr(i_ap,2)/0.622, num2str(i_ap), 'color', 'r')
end
plot(caps_sr(distal_sr,1)/0.622, caps_sr(distal_sr,2)/0.622, 'g.');

subplot(1,2,2);
imgray(raw_frame);
for i_ap = 1:50
    text(caps_lo(i_ap,1)/0.933, caps_lo(i_ap,2)/0.933, num2str(i_ap), 'color', 'r')
end
plot(caps_sr(distal_sr,1)/(1.5*0.622)+25, caps_sr(distal_sr,2)/(1.5*0.622)+4, 'g.');
%
figure;
subplot(1,2,1);
imgray(avg_im_man);
plot(caps_sr(distal_sr,1)/0.622, caps_sr(distal_sr,2)/0.622, 'g.');

subplot(1,2,2);
imgray(raw_frame);
plot(caps_lo(1:8,1)/(0.933), caps_lo(1:8,2)/(0.933), 'g.');
plot(caps_lo(9:21,1)/(0.933), caps_lo(9:21,2)/(0.933), 'r.');
%%
figure;
imgray(avg_im_man);
plot(caps_sr(distal_sr,1)/0.622, caps_sr(distal_sr,2)/0.622, 'g^', 'markerfacecolor', 'g', 'markersize', 8);

figure;
imgray(raw_frame);
plot(caps_lo(1:8,1)/(0.933), caps_lo(1:8,2)/(0.933), 'g^', 'markerfacecolor', 'g', 'markersize', 8);
plot(caps_lo([9 10 11 14 16 17 21],1)/(0.933), caps_lo([9 10 11 14 16 17 21],2)/(0.933), 'co', 'markerfacecolor', 'c', 'markersize', 8);
plot(caps_lo([12 13 15 18 19 20],1)/(0.933), caps_lo([12 13 15 18 19 20],2)/(0.933), 'ro', 'markerfacecolor', 'r', 'markersize', 8);
%%
frames_dir = 'C:\isbe\nailfold\camera_capture\cic_2016\004cic_2016\2016_12_06\L4_12_25_04\';
K = 11;

for i_fr = 1:K
    frame = imresize(rgb2gray(rot90(imread([frames_dir 'frame' zerostr(i_fr + 177, 5) '.bmp']),2)), 0.5);
    if i_fr == 1
        frames = zeros(size(frame,1), size(frame,2), K);
    end
    frames(:,:,i_fr) = double(frame);
    figure; imgray(frames(:,:,i_fr));
end
%%
n_pts = 10;
corr_points_x = zeros(n_pts, K);
corr_points_y = zeros(n_pts, K);

figure;
for i_k = 1:K
    clf;
    imgray(frames(:,:,i_k));
    if (i_k > 1)
        for i_pt = 1:n_pts
        	text(x(i_pt), y(i_pt), num2str(i_pt), 'color', 'r');
        end
    end
    
    title(i_k);
    
    x = [];
    y = [];
    while length(x) ~= n_pts
        [x, y, ~] = impixel();
    end
    corr_points_x(:,i_k) = x(:);
    corr_points_y(:,i_k) = y(:);
end
figure; axis image xy; plot(corr_points_x', corr_points_y', '.');
%%
zm = 1.0;
xy_gt = zm * [mean(corr_points_x,2) mean(corr_points_y,2)];

theta_man = zeros(8, K);
for i_k = 1:K

    xy_lr = [corr_points_x(:,i_k) corr_points_y(:,i_k)];

    tform = fitgeotrans(xy_lr,xy_gt,'projective');
    [x_lrt y_lrt] = transformPointsForward(tform, xy_lr(:,1), xy_lr(:,2));
    
    theta_man(:, i_k) = tform.T(1:8);
   
end
%%
[ni_lo nj_lo, K] = size(frames);
ni_hi = zm*ni_lo;
nj_hi = zm*nj_lo;
gamma = (0.4^2)*ones(1, K);
[avg_im_man] = getAvim(ni_hi, nj_hi, ni_lo, nj_lo, frames(:), theta_man, gamma);
mask = avg_im_man > 0;

figure; imgray(avg_im_man); caxis([min(avg_im_man(mask)) max(avg_im_man(mask))]);

imwrite(uint8(avg_im_man), [frames_dir 'sr_mosaic2.png']);
imwrite(uint8(avg_im_man(1:400,:)), [frames_dir 'sr_mosaic2_cropped.png']);

caps_usb = load([frames_dir 'detected_capillaries_cxx\sr_mosaic2_croppedapex_candidates.txt']);
%%   
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
frames_dir = 'C:\isbe\nailfold\camera_capture\general_imaging\test001\2017_03_22\L4_11_59_50\';
[frames_rgb] = load_hdf5_frames([frames_dir 'frames.h5'], [], [], 2, [frames_dir 'png\']);

%%
K = 15;
frames = zeros(480, 640, K);
for i_fr = 1:K
    frames(:,:,i_fr) = double(rgb2gray(frames_rgb(:, :, :, i_fr + 15)));    
    figure; imgray(frames(:,:,i_fr));
end
%%
n_pts = 10;
corr_points_x = zeros(n_pts, K);
corr_points_y = zeros(n_pts, K);
%%
figure;
for i_k = 11:K
    clf;
    imgray(frames(:,:,i_k));
    if (i_k > 1)
        for i_pt = 1:n_pts
        	text(x(i_pt), y(i_pt), num2str(i_pt), 'color', 'r');
        end
    end
    
    title(i_k);
    
    x = [];
    y = [];
    while length(x) ~= n_pts
        [x, y, ~] = impixel();
    end
    corr_points_x(:,i_k) = x(:);
    corr_points_y(:,i_k) = y(:);
end
figure; axis image xy; plot(corr_points_x', corr_points_y', '.');
%%
zm = 1.5;
xy_gt = zm * [mean(corr_points_x,2) mean(corr_points_y,2)];

theta_man = zeros(8, K);
for i_k = 1:K

    xy_lr = [corr_points_x(:,i_k) corr_points_y(:,i_k)];

    tform = fitgeotrans(xy_lr,xy_gt,'projective');
    [x_lrt y_lrt] = transformPointsForward(tform, xy_lr(:,1), xy_lr(:,2));
    
    theta_man(:, i_k) = tform.T(1:8);
   
end
%%
[ni_lo nj_lo, K] = size(frames);
ni_hi = zm*ni_lo;
nj_hi = zm*nj_lo;
gamma = (0.4^2)*ones(1, K);
[avg_im_man] = getAvim(ni_hi, nj_hi, ni_lo, nj_lo, frames(:), theta_man, gamma);
mask = avg_im_man > 0;

figure; imgray(avg_im_man); caxis([min(avg_im_man(mask)) max(avg_im_man(mask))]);

imwrite(uint8(avg_im_man), [frames_dir 'sr_mosaic.png']);
%%
frames_dir = 'C:\isbe\nailfold\camera_capture\cic_2016\004cic_2016\2016_12_06\L4_12_27_34\';
K = 15;

for i_fr = 1:K
    frame = rgb2gray(rot90(imread([frames_dir 'frame' zerostr(i_fr + 65, 5) '.bmp']),2));
    if i_fr == 1
        frames = zeros(size(frame,1), size(frame,2), K);
    end
    frames(:,:,i_fr) = double(frame);
    figure; imgray(frames(:,:,i_fr));
end
%%
n_pts = 10;
corr_points_x = zeros(n_pts, K);
corr_points_y = zeros(n_pts, K);

figure;
for i_k = 1:K
    clf;
    imgray(frames(:,:,i_k));
    if (i_k > 1)
        for i_pt = 1:n_pts
        	text(x(i_pt), y(i_pt), num2str(i_pt), 'color', 'r');
        end
    end
    
    title(i_k);
    
    x = [];
    y = [];
    while length(x) ~= n_pts
        [x, y, ~] = impixel();
    end
    corr_points_x(:,i_k) = x(:);
    corr_points_y(:,i_k) = y(:);
end
figure; axis image xy; plot(corr_points_x', corr_points_y', '.');
%%
zm = 4;
xy_gt = zm * [mean(corr_points_x,2) mean(corr_points_y,2)];

theta_man = zeros(8, K);
for i_k = 1:K

    xy_lr = [corr_points_x(:,i_k) corr_points_y(:,i_k)];

    tform = fitgeotrans(xy_lr,xy_gt,'projective');
    [x_lrt y_lrt] = transformPointsForward(tform, xy_lr(:,1), xy_lr(:,2));
    
    theta_man(:, i_k) = tform.T(1:8);
   
end
%%
[ni_lo nj_lo, K] = size(frames);
ni_hi = zm*ni_lo;
nj_hi = zm*nj_lo;
gamma = (0.4^2)*ones(1, K);
[avg_im_man] = getAvim(ni_hi, nj_hi, ni_lo, nj_lo, frames(:), theta_man, gamma);
[avg_im_mask] = getAvim(ni_hi, nj_hi, ni_lo, nj_lo, frames(:), theta_man, gamma, frames(:) > 175);
mask = avg_im_man > 0;

figure; imgray(avg_im_man); caxis([min(avg_im_man(mask)) max(avg_im_man(mask))]);
%%
imwrite(uint8(avg_im_man), [frames_dir 'sr_mosaic.png']);
imwrite(uint8(avg_im_man(1:800,:)), [frames_dir 'sr_mosaic_cropped.png']);
%%
frames_dir = 'C:\isbe\nailfold\camera_capture\cic_2016\013cic_2016\2016_12_13\L4_11_21_26\';
frames_idx = 28:42;
[sr_im, frames] = usb_manual_super_res(frames_dir, frames_idx);
%%
frames_dir = 'C:\isbe\nailfold\camera_capture\cic_2016\004cic_2016\2016_12_06\L4_12_27_34\';
frames_idx = 66:80;
[sr_im, frames] = usb_manual_super_res(frames_dir, frames_idx);
imwrite(uint8(sr_im(1:800,:)), fullfile([frames_dir 'sr_pts'], 'sr_mosaic_cropped.png'));

frame_enlarged = imresize(frames(:,:,15), 4, 'lanczos3');

imwrite(uint8(sr_im(1:800,:)), fullfile([frames_dir 'sr_pts'], 'sr_mosaic_cropped.png'));
imwrite(uint8(frame_enlarged(1:800,:)), fullfile([frames_dir 'sr_pts'], 'frame_mosaic_cropped.png'));
%%
mosaic_usb = imread([frames_dir 'sr_pts\sr_mosaic_cropped.png']);
caps_usb = load([frames_dir 'sr_pts\detected_capillaries_cxx\sr_mosaic_croppedapex_candidates.txt']);
distal_sr = (caps_usb(:,11) > 0 | caps_usb(:,6) > 0.8) & caps_usb(:,1)/0.622 < 2000;
nondistal_sr = ~distal_sr & caps_usb(:,6) > 0.6 & caps_usb(:,8) > 0 & caps_usb(:,1)/0.622 < 2000;
figure(...
    'windowstyle', 'normal',...
    'units', 'centimeters',...
    'position', [2 2 30, 10]); imgray(mosaic_usb(201:end, 1:2000));
plot(caps_usb(nondistal_sr,1) / 0.622, caps_usb(nondistal_sr,2) / 0.622 - 200, 'go', 'markerfacecolor', 'g', 'markersize', 4);
plot(caps_usb(distal_sr,1) / 0.622, caps_usb(distal_sr,2) / 0.622 - 200, 'go', 'markerfacecolor', 'g', 'markersize', 4);
axis off;
exportfig([fig_dir 'mosaic_usb_with_overlay']);
%%
fig_dir = 'C:\Users\momeemb2\Dropbox (The University of Manchester)\nailfold\funding\mrc_dpfs\MArch 2017\figs\';

frame_usb = imread([frames_dir 'sr_pts\frame_mosaic_cropped.png']);
caps_usb = load([frames_dir 'sr_pts\detected_capillaries_cxx\frame_mosaic_croppedapex_candidates.txt']);
distal_sr = caps_usb(:,11) > 0 | caps_usb(:,6) > 0.8;
nondistal_sr = ~distal_sr & caps_usb(:,6) > 0.6 & caps_usb(:,8) > 0;
figure('windowstyle', 'normal'); imgray(frame_usb);
plot(caps_usb(nondistal_sr,1) / 0.622, caps_usb(nondistal_sr,2) / 0.622, 'mo', 'markerfacecolor', 'm', 'markersize', 4);
plot(caps_usb(distal_sr,1) / 0.622, caps_usb(distal_sr,2) / 0.622, 'go', 'markerfacecolor', 'g', 'markersize', 8);



%%
frames_dir = 'C:\isbe\nailfold\camera_capture\cic_2016\004cic_2016\2016_12_06\L4_12_25_04\';
frames_idx = 178:188;
%%
[sr_im] = optilia_manual_super_res(frames_dir, frames_idx);
imwrite(uint8(sr_im(1:400,:)), fullfile([frames_dir 'sr_pts'], 'sr_mosaic_cropped.png'));
%%
mosaic_opt = imread([frames_dir 'sr_pts\mosaic_opt.png']);
caps_opt = load([frames_dir 'sr_pts\detected_capillaries_cxx\mosaic_optapex_candidates.txt']);
distal_sr = (caps_opt(:,11) > 0 | caps_opt(:,6) > 0.8);
nondistal_sr = ~distal_sr & caps_opt(:,6) > 0.6 & caps_opt(:,8) > 0;

%roi_opt = mosaic_opt(401:1200,:,1);
%write_im_from_colormap(roi_opt, [frames_dir 'sr_pts\frame00165_roi.png'], gray(256));

figure(...
    'windowstyle', 'normal',...
    'units', 'centimeters',...
    'position', [2 2 30, 10]);

imgray(mosaic_opt);
plot(caps_opt(nondistal_sr,1), caps_opt(nondistal_sr,2), 'go', 'markerfacecolor', 'g', 'markersize', 4);
plot(caps_opt(distal_sr,1), caps_opt(distal_sr,2), 'go', 'markerfacecolor', 'g', 'markersize', 4);
axis off;
exportfig([fig_dir 'mosaic_opt_with_overlay']);
%%
frames_dir = 'C:\isbe\nailfold\camera_capture\cic_2016\004cic_2016\2016_12_06\L4_12_18_20\';
mosaic_wel = imread([frames_dir 'sequence_data\full_mosaic.png']);
bg_mask = ~mosaic_wel;
mosaic_wel(bg_mask) = max(mosaic_wel(~bg_mask));

caps_wel = load([frames_dir 'sequence_data\detected_capillaries_cxx\full_mosaicapex_candidates.txt']);
distal_sr = (caps_wel(:,11) > 0 | caps_wel(:,6) > 0.8) & caps_wel(:,1)/0.622 < 4000 & caps_wel(:,1)/0.622 > 1000;
nondistal_sr = ~distal_sr & caps_wel(:,6) > 0.6 & caps_wel(:,8) > 0 & caps_wel(:,1)/0.622 < 4000 & caps_wel(:,1)/0.622 > 1000;
figure(...
    'windowstyle', 'normal',...
    'units', 'centimeters',...
    'position', [2 2 30, 10]); imgray(mosaic_wel(701:end, 1001:4000));
axis off;

plot(caps_wel(nondistal_sr,1) / 0.622 - 1000, caps_wel(nondistal_sr,2) / 0.622 - 700, 'go', 'markerfacecolor', 'g', 'markersize', 4);
plot(caps_wel(distal_sr,1) / 0.622 - 1000, caps_wel(distal_sr,2) / 0.622 - 700, 'go', 'markerfacecolor', 'g', 'markersize', 4);
exportfig([fig_dir 'mosaic_wel_with_overlay']);
%%

frames_dir = 'C:\isbe\nailfold\camera_capture\cic_2016\004cic_2016\2016_12_06\L4_12_27_34\';
load(fullfile([frames_dir 'sr_pts'], 'sr_pts.mat'), 'corr_points*');
%%
i_pt = 9;
for i_k = 1:K; 
    roi = frames(corr_points_y(i_pt,i_k) + (-31:64), corr_points_x(i_pt,i_k) + (-31:40), i_k);
    figure; imgray(roi); 
    write_im_from_colormap(kron(roi, ones(4)), ...
        fullfile([frames_dir 'sr_pts'], ['frame_roi' zerostr(i_k,2) '.png']), gray(256));
end

roi_usb = mosaic_usb(453 + (-127:256), 1577 + (-127:160));
figure; imgray(roi_usb);

write_im_from_colormap(roi_usb, ...
    fullfile([frames_dir 'sr_pts'], 'sr_roi.png'), gray(256));
%%    
for i_k = [1 7 15]; 
    frame_enlarged = uint8(imresize(frames(:,:,i_k), 4, 'lanczos3'));
    imwrite(frame_enlarged, ...
        fullfile([frames_dir 'sr_pts'], ['frame' zerostr(i_k,2) '.png']));
end
%%
frame = rot90(rgb2gray(imread([frames_dir 'frame00165.bmp'])),2);
frame2 = rot90(rgb2gray(imread([frames_dir 'frame00188.bmp'])),2);
frames_opt = double(cat(3, imresize(frame, 0.5, 'lanczos3'), imresize(frame2, 0.5, 'lanczos3')));

[init_reg, match_counts] = ...
        register_tiles_features(frames_opt, ...
                            'ref_type', 'previous',...
                            'region', 'centre',...
                            'theta_range', -5:0.5:5, ...
                            'offset_lim', 10, ...
                            'mosaic', [],...
                            'sigma', 6,...
                            'tile_masks', [],...
                            'offset_centres', [548 222; 0 0],...
                            'debug', 1);
                        
%init_reg = cat(3, eye(3), [1 0 548; 0 1 222; 0 0 1]);                       
                        
[init_avg, ~, t2m_reg ] = create_mosaic(frames_opt, init_reg);
figure; imgray(init_avg(231:750,1:end-50));
imwrite(uint8(init_avg(231:750,1:end-50)), [frames_dir '/sr_pts/mosaic_opt.png']);


    




    


    



