%--------------------------------------------------------------------------
% Script for producing results for journal submission
% Experiment predicting orientation of synthetic lines in the presence of increasing noise 

%--------------------------------------------------------------------------
%%
%0. Generic arguments to initialise
clear
data_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\increasing_noise_exp\';
exp_dir = 'C:\isbe\asymmetry_project\experiments\synthetic_lines\increasing_noise_exp\';

sz_x = 64;
sz_y = 64;
s_noise = [0.25 0.50 1.00 2.00];
warning('off', 'ASYM:unexpectedArgument');
%%
%1. Train a set of random forests on line images with varying degrees of
%rician and gaussian noise

mkdir([data_dir 'bg']);
bg = ones(sz_y, sz_x);
save([data_dir 'bg\bg00001.mat'], 'bg');
%  
% Commands to commit jobs to the CSF to train
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[0.25]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[0.50]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="gaussian" NOISE_PARAMS="[0.25]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="gaussian" NOISE_PARAMS="[0.50]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="gaussian" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="gaussian" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="g2d" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh

NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="dt" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
NUM_TREES=10 OUTPUT_TYPE="centre_orientation" IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[2.00]" MODEL_ROOT="models/synthetic_lines" PREDICTION_TYPE="rf_regression" DECOMP_TYPE="dt" WIN_SIZE=1 qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh

MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5671'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5672'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5763'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5764'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5675'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5676'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5677'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5678'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'5679'" MAKE_SAMPLED_MAPS=0 qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh

MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'10210'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10210 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_orientation/rf_regression" MODEL_PATH="'10211'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10211 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh

%
%Get the forests from the CSF and change their tree root
for ii = [5671:5679 10210:10211]
    predictor = u_load([data_dir 'rfs\' num2str(ii) '\predictor.mat']);
    predictor.tree_root = [data_dir 'rfs\'];
    save([data_dir 'rfs\' num2str(ii) '\predictor.mat'], 'predictor');
end
    
    
%%
%2. Generate lots of test images formed using the same model, and compare
%forest predictions with the analytic method

%Define params for testing
sigma = [1 2 4 8];
pts_per_im = 10;
num_ims = 2000;

%Load in the predictors
signal_predictor = u_load([data_dir 'rfs\5671\predictor.mat']);
rician_predictor{1} = u_load([data_dir 'rfs\5672\predictor.mat']);
rician_predictor{2} = u_load([data_dir 'rfs\5673\predictor.mat']);
rician_predictor{3} = u_load([data_dir 'rfs\5674\predictor.mat']);
rician_predictor{4} = u_load([data_dir 'rfs\5675\predictor.mat']);
gaussian_predictor{1} = u_load([data_dir 'rfs\5676\predictor.mat']);
gaussian_predictor{2} = u_load([data_dir 'rfs\5677\predictor.mat']);
gaussian_predictor{3} = u_load([data_dir 'rfs\5678\predictor.mat']);
gaussian_predictor{4} = u_load([data_dir 'rfs\5679\predictor.mat']);

D_g2 = signal_predictor.D;
job_args = u_load([data_dir 'rfs\5671\job_args.mat']);
d_args_g2 = job_args.decomposition_args;

dt_predictor = u_load([data_dir 'rfs\10210\predictor.mat']);
D_dt = dt_predictor.D;
job_args = u_load([data_dir 'rfs\10210\job_args.mat']);
d_args_dt = job_args.decomposition_args;

%Pre-allocate space for the true line parameters
true_oris = zeros(num_ims*pts_per_im,1);
true_cons = zeros(num_ims*pts_per_im,1);
true_widths = zeros(num_ims*pts_per_im,1);

%Pre-allocate space for analytic errors
ori_errors_signal = zeros(num_ims*pts_per_im,1);
ori_errors_gaussian = zeros(num_ims*pts_per_im,4);
ori_errors_rician = zeros(num_ims*pts_per_im,4);

%Pre-allocate psace for the responses to input to the predictors
signal_data_g2 = zeros(num_ims*pts_per_im, D_g2, 1);
gaussian_data_g2 = zeros(num_ims*pts_per_im, D_g2, 4);
rician_data_g2 = zeros(num_ims*pts_per_im, D_g2, 4);

signal_data_dt = zeros(num_ims*pts_per_im, D_dt, 1);
gaussian_data_dt = zeros(num_ims*pts_per_im, D_dt, 4);
rician_data_dt = zeros(num_ims*pts_per_im, D_dt, 4);

for ii = 1:num_ims
    
    display(['Testing image ' num2str(ii)]);
    
    %Sample properties of line
    line_width = sample_uniform([1 8]);
    line_contrast = sample_uniform([1 8]);
    line_ori = sample_uniform([0 180]);
    line_rad = pi * line_ori / 180;
    
    %Generate line
    [line, label, label_centre] =...
        create_ellipse_bar(line_width/2, line_contrast, line_ori, sz_y, sz_x, sz_x/2, sz_y/2);
    signal = 1 + line;
    
    %Select some random pixels in the signal image
    shuffle = randperm(sum(label_centre(:)));
    idx = find(label_centre);
    r_idx = idx(shuffle(1:pts_per_im));
    [r_rows r_cols] = ind2sub([sz_y sz_x], r_idx);
    sample_idx = ((ii-1)*pts_per_im + 1):(ii*pts_per_im);
    
    %Save the line parameters
    true_oris(sample_idx,:) = line_rad;
    true_cons(sample_idx,:) = line_contrast;
    true_widths(sample_idx,:) = line_width;
    
    %Predict orientation in signal image using analytic G" method    
    [~, ori_map] = gaussian_2nd_derivative_line(signal, sigma);
    ori_errors_signal(sample_idx,:) =...
        ori_error(ori_map(r_idx), line_rad);

    %Compute signal responses for RF
    signal_responses = compute_filter_responses(signal, d_args_g2);
    signal_data_g2(sample_idx, :, 1) =...
        sample_image_features(signal_responses, r_rows(:), r_cols(:), d_args_g2);
    
    signal_responses = compute_filter_responses(signal, d_args_dt);
    signal_data_dt(sample_idx, :, 1) =...
        sample_image_features(signal_responses, r_rows(:), r_cols(:), d_args_dt);

    for ss = 1:4
        
        %Corrupt the signal with Rician noise and recompute G" errors
        test_image = ricernd(signal, s_noise(ss));
        [~, ori_map] = gaussian_2nd_derivative_line(test_image, sigma);
        ori_errors_rician(sample_idx,ss) = ...
            ori_error(ori_map(r_idx), line_rad);
        
        %Compute filter responses for RF
        test_responses = compute_filter_responses(test_image, d_args_g2);
        rician_data_g2(sample_idx, :, ss) =...
            sample_image_features(test_responses, r_rows(:), r_cols(:), d_args_g2);
        test_responses = compute_filter_responses(test_image, d_args_dt);
        rician_data_dt(sample_idx, :, ss) =...
            sample_image_features(test_responses, r_rows(:), r_cols(:), d_args_dt);
    
        %Corrupt signal with Gaussian noise
        test_image = signal + reshape(sample_from_normal(0, s_noise(ss)^2, sz_x*sz_y), sz_y, sz_x);
        [~, ori_map] = gaussian_2nd_derivative_line(test_image, sigma);
        ori_errors_gaussian(sample_idx,ss) =...
            ori_error(ori_map(r_idx), line_rad);
        
        %Compute filter responses for RF
        test_responses = compute_filter_responses(test_image, d_args_g2);
        gaussian_data_g2(sample_idx, :, ss) =...
            sample_image_features(test_responses, r_rows(:), r_cols(:), d_args_g2);
        test_responses = compute_filter_responses(test_image, d_args_dt);
        gaussian_data_dt(sample_idx, :, ss) =...
            sample_image_features(test_responses, r_rows(:), r_cols(:), d_args_dt);
        
    end
end
%Save the workspace - it'll be big, but may be useful later
save([exp_dir 'workspace.mat']);
%%
load([exp_dir 'workspace.mat']);
%
%Predict orientations for each rf and compute the errors

%signal
ori_rf = random_forest_reg_predict(signal_predictor, signal_data_g2, 1);
ori_errors_signal_g2 = ori_error(ori_rf, true_oris);

ori_errors_gaussian_g2 = zeros(num_ims*pts_per_im,4);
ori_errors_rician_g2 = zeros(num_ims*pts_per_im,4);

for ss = 1:4
    %Rician
    ori_rf = random_forest_reg_predict(rician_predictor{ss}, rician_data_g2(:,:,ss), 1);
    ori_errors_rician_g2(:,ss) = ori_error(ori_rf, true_oris);
    %Gaussian
    ori_rf = random_forest_reg_predict(gaussian_predictor{ss}, gaussian_data_g2(:,:,ss), 1);
    ori_errors_gaussian_g2(:,ss) = ori_error(ori_rf, true_oris);
end

%Save all the results
save([exp_dir 'ori_errors.mat'], 'ori_errors_*');
save([exp_dir 'line_params.mat'], 'true_*');
%%
% 3. Analyse the results

%Load in the results and compute and plot PDFs for the various predictions
load([exp_dir 'ori_errors.mat'], 'ori_errors_*');
load([exp_dir 'line_params.mat'], 'true_*');

sorted_errors_signal = sort(abs(ori_errors_signal)) * 180 / pi;
sorted_errors_signal_g2 = sort(abs(ori_errors_signal_g2)) * 180 / pi;
sorted_errors_gaussian = sort(abs(ori_errors_gaussian)) * 180 / pi;
sorted_errors_gaussian_g2 = sort(abs(ori_errors_gaussian_g2)) * 180 / pi;
sorted_errors_rician = sort(abs(ori_errors_rician)) * 180 / pi;
sorted_errors_rician_g2 = sort(abs(ori_errors_rician_g2)) * 180 / pi;

xi = [1 100:100:20000];
xp = xi / 20000;
%
figure; hold all; title('Noise-free errors CDF');
plot([sorted_errors_signal(xi) sorted_errors_signal_g2(xi)], xp);
axis([0 45 0 1]);
exportfig([exp_dir 'noise_free_errors_an_v_rf.pdf']);

for ss = 1:4
    figure; hold all; title(['Gaussian errors (' num2str(s_noise(ss)) ') CDF']);
    plot([sorted_errors_gaussian(xi,ss) sorted_errors_gaussian_g2(xi,ss)], xp);
    axis([0 45 0 1]);
    exportfig([exp_dir 'gaussian_errors_an_v_rf' num2str(ss) '.pdf']);
    
    figure; hold all; title(['Rician errors (' num2str(s_noise(ss)) ') CDF']);
    plot([sorted_errors_rician(xi,ss) sorted_errors_rician_g2(xi,ss)], xp);
    axis([0 45 0 1]);
    exportfig([exp_dir 'rician_errors_an_v_rf' num2str(ss) '.pdf']);
end

figure; hold all; plot(sorted_errors_gaussian(xi,:), xp); axis([0 45 0 1]);
title('Errors for analytic prediction - Gaussian noise');
legend({'\sigma = 0.25', '\sigma = 0.50', '\sigma = 1.00', '\sigma = 2.00'});
exportfig([exp_dir 'gaussian_errors_an.pdf']);

figure; hold all; plot(sorted_errors_gaussian_g2(xi,:), xp); axis([0 45 0 1]);
title('Errors for RF prediction - Gaussian noise');
legend({'\sigma = 0.25', '\sigma = 0.50', '\sigma = 1.00', '\sigma = 2.00'});
exportfig([exp_dir 'gaussian_errors_rf.pdf']);

figure; hold all; plot(sorted_errors_rician(xi,:), xp); axis([0 45 0 1]);
title('Errors for analytic prediction - Rician noise');
legend({'SNR = 0.25', 'SNR = 0.50', 'SNR = 1.00', 'SNR = 2.00'});
exportfig([exp_dir 'rician_errors_an.pdf']);

figure; hold all; plot(sorted_errors_rician_g2(xi,:), xp); axis([0 45 0 1]);
title('Errors for RF prediction - Rician noise');
legend({'SNR = 0.25', 'SNR = 0.50', 'SNR = 1.00', 'SNR = 2.00'});
exportfig([exp_dir 'rician_errors_rf.pdf']);
%%
%Compute kernel smoothed errors as a function of line contrast, width and
%orientation
smoother_width = 500;
smoother_pts = 1000;

[widths_x, widths_y] = moving_average_smoother(true_widths, abs(ori_errors_signal)*180/pi, smoother_width, smoother_pts);
[cons_x, cons_y] = moving_average_smoother(true_cons, abs(ori_errors_signal)*180/pi, smoother_width, smoother_pts);
[oris_x, oris_y] = moving_average_smoother(true_oris, abs(ori_errors_signal)*180/pi, smoother_width, smoother_pts);
[widths_xrf, widths_yrf] = moving_average_smoother(true_widths, abs(ori_errors_signal_g2)*180/pi, smoother_width, smoother_pts);
[cons_xrf, cons_yrf] = moving_average_smoother(true_cons, abs(ori_errors_signal_g2)*180/pi, smoother_width, smoother_pts);
[oris_xrf, oris_yrf] = moving_average_smoother(true_oris, abs(ori_errors_signal_g2)*180/pi, smoother_width, smoother_pts);

fw1g = figure; hold all; axis([1 8 0 30]); title('Error as a fucntion of line width, analytic prediction, Gaussian noise');
plot(widths_x, widths_y);
fw2g = figure; hold all; axis([1 8 0 30]); title('Error as a fucntion of line width, RF prediction, Gaussian noise');
plot(widths_xrf, widths_yrf);

fc1g = figure; hold all; axis([1 8 0 30]); title('Error as a fucntion of line contrast, analytic prediction, Gaussian noise');
plot(cons_x, cons_y);
fc2g = figure; hold all; axis([1 8 0 30]);title('Error as a fucntion of line contrast, RF prediction, Gaussian noise');
plot(cons_xrf, cons_yrf);

fo1g = figure; hold all; axis([0 pi 0 30]); title('Error as a fucntion of line orientation, analytic prediction, Gaussian noise');
plot(oris_x, oris_y);
fo2g = figure; hold all; axis([0 pi 0 30]); title('Error as a fucntion of line orientation, RF prediction, Gaussian noise');
plot(oris_xrf, oris_yrf);

fw1r = figure; hold all; axis([1 8 0 30]); title('Error as a fucntion of line width, analytic prediction, Rician noise');
plot(widths_x, widths_y);
fw2r = figure; hold all; axis([1 8 0 30]); title('Error as a fucntion of line width, RF prediction, Rician noise');
plot(widths_xrf, widths_yrf);

fc1r = figure; hold all; axis([1 8 0 30]); title('Error as a fucntion of line contrast, analytic prediction, Rician noise');
plot(cons_x, cons_y);
fc2r = figure; hold all; axis([1 8 0 30]);title('Error as a fucntion of line contrast, RF prediction, Rician noise');
plot(cons_xrf, cons_yrf);

fo1r = figure; hold all; axis([0 pi 0 30]); title('Error as a fucntion of line orientation, analytic prediction, Rician noise');
plot(oris_x, oris_y);
fo2r = figure; hold all; axis([0 pi 0 30]); title('Error as a fucntion of line orientation, RF prediction, Rician noise');
plot(oris_xrf, oris_yrf);
%
for ss = 1:4

    [widths_x, widths_y] = moving_average_smoother(true_widths, abs(ori_errors_gaussian(:,ss))*180/pi, smoother_width, smoother_pts);
    [cons_x, cons_y] = moving_average_smoother(true_cons, abs(ori_errors_gaussian(:,ss))*180/pi, smoother_width, smoother_pts);
    [oris_x, oris_y] = moving_average_smoother(true_oris, abs(ori_errors_gaussian(:,ss))*180/pi, smoother_width, smoother_pts);
    
    [widths_xrf, widths_yrf] = moving_average_smoother(true_widths, abs(ori_errors_gaussian_g2(:,ss))*180/pi, smoother_width, smoother_pts);
    [cons_xrf, cons_yrf] = moving_average_smoother(true_cons, abs(ori_errors_gaussian_g2(:,ss))*180/pi, smoother_width, smoother_pts);
    [oris_xrf, oris_yrf] = moving_average_smoother(true_oris, abs(ori_errors_gaussian_g2(:,ss))*180/pi, smoother_width, smoother_pts);
    
    figure(fw1g); plot(widths_x, widths_y);
    figure(fw2g); plot(widths_xrf, widths_yrf);
    
    figure(fc1g); plot(cons_x, cons_y);
    figure(fc2g); plot(cons_xrf, cons_yrf);
    
    figure(fo1g); plot(oris_x, oris_y);
    figure(fo2g); plot(oris_xrf, oris_yrf);
    
    [widths_x, widths_y] = moving_average_smoother(true_widths, abs(ori_errors_rician(:,ss))*180/pi, smoother_width, smoother_pts);
    [cons_x, cons_y] = moving_average_smoother(true_cons, abs(ori_errors_rician(:,ss))*180/pi, smoother_width, smoother_pts);
    [oris_x, oris_y] = moving_average_smoother(true_oris, abs(ori_errors_rician(:,ss))*180/pi, smoother_width, smoother_pts);
    
    [widths_xrf, widths_yrf] = moving_average_smoother(true_widths, abs(ori_errors_rician_g2(:,ss))*180/pi, smoother_width, smoother_pts);
    [cons_xrf, cons_yrf] = moving_average_smoother(true_cons, abs(ori_errors_rician_g2(:,ss))*180/pi, smoother_width, smoother_pts);
    [oris_xrf, oris_yrf] = moving_average_smoother(true_oris, abs(ori_errors_rician_g2(:,ss))*180/pi, smoother_width, smoother_pts);
    
    figure(fw1r); plot(widths_x, widths_y);
    figure(fw2r); plot(widths_xrf, widths_yrf);
    
    figure(fc1r); plot(cons_x, cons_y);
    figure(fc2r); plot(cons_xrf, cons_yrf);
    
    figure(fo1r); plot(oris_x, oris_y);
    figure(fo2r); plot(oris_xrf, oris_yrf);
end

figure(fw1g); exportfig([exp_dir 'loc_av_width_gaussian_an.pdf']);
figure(fw2g); exportfig([exp_dir 'loc_av_width_gaussian_rf.pdf']);
figure(fc1g); exportfig([exp_dir 'loc_av_con_gaussian_an.pdf']);
figure(fc2g); exportfig([exp_dir 'loc_av_con_gaussian_rf.pdf']);
figure(fo1g); exportfig([exp_dir 'loc_av_ori_gaussian_an.pdf']);
figure(fo2g); exportfig([exp_dir 'loc_av_ori_gaussian_rf.pdf']);

figure(fw1r); exportfig([exp_dir 'loc_av_width_rician_an.pdf']);
figure(fw2r); exportfig([exp_dir 'loc_av_width_rician_rf.pdf']);
figure(fc1r); exportfig([exp_dir 'loc_av_con_rician_an.pdf']);
figure(fc2r); exportfig([exp_dir 'loc_av_con_rician_rf.pdf']);
figure(fo1r); exportfig([exp_dir 'loc_av_ori_rician_an.pdf']);
figure(fo2r); exportfig([exp_dir 'loc_av_ori_rician_rf.pdf']);
%%
kernel_pts = 100;

[widths_x, widths_y] = kernel_smoother(true_widths, abs(ori_errors_signal)*180/pi, kernel_pts);
[cons_x, cons_y] = kernel_smoother(true_cons, abs(ori_errors_signal)*180/pi, kernel_pts);
[oris_x, oris_y] = kernel_smoother(true_oris, abs(ori_errors_signal)*180/pi, kernel_pts);
[widths_xrf, widths_yrf] = kernel_smoother(true_widths, abs(ori_errors_signal_g2)*180/pi, kernel_pts);
[cons_xrf, cons_yrf] = kernel_smoother(true_cons, abs(ori_errors_signal_g2)*180/pi, kernel_pts);
[oris_xrf, oris_yrf] = kernel_smoother(true_oris, abs(ori_errors_signal_g2)*180/pi, kernel_pts);

fw1g = figure; hold all; axis([1 8 0 30]); title('Error as a fucntion of line width, analytic prediction, Gaussian noise');
plot(widths_x, widths_y);
fw2g = figure; hold all; axis([1 8 0 30]); title('Error as a fucntion of line width, RF prediction, Gaussian noise');
plot(widths_xrf, widths_yrf);

fc1g = figure; hold all; axis([1 8 0 30]); title('Error as a fucntion of line contrast, analytic prediction, Gaussian noise');
plot(cons_x, cons_y);
fc2g = figure; hold all; axis([1 8 0 30]);title('Error as a fucntion of line contrast, RF prediction, Gaussian noise');
plot(cons_xrf, cons_yrf);

fo1g = figure; hold all; axis([0 pi 0 30]); title('Error as a fucntion of line orientation, analytic prediction, Gaussian noise');
plot(oris_x, oris_y);
fo2g = figure; hold all; axis([0 pi 0 30]); title('Error as a fucntion of line orientation, RF prediction, Gaussian noise');
plot(oris_xrf, oris_yrf);

fw1r = figure; hold all; axis([1 8 0 30]); title('Error as a fucntion of line width, analytic prediction, Rician noise');
plot(widths_x, widths_y);
fw2r = figure; hold all; axis([1 8 0 30]); title('Error as a fucntion of line width, RF prediction, Rician noise');
plot(widths_xrf, widths_yrf);

fc1r = figure; hold all; axis([1 8 0 30]); title('Error as a fucntion of line contrast, analytic prediction, Rician noise');
plot(cons_x, cons_y);
fc2r = figure; hold all; axis([1 8 0 30]);title('Error as a fucntion of line contrast, RF prediction, Rician noise');
plot(cons_xrf, cons_yrf);

fo1r = figure; hold all; axis([0 pi 0 30]); title('Error as a fucntion of line orientation, analytic prediction, Rician noise');
plot(oris_x, oris_y);
fo2r = figure; hold all; axis([0 pi 0 30]); title('Error as a fucntion of line orientation, RF prediction, Rician noise');
plot(oris_xrf, oris_yrf);
%
for ss = 1:4

    [widths_x, widths_y] = kernel_smoother(true_widths, abs(ori_errors_gaussian(:,ss))*180/pi, kernel_pts);
    [cons_x, cons_y] = kernel_smoother(true_cons, abs(ori_errors_gaussian(:,ss))*180/pi, kernel_pts);
    [oris_x, oris_y] = kernel_smoother(true_oris, abs(ori_errors_gaussian(:,ss))*180/pi, kernel_pts);
    
    [widths_xrf, widths_yrf] = kernel_smoother(true_widths, abs(ori_errors_gaussian_g2(:,ss))*180/pi, kernel_pts);
    [cons_xrf, cons_yrf] = kernel_smoother(true_cons, abs(ori_errors_gaussian_g2(:,ss))*180/pi, kernel_pts);
    [oris_xrf, oris_yrf] = kernel_smoother(true_oris, abs(ori_errors_gaussian_g2(:,ss))*180/pi, kernel_pts);
    
    figure(fw1g); plot(widths_x, widths_y);
    figure(fw2g); plot(widths_xrf, widths_yrf);
    
    figure(fc1g); plot(cons_x, cons_y);
    figure(fc2g); plot(cons_xrf, cons_yrf);
    
    figure(fo1g); plot(oris_x, oris_y);
    figure(fo2g); plot(oris_xrf, oris_yrf);
    
    [widths_x, widths_y] = kernel_smoother(true_widths, abs(ori_errors_rician(:,ss))*180/pi, kernel_pts);
    [cons_x, cons_y] = kernel_smoother(true_cons, abs(ori_errors_rician(:,ss))*180/pi, kernel_pts);
    [oris_x, oris_y] = kernel_smoother(true_oris, abs(ori_errors_rician(:,ss))*180/pi, kernel_pts);
    
    [widths_xrf, widths_yrf] = kernel_smoother(true_widths, abs(ori_errors_rician_g2(:,ss))*180/pi, kernel_pts);
    [cons_xrf, cons_yrf] = kernel_smoother(true_cons, abs(ori_errors_rician_g2(:,ss))*180/pi, kernel_pts);
    [oris_xrf, oris_yrf] = kernel_smoother(true_oris, abs(ori_errors_rician_g2(:,ss))*180/pi, kernel_pts);
    
    figure(fw1r); plot(widths_x, widths_y);
    figure(fw2r); plot(widths_xrf, widths_yrf);
    
    figure(fc1r); plot(cons_x, cons_y);
    figure(fc2r); plot(cons_xrf, cons_yrf);
    
    figure(fo1r); plot(oris_x, oris_y);
    figure(fo2r); plot(oris_xrf, oris_yrf);
end

figure(fw1g); exportfig([exp_dir 'kernel_width_gaussian_an.pdf']);
figure(fw2g); exportfig([exp_dir 'kernel_width_gaussian_rf.pdf']);
figure(fc1g); exportfig([exp_dir 'kernel_con_gaussian_an.pdf']);
figure(fc2g); exportfig([exp_dir 'kernel_con_gaussian_rf.pdf']);
figure(fo1g); exportfig([exp_dir 'kernel_ori_gaussian_an.pdf']);
figure(fo2g); exportfig([exp_dir 'kernel_ori_gaussian_rf.pdf']);

figure(fw1r); exportfig([exp_dir 'kernel_width_rician_an.pdf']);
figure(fw2r); exportfig([exp_dir 'kernel_width_rician_rf.pdf']);
figure(fc1r); exportfig([exp_dir 'kernel_con_rician_an.pdf']);
figure(fc2r); exportfig([exp_dir 'kernel_con_rician_rf.pdf']);
figure(fo1r); exportfig([exp_dir 'kernel_ori_rician_an.pdf']);
figure(fo2r); exportfig([exp_dir 'kernel_ori_rician_rf.pdf']);

%%
%Finally load in the naive predictor of orientation
naive_predictor = u_load([data_dir '\rfs\6098\predictor.mat']);
naive_predictor.tree_root = 'C:\isbe\asymmetry_project\data\synthetic_lines\increasing_noise_exp\rfs\';
save([data_dir '\rfs\6098\predictor.mat'], 'naive_predictor');
load([exp_dir 'workspace.mat'], 'rician_data_g2', 'true_*');
load([exp_dir 'ori_errors.mat']);
%
%Compute errors for rician data with s = 1.0
ori_rf = random_forest_reg_predict(naive_predictor, rician_data_g2(:,:,3), 1);
ori_errors_naive_rf = ori_error(ori_rf, true_oris);
save([exp_dir 'ori_errors.mat'], 'ori_errors_*');
%
%Plot error as a function of orientation, showing how the naive RF fails at
%predicting angle at the 0/180 cut point
[oris_xn, oris_yn] = kernel_smoother(true_oris*180/pi, abs(ori_errors_naive_rf)*180/pi, 100);
[oris_xr, oris_yr] = kernel_smoother(true_oris*180/pi, abs(ori_errors_rician_g2(:,3))*180/pi, 100);

figure; hold all;
plot(oris_xn, oris_yn);
plot(oris_xr, oris_yr);
legend({'Naive regression', 'Circular regression'});
title({'Orientation error as a function of predicted orientation - naive vs circular RF regression';...
    'Rician noise, SNR = 1.00'});
exportfig([exp_dir 'naive_vs_circular_rf.pdf']);
%%
%Do a comparison to the dual-tree, we'll probably move this later
dt_predictor = cell(2,1);
dt_predictor{1} = [data_dir 'rfs\10210\predictor.mat'];
dt_predictor{2} = [data_dir 'rfs\10211\predictor.mat'];
load([exp_dir 'workspace.mat'], 'rician_data_dt', 'true_*');
load([exp_dir 'ori_errors.mat']);

%Compute errors for rician data with s = 1.0
ori_errors_rician_dt = nan(size(ori_errors_rician_g2));
for ii = 3:4
    predictor = u_load(dt_predictor{ii-2});
    predictor.tree_root = [data_dir 'rfs\'];
    save(dt_predictor{ii-2}, 'predictor');
    
    ori_rf = random_forest_reg_predict(predictor, rician_data_dt(:,:,ii), 1);
    ori_errors_rician_dt(:,ii) = ori_error(ori_rf, true_oris);
end
save([exp_dir 'ori_errors.mat'], 'ori_errors_*');
%%
load([exp_dir 'workspace.mat'], 'true_*');
load([exp_dir 'ori_errors.mat']);
sorted_errors_rician_g2 = sort(abs(ori_errors_rician_g2(:,3:4))) * 180 / pi;
sorted_errors_rician_dt = sort(abs(ori_errors_rician_dt(:,3:4))) * 180 / pi;

xi = [1 100:100:20000];
xp = xi / 20000;
%
figure; hold all; title('DT vs G", CDF');
plot([sorted_errors_rician_g2(xi,:) sorted_errors_rician_dt(xi,:)], xp);
axis([0 45 0 1]);
legend({'G", Rician noise (1.0)', 'G", Rician noise (2.0)', 'DT, Rician noise (1.0)', 'DT, Rician noise (2.0)'},...
    'location', 'southeast');
%%
kernel_pts = 100;

for ii = 1:2
    [widths_g2_x(ii,:), widths_g2_y(ii,:)] = kernel_smoother(true_widths, abs(ori_errors_rician_g2(:,ii+2))*180/pi, kernel_pts);
    [cons_g2_x(ii,:), cons_g2_y(ii,:)] = kernel_smoother(true_cons, abs(ori_errors_rician_g2(:,ii+2))*180/pi, kernel_pts);
    [oris_g2_x(ii,:), oris_g2_y(ii,:)] = kernel_smoother(true_oris, abs(ori_errors_rician_g2(:,ii+2))*180/pi, kernel_pts);

    [widths_dt_x(ii,:), widths_dt_y(ii,:)] = kernel_smoother(true_widths, abs(ori_errors_rician_dt(:,ii+2))*180/pi, kernel_pts);
    [cons_dt_x(ii,:), cons_dt_y(ii,:)] = kernel_smoother(true_cons, abs(ori_errors_rician_dt(:,ii+2))*180/pi, kernel_pts);
    [oris_dt_x(ii,:), oris_dt_y(ii,:)] = kernel_smoother(true_oris, abs(ori_errors_rician_dt(:,ii+2))*180/pi, kernel_pts);
end
%%
figure; hold all; axis([1 8 0 30]); title('Error as a function of line width');
plot(widths_g2_x', widths_g2_y');
plot(widths_dt_x', widths_dt_y');
legend({'G", Rician noise (1.0)', 'G", Rician noise (2.0)', 'DT, Rician noise (1.0)', 'DT, Rician noise (2.0)'},...
    'location', 'northeast');
figure; hold all; axis([1 8 0 30]); title('Error as a function of line contrast');
plot(cons_g2_x', cons_g2_y');
plot(cons_dt_x', cons_dt_y');
legend({'G", Rician noise (1.0)', 'G", Rician noise (2.0)', 'DT, Rician noise (1.0)', 'DT, Rician noise (2.0)'},...
    'location', 'northeast');
figure; hold all; axis([0 pi 0 30]); title('Error as a function of line orientation');
plot(oris_g2_x', oris_g2_y');
plot(oris_dt_x', oris_dt_y');
legend({'G", Rician noise (1.0)', 'G", Rician noise (2.0)', 'DT, Rician noise (1.0)', 'DT, Rician noise (2.0)'},...
    'location', 'northeast');


