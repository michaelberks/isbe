%--------------------------------------------------------------------------
% Script for producing results for journal submission
% Experiment predicting orientation of synthetic lines in the presence of increasing noise 

%--------------------------------------------------------------------------
%%
%0. Generic arguments to initialise
data_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\edge_vs_line\';
exp_dir = 'C:\isbe\asymmetry_project\experiments\synthetic_lines\edge_vs_line\';

dx = 64;
dy = 64;
cx = 32;
cy = 32;
snr = 1.00;
warning('off', 'ASYM:unexpectedArgument');
%%
%1. Train a set of random forests on line images
%%
% DECOMP_TYPE="g2d" WIN_SIZE=1 OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" BG_TYPE="edge" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
% DECOMP_TYPE="g2d" WIN_SIZE=3 OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" BG_TYPE="edge" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
% DECOMP_TYPE="dt" WIN_SIZE=1 OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" BG_TYPE="edge" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
% DECOMP_TYPE="dt" WIN_SIZE=3 OUTPUT_TYPE="orientation" PREDICTION_TYPE="rf_regression" BG_TYPE="edge" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh

DECOMP_TYPE="g2d" WIN_SIZE=1 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="edge" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="g2d" WIN_SIZE=3 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="edge" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="dt" WIN_SIZE=1 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="edge" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="dt" WIN_SIZE=3 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="edge" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh

DECOMP_TYPE="g2d" WIN_SIZE=1 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="g2d" WIN_SIZE=3 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="dt" WIN_SIZE=1 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="dt" WIN_SIZE=3 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="flat" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh

DECOMP_TYPE="g12d" WIN_SIZE=3 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="edge" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh
DECOMP_TYPE="g2dh" WIN_SIZE=3 OUTPUT_TYPE="centre_detection" PREDICTION_TYPE="rf_classification" BG_TYPE="edge" NUM_TREES=10 IMAGE_TYPE="line" NUM_SAMPLES=100000 PTS_PER_IMAGE=50 NOISE_TYPE="rician" NOISE_PARAMS="[1.00]" MODEL_ROOT="models/synthetic_lines" qsub -V -l twoday -t 1-10 matlab_code/trunk/hydra/cuc/build_predictor.sh

MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10273'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10273 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10274'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10274 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10275'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10275 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10276'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10276 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10277'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10277 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10278'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10278 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10279'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10279 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10280'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10280 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10449'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10449 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10450'" MAKE_SAMPLED_MAPS=0 qsub -V -hold_jid 10450 -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh

%%
%2. Make a set of test images
mkdir([data_dir 'images']);
mkdir([data_dir 'orientations']);
mkdir([data_dir 'widths']);
mkdir([data_dir 'line_labels']);
mkdir([data_dir 'edge_labels']);
mkdir([data_dir 'fov_masks']);

for ii = 1:100
    edge_width = sample_uniform([1 8]);
    line_width = sample_uniform([1 8]);
    edge_contrast = sample_uniform([4 8]);
    line_contrast = sample_uniform([4 8]);
    line_ori = sample_uniform([0 180]);
    edge_ori = sample_from_normal(line_ori + 90, 30^2, 1);

    [line, line_label, ~, orientation_map] =...
        create_ellipse_bar(line_width/2, line_contrast, line_ori, dy, dx, cx, cy);
    width_map = line_label * line_width;
    
    [edge, edge_label] = create_sin_step(edge_width/2, edge_contrast, edge_ori, dy, dx, cx, cy);
    
    signal = 1 + line + edge;
    test_image = ricernd(signal, snr);
    fov_mask = edge_label | line_label;
    
    save([data_dir 'images\image' zerostr(ii, 4) '.mat'],...
        'test_image');
    save([data_dir 'orientations\image' zerostr(ii, 4) '_ori.mat'],...
        'orientation_map');
    save([data_dir 'widths\image' zerostr(ii, 4) '_width.mat'],...
        'width_map');
    save([data_dir 'line_labels\image' zerostr(ii, 4) '_label.mat'],...
        'line_label');
    save([data_dir 'edge_labels\image' zerostr(ii, 4) '_label.mat'],...
        'edge_label');
    save([data_dir 'fov_masks\image' zerostr(ii, 4) '_mask.mat'],...
        'fov_mask');
end
%%
%3. Predict line presence on the images using the CSF
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10273'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10290 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10274'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10291 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10275'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10292 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10276'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10293 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10277'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10294 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10278'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10295 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10279'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10296 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10280'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10297 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10449'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10451 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh
MODEL_ROOT="models/synthetic_lines/centre_detection/rf_classification" MODEL_PATH="'10450'" NUM_JOBS=10 IMAGE_ROOT="synthetic_lines/edge_vs_line" USE_SAMPLED_MAPS=0 qsub -V -hold_jid 10452 -t 1-10 -l twoday matlab_code/trunk/hydra/cuc/predict_image_set.sh

%%
%4. Analyse results
pred_dir = [data_dir 'predictions\centre_detection\rf_classification\'];
label_dir = [data_dir 'line_labels\'];
fov_mask_dir = [data_dir 'fov_masks\'];

rf_codes = cell(8, 4);
rf_codes( 1,:) = {'10273', 'g2', '1', 'with edges'};
rf_codes( 2,:) = {'10274', 'g2', '3', 'with edges'};
rf_codes( 3,:) = {'10275', 'dt', '1', 'with edges'};
rf_codes( 4,:) = {'10276', 'dt', '3', 'with edges'};
rf_codes( 5,:) = {'10277', 'g2', '1', 'lines only'};
rf_codes( 6,:) = {'10278', 'g2', '3', 'lines only'};
rf_codes( 7,:) = {'10279', 'dt', '1', 'lines only'};
rf_codes( 8,:) = {'10280', 'dt', '3', 'lines only'};
rf_codes( 9,:) = {'10449', 'gh2', '3', 'lines only'};
rf_codes(10,:) = {'10450', 'gh2', '3', 'lines only'};
%%
%
%4.a compute ROC curves and Az values for the whole images
roc_pts = [];
auc = [];
for ii = 1:10
    [roc_pts(:,:,ii), auc(ii,1)] = compute_roc_image_set([pred_dir rf_codes{ii,1} '\'], label_dir, [], 1); %#ok
end

%4.b compute ROC curves and Az values for just the areas around the line and edge
roc_pts_roi = [];
auc_roi = [];
for ii = 1:10
    [roc_pts_roi(:,:,ii), auc_roi(ii,1)] = compute_roc_image_set([pred_dir rf_codes{ii,1} '\'], label_dir, fov_mask_dir, 1); %#ok
end
%%
%4.c Visually inspect the prediction images
for ii = 1:20
    %test_im = u_load([data_dir 'images\image' zerostr(ii,4) '.mat']);
    figure;
    for jj = 1:8
        pred_im = load_uint8([data_dir 'predictions\centre_detection\rf_classification\' rf_codes{jj,1} '\image' zerostr(ii,4) '_pred.mat']);
        subplot(2,4,jj); imgray(pred_im); caxis([0 1]);
        title(['Decomp: ' rf_codes{jj,2} ', w = ' rf_codes{jj,3} ', ' rf_codes{jj,4}]); 
    end
    %colormap(jet(256));
    %exportfig([exp_dir 'predictions ' zerostr(ii,3) '_mp.png']);
end
%%
ii = 7;
test_im = u_load([data_dir 'images\image' zerostr(ii,4) '.mat']);
write_im_from_colormap(test_im, 'K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\line_edge.png');
for jj = [2 4]
    pred_im = load_uint8([data_dir 'predictions\centre_detection\rf_classification\' rf_codes{jj,1} '\image' zerostr(ii,4) '_pred.mat']);
    write_im_from_colormap(pred_im, ['K:\isbe\conferences_and_symposia\Monday Meeting 12-10-29\figures\line_edge' rf_codes{jj,2} '.png']);
end
%%
%5. Inspect tree outputs
%
counts = zeros(1,48);
for ii = 1:100
    load(['C:\isbe\asymmetry_project\data\models\synthetic_lines\detection\rf_classification\dt\01_trees\rf_tree' zerostr(ii,4) '.mat']);
    c = hist(tree.var(tree.var > 0), 1:48);
    counts = counts + c;
end
figure; bar(1:48, counts);
%%
