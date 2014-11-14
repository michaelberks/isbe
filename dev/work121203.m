%%
%It's something to do with mage size I think...
%%
retroot = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\';

ret1 = u_load([retroot 'images\01_test.mat']);
f_mask1 = u_load([retroot 'fov_masks\01_test_f_mask.mat']);
v_mask1 = u_load([retroot 'vessel_masks\01_test_v_mask.mat']);
pred1 = u_load([retroot 'predictions\detection\rf_classification\9162\01_test_pred.mat']);

rf = u_load('C:\isbe\asymmetry_project\data\models\vessel\detection\rf_classification\9162\predictor.mat');
rf.tree_root = 'C:\isbe\asymmetry_project\data\models\vessel\detection\rf_classification\';
job_args = u_load('C:\isbe\asymmetry_project\data\models\vessel\detection\rf_classification\9162\job_args.mat');
job_args.decomposition_args.decomp_type = {job_args.decomposition_args.decomp_type};


%%
f_mask2 = f_mask1;
f_mask2(end+1,end+1) = 0;
f_mask2(1,:) = []; f_mask2(:,1) = [];

v_mask2 = v_mask1;
v_mask2(end+1,end+1) = 0;
v_mask2(1,:) = []; v_mask2(:,1) = [];

ret2 = ret1;
ret2(end+1,end+1,:) = 0;
ret2(1,:,:) = []; ret2(:,1,:) = [];
%%
f_mask3 = padarray(f_mask1, [1 1], 'pre');
f_mask3(end,:) = []; f_mask3(:,end) = [];

v_mask3 = padarray(v_mask1, [1 1], 'pre');
v_mask3(end,:) = []; v_mask3(:,end) = [];

ret3 = padarray(ret1, [1 1], 'pre');
ret3(end,:,:) = []; ret3(:,end,:) = [];

figure;
subplot(1,2,1); imgray(ret1);
subplot(1,2,2); imgray(ret3);
%%
% f_mask4 = padarray(f_mask1, [1 1], 'pre');
% v_mask4 = padarray(v_mask1, [1 1], 'pre');
%ret4 = padarray(ret1, [1 1], 'pre');
clear ret4 f_mask4 v_mask4
f_mask4 = [zeros(586,1) [zeros(1,565); f_mask1; zeros(1,565)] zeros(586,1)];
v_mask4 = [zeros(586,1) [zeros(1,565); v_mask1; zeros(1,565)] zeros(586,1)];
ret4(:,:,1) = [zeros(586,1) [zeros(1,565); ret1(:,:,1); zeros(1,565)] zeros(586,1)];
ret4(:,:,2) = [zeros(586,1) [zeros(1,565); ret1(:,:,2); zeros(1,565)] zeros(586,1)];
ret4(:,:,3) = [zeros(586,1) [zeros(1,565); ret1(:,:,3); zeros(1,565)] zeros(586,1)];
figure;
subplot(1,2,1); imgray(ret1);
subplot(1,2,2); imgray(ret4);
%%
f_mask5 = f_mask1;
f_mask5(end,:) = []; f_mask5(:,end) = [];

v_mask5 = v_mask1;
v_mask5(end,:) = []; v_mask5(:,end) = [];

ret5 = ret1;
ret5(end,:,:) = []; ret5(:,end,:) = [];

figure;
subplot(1,2,1); imgray(ret1);
subplot(1,2,2); imgray(ret5);
%%


%%
[pred4] = predict_image(...
    'image_in', ret4,...
    'decomposition_args', job_args.decomposition_args,...
    'predictor', rf, ...
    'prediction_type', 'rf_classification',...
    'output_type', 'detection',...
    'use_probs', 0,...
    'mask', f_mask4,...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 128,...
    'incremental_results', 0);


%
[roc_pts1 auc1 tp_count1 fp_count1 auc_se1 n_pos1 n_neg1] = calculate_roc_image(...
    pred1, v_mask1, (-1:100)/100, f_mask1);
[roc_pts2 auc2 tp_count2 fp_count2 auc_se2 n_pos2 n_neg2] = calculate_roc_image(...
    pred2, v_mask2, (-1:100)/100, f_mask2);
[roc_pts3 auc3 tp_count3 fp_count3 auc_se3 n_pos3 n_neg3] = calculate_roc_image(...
    pred3, v_mask3, (-1:100)/100, f_mask3);
[roc_pts4 auc4 tp_count4 fp_count4 auc_se4 n_pos4 n_neg4] = calculate_roc_image(...
    pred4, v_mask4, (-1:100)/100, f_mask4);
[roc_pts5 auc5 tp_count5 fp_count5 auc_se5 n_pos5 n_neg5] = calculate_roc_image(...
    pred5, v_mask5, (-1:100)/100, f_mask5);
%%
warning('off', 'ASYM:unexpectedArgument');
f_mask_small = false(size(f_mask1));
f_mask_small(129:256, 129:256) = 1;

for i_r = 16%:4
    for i_c = 0
        f_maskp = f_mask_small;
        f_maskp(end+i_r,end+i_c) = 0;
        
        v_maskp = v_mask1;
        v_maskp(end+i_r,end+i_c) = 0;
        
        retp = ret1;
        retp(end+i_r,end+i_c,:) = 0;
        
        [predp] = predict_image(...
            'image_in', retp,...
            'decomposition_args', job_args.decomposition_args,...
            'predictor', rf, ...
            'prediction_type', 'rf_classification',...
            'output_type', 'detection',...
            'use_probs', 0,...
            'mask', f_maskp,...
            'tree_mask', [], ...
            'num_trees', [], ...
            'max_size', 128,...
            'incremental_results', 0);

        %
        [roc_ptsp aucp] = calculate_roc_image(...
            predp, v_maskp, (-1:100)/100, f_maskp);
        
        display(['R = ' num2str(i_r) ', C = ' num2str(i_c) ', auc = ' num2str(aucp)]);
    end
end
%%
d_args{1}.decomp_type = {'dt'};
d_args{1}.feature_shape = 'rect';
d_args{1}.feature_type = 'complex';
d_args{1}.do_max = 0;
d_args{1}.rotate = 0;
d_args{1}.use_nag = 0;
d_args{1}.levels = 1:3;
d_args{1}.win_size = 1;
d_args{1}.normalise = 0;
d_args{1}.pca = [];
for off = 1:8
    a1 = zeros(64+off,64);
    r1 = compute_filter_responses(a1, d_args{1});
    sample_image_features(r1, 32, 32, d_args{1});
end
%%
warning('off', 'ASYM:unexpectedArgument');
d_args{1}.decomp_type = {'dt'};
d_args{1}.feature_shape = 'rect';
d_args{1}.feature_type = 'conj';
d_args{1}.do_max = 0;
d_args{1}.rotate = 0;
d_args{1}.use_nag = 0;
d_args{1}.levels = 1:3;
d_args{1}.unwrap_phase = 1;
d_args{1}.interp_mag_phase = 0;

d_args{2}.decomp_type = {'gabor'};
d_args{2}.num_angles = 6;
d_args{2}.sigma_range = [1 2 4];	
d_args{2}.do_max = 0;
d_args{2}.rotate = 0;
d_args{2}.feature_type = 'conj';

for i_d = 1:2
    d_args{i_d}.win_size = 3;
    d_args{i_d}.normalise = 0;
    d_args{i_d}.pca = [];
end

a1 = create_gauss_bar(4, 1, 22.5, 128, 128, 64, 64);
%%
times = zeros(8,2);
mean_diffs = zeros(8,2);
max_diffs = zeros(8,2);
median_diffs = zeros(8,2);
do_plot = 0; 
for off = 1:8
    %a2 = padarray(a1, [off off]);
    a2 = [[a1(off+1:end,off+1:end); zeros(off,128-off)] zeros(128,off)];

    if do_plot
        figure;
    end
    for i_d = 1:2
        r1 = compute_filter_responses(a1, d_args{i_d});
        r2 = compute_filter_responses(a2, d_args{i_d});

        tic;
        f1 = sample_image_features(r1, 64, 64, d_args{i_d});        
        f2 = sample_image_features(r2, 64-off, 64-off, d_args{i_d});
        times(off,i_d) = toc;
        diffs = abs(f1 - f2);
        mean_diffs(off,i_d) = mean(diffs);
        max_diffs(off,i_d) = max(diffs);
        median_diffs(off,i_d) = median(diffs);
        
        if do_plot
            subplot(1,2,1); plot(f1(1:162), f2(1:162), 'x'); axis equal; hold all;
            subplot(1,2,2); plot(f1(163:324), f2(163:324), 'x'); axis equal; hold all;
        end
                
%         i_p = 1;
%         for i_x = -1:1
%             for i_y = -1:1
%                 f2 = sample_image_features(r2, 64-off+i_y, 64-off+i_x, d_args{i_d});
%                 subplot(3,3,i_p); plot(angle(f1), angle(f2), 'x'); hold all;
%                 i_p = i_p + 1;
%             end
%         end
    end
    if do_plot
        legend({'DT samples', 'Gabor samples'});
    end
end
%%
dt_args{1}.decomp_type = {'dt'};
dt_args{1}.feature_shape = 'rect';
dt_args{1}.feature_type = 'complex';
dt_args{1}.do_max = 0;
dt_args{1}.rotate = 0;
dt_args{1}.use_nag = 0;
dt_args{1}.levels = 1:2;
dt_args{1}.win_size = 1;
dt_args{1}.normalise = 0;
dt_args{1}.pca = [];

a1 = create_gauss_bar(4, 1, 22.5, 64, 64, 32, 32);
[cc rr] = meshgrid(1:64, 1:64);
r1 = compute_filter_responses(a1, dt_args{1});
f1 = sample_image_features(r1, rr(:), cc(:), dt_args{1});
dt = reshape(f1, 64, 64, 6, 2);
%%
for i_l = 1:2
    figure;
    for i_b = 1:6
        subplot(2,3,i_b); image(complex2rgb(dt(:,:,i_b,i_l),[-pi pi], max(abs(f1(:)))));
    end
end
%%
differences = zeros(20);
for ii = 1:20
    ret1 = u_load(['C:\isbe\asymmetry_project\data\sanity_check\images\' zerostr(ii,2) '_test.mat']);
    for jj = 21:40
        ret2 = u_load(['C:\isbe\asymmetry_project\data\sanity_check\images\' zerostr(jj,2) '_training.mat']);
        
        differences(ii,jj-20) = mean(abs(ret1(:) - ret2(:)));
    end
end
%%
rf_codes = cell(4,2);
rf_codes( 1,:) = {'rf_regression\88171', 'Gaussian (G,H) + RF (mae='};
rf_codes( 2,:) = {'rf_regression\24308', 'Gabor (Re,Im) + RF (mae='};
rf_codes( 3,:) = {'rf_regression\42813',  'DT-CWT (Re,Im) + RF (mae='};
rf_codes( 4,:) = {'rf_regression\13295', 'Monogenic + RF (mae='};

retroot = 'C:\isbe\asymmetry_project\data\retinograms\DRIVE_clean\test\';
pred_dir = [retroot 'predictions\orientation\'];

for ii = 3
    figure;
    for jj = 1:4
        ori_map = u_load([pred_dir rf_codes{jj,1} '\' zerostr(ii,2) '_test_pred.mat']);
        subplot(2,2,jj); imgray(complex2rgb(ori_map));
    end
end
%%
R = 4, C = 4, auc = 0.91634
R = 0, C = 0, auc = 0.97256
R = 0, C = 1, auc = 0.97256
R = 1, C = 0, auc = 0.88255
R = 1, C = 1, auc = 0.88255
R = 2, C = 0, auc = 0.88255
R = 2, C = 1, auc = 0.88255
R = 3, C = 0, auc = 0.94543
R = 3, C = 1, auc = 0.94543
R = 4, C = 0, auc = 0.94543
R = 4, C = 1, auc = 0.94543
R = 0, C = 2, auc = 0.96778
R = 0, C = 3, auc = 0.96778
R = 0, C = 4, auc = 0.95225
R = 0, C = 5, auc = 0.95225
R = 0, C = 6, auc = 0.96653
R = 0, C = 7, auc = 0.96653
R = 0, C = 8, auc = 0.9713
%%
tt = zeros(100,4);
cc = 1;
for aa = 0:1
    d_args{1}.unwrap_phase = aa;
	
    for bb = 0:1
        d_args{1}.interp_mag_phase = bb;
        for ii = 1:100
            tic;
            f1 = sample_image_features(r1, 64, 64, d_args{1});
            tt(ii,cc) = toc;
        end
        cc = cc + 1;
        
    end
end
        
