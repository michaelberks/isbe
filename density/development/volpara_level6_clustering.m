base_dir = 'A:\D_2.2_Datasets\125_375_Single_MLO_Dataset_VOLPARA\';
control_list = dir([base_dir 'CONTROLS\ASSURE_CONTROLS*']);
cancer_list = dir([base_dir 'CANCERS\ASSURE_CANCERS*']);

num_controls = length(control_list);
num_cancers = length(cancer_list);
%
control_im_names = cell(num_controls,1);
control_mask_names = cell(num_controls,1);
for i_case = 1:num_controls
    im_name = dir([base_dir 'CONTROLS\' control_list(i_case).name '\*ML*hint_densityMap*.pgm']);
    mask_name = dir([base_dir 'CONTROLS\' control_list(i_case).name '\*ML*hint_Segmentation.pgm']);
    
    control_im_names{i_case} = [base_dir 'CONTROLS\' control_list(i_case).name '\' im_name(end).name];
    control_mask_names{i_case} = [base_dir 'CONTROLS\' control_list(i_case).name '\' mask_name(end).name];
end

cancer_im_names = cell(num_cancers,1);
cancer_mask_names = cell(num_cancers,1);
for i_case = 1:num_cancers
    im_name = dir([base_dir 'CANCERS\' cancer_list(i_case).name '\*ML*hint_densityMap*.pgm']);
    mask_name = dir([base_dir 'CANCERS\' cancer_list(i_case).name '\*ML*hint_Segmentation.pgm']);
    
    cancer_im_names{i_case} = [base_dir 'CANCERS\' cancer_list(i_case).name '\' im_name(end).name];
    cancer_mask_names{i_case} = [base_dir 'CANCERS\' cancer_list(i_case).name '\' mask_name(end).name];
end
    
case_im_names = [control_im_names; cancer_im_names];
case_mask_names = [control_mask_names; cancer_mask_names];

num_cases = length(case_im_names);
rand_i = randperm(num_cases);

num_train = num_cases*0.8;
num_test = num_cases - num_train;
train_case_idx = rand_i(1:num_train);
test_case_idx = rand_i(num_train+1:end);


%%
pts_per_img = 100;
total_samples = num_train*pts_per_img;

train_responses = zeros(total_samples, 12, 2);
train_xy = zeros(total_samples, 2, 2);
train_thickness = zeros(total_samples, 1);
train_s = zeros(total_samples, 1, 2);
train_phi = zeros(total_samples, 1, 2);
curr_sample = 0;

decomposition_args.decomp_type = 'dt';      %Use the dual-tree
decomposition_args.win_size = 1;            %Window size about pixel to sample features
decomposition_args.levels = 6;            %Levels to include
decomposition_args.feature_shape = 'rect';  %Keep this as rect, although there are other things you can play with (see sample_dt_data)
decomposition_args.feature_type = 'conj';   %Again, 'conj' works best, but can use 'mag' for magnitude only, or 'phase' (see convert_complex_representation)
decomposition_args.do_max = 0;              %Keep all 6 sub-bands (0) or only keep the maximum (1)
decomposition_args.rotate = 0;              %Try and make rotation invariant (best left switched off)
decomposition_args.use_nag = 0;             %Use the NAG toolbox to interpolate (best left swicthed off)
decomposition_args.normalise = 0;
labelBackground = 0;
labelBreast = 2;
labelPectoral = 3;

for i_idx = 1:num_train
    display(['Sampling from image ' num2str(i_idx) ' of ' num2str(num_train)]);
    
    i_case = train_case_idx(i_idx);
    
    mam_vdm = load_volpara_vdm(case_im_names{i_case});
    mask = imread(case_mask_names{i_case});
    mask2 = mask;
    mask2(mask==4) = 2;
    
    dual_tree = compute_dual_tree(mam_vdm, 6);       

    objLM = landmarkDetection(mam_vdm, mask2, labelBackground, labelBreast, labelPectoral); 
    objLM = objLM.detectCSLandmarks(0); 
    objCS1 = BreastCoodSystem();  
    objCS1 = objCS1.setLabels(labelBackground, labelBreast, labelPectoral); 
    objCS1 = objCS1.mask2BreastCoodSystem(objLM.mask, objLM.nipple, objLM.pctLine);
        
    sample_idx = curr_sample+(1:pts_per_img);
    curr_sample = curr_sample + pts_per_img;
    
    for i_exp = 1:2
        
        if i_exp == 1
            [rows cols] = find(mask>-1);
        else
            [rows cols] = find(mask==4);
        end
        
        rand_i = randperm(length(rows),pts_per_img);
        rows = rows(rand_i);
        cols = cols(rand_i);       

        [responses] = compute_filter_responses(mam_vdm, decomposition_args);
        train_responses(sample_idx,:,i_exp) = sample_image_features(responses, rows, cols, decomposition_args); 
        train_xy(sample_idx,1,i_exp) = cols;
        train_xy(sample_idx,2,i_exp) = rows;
        
        [train_s(sample_idx,1,i_exp), train_phi(sample_idx,1,i_exp)] =...
            objCS1.cart2CoordSystem(rows, cols);
    end
    
    [~, fname, ~] = fileparts(case_im_names{i_case});
    s_idx = strfind(fname, 'Map_H_') + 6;
    e_idx = strfind(fname, 'mm') - 1;
    
    train_thickness(sample_idx) = str2double(fname(s_idx:e_idx));
end
%
opts = statset('MaxIter', 1000);
num_k = 8;

output_dir = 'C:\isbe\density\assure\level6_clustering\norm\';

[k_idx, k_centres, sumd] = kmeans(train_responses(:,:,1), num_k, 'EmptyAction', 'drop', 'Replicates', 10, 'Options', opts);
save([output_dir 'k_all.mat'], 'k_idx', 'k_centres', 'sumd');

[k_idx, k_centres, sumd] = kmeans(train_responses(:,:,2), num_k, 'EmptyAction', 'drop', 'Replicates', 10, 'Options', opts);
save([output_dir 'k_compressed.mat'], 'k_idx', 'k_centres', 'sumd'); 

save([output_dir 'k_samples.mat'], 'train_case_idx', 'test_case_idx',...
    'train_responses', 'train_xy', 'train_thickness', 'train_s', 'train_phi');
%%
k_data = load([output_dir 'k_all.mat']);
k_data(2) = load([output_dir 'k_compressed.mat']);
%%
band_angles = pi*(-75:30:75)/180;
band_x = cos(band_angles);
band_y = sin(band_angles);

colors = jet(8);
for i_exp = 1:2
    axis_lim = max(max(k_data(i_exp).k_centres(:,1:6)));
    figure;
    for i_k = 1:8
        subplot(2,4,i_k);
        axis equal ij; 
        axis([0 axis_lim -axis_lim axis_lim]); 
        axis off; hold all;
        plot([0 axis_lim], [0 0], 'k');
        plot([0 0], [-axis_lim axis_lim], 'k');

        for i_band = 1:6
            plot(...
                [0 band_x(i_band)*k_data(i_exp).k_centres(i_k,i_band)],...
                [0 band_y(i_band)*k_data(i_exp).k_centres(i_k,i_band)],...
                'color', colors(i_k,:), 'linewidth', 3);
        end
    end
end
%%
labelBackground = 0;
labelBreast = 2;
labelPectoral = 3;
        
colors = jet(num_k);
make_pictures = 1;
save_pictures = 1;
num_pictures = 10;

pts_per_img = 100;
test_samples = num_test*pts_per_img;
curr_sample = 0;

test_responses = zeros(test_samples, 12, 2);
test_xy = zeros(test_samples, 2, 2);
test_thickness = zeros(test_samples, 1);
test_s = zeros(test_samples, 1, 2);
test_phi = zeros(test_samples, 1, 2);
%
for i_idx = 1:num_test
    
    display(['Processing image ' num2str(i_idx) ' of ' num2str(num_test)]);
    
    i_case = test_case_idx(i_idx);
    
    mam_vdm = load_volpara_vdm(case_im_names{i_case});
    mask2 = mask;
    mask2(mask==4) = 2;
    
    dual_tree = compute_dual_tree(mam_vdm, 6);       

    objLM = landmarkDetection(mam_vdm, mask2, labelBackground, labelBreast, labelPectoral);
    objLM = objLM.detectCSLandmarks(0); 
    objCS1 = BreastCoodSystem();  
    objCS1 = objCS1.setLabels(labelBackground, labelBreast, labelPectoral); 
    objCS1 = objCS1.mask2BreastCoodSystem(objLM.mask, objLM.nipple, objLM.pctLine);

    sample_idx = curr_sample+(1:pts_per_img);
    curr_sample = curr_sample + pts_per_img;
        
    for i_exp = 1:2
        
        if i_exp == 1
            [rows cols] = find(mask>-1);
        else
            [rows cols] = find(mask==4);
        end
        
        [responses] = compute_filter_responses(mam_vdm, decomposition_args);
        sampled_features = sample_image_features(responses, rows, cols, decomposition_args);
        %
        if make_pictures && i_idx <= num_pictures
            k_dists = zeros(size(sampled_features,1), num_k);
            for i_k = 1:num_k
                k_dists(:,i_k) = sum(bsxfun(@minus, sampled_features, k_data(i_exp).k_centres(i_k,:)).^2,2);
            end
            [~,assigned_k] = min(k_dists,[],2);

            b_inner = bwboundaries(mask == 4);
            b_outer = bwboundaries(mask > 0);

            %
            figure; 
            subplot(1,2,1); imgray(mam_vdm);
            subplot(1,2,2); axis equal ij; hold all;
            for i_k = 1:num_k
                idx = assigned_k == i_k;
                plot(cols(idx), rows(idx), '.', 'markeredgecolor', colors(i_k,:));
            end
            plot(b_inner{1}(:,2), b_inner{1}(:,1), 'k', 'linewidth', 2);
            plot(b_outer{1}(:,2), b_outer{1}(:,1), 'k', 'linewidth', 2);

            set(gca, 'xlim', [0 size(mam_vdm,2)], 'ylim', [0 size(mam_vdm,1)]);
            
            if save_pictures
                if i_exp == 1
                    exportfig([output_dir 'k_all\test' zerostr(i_idx,2) '.png']);
                else
                    exportfig([output_dir 'k_compressed\test' zerostr(i_idx,2) '.png']);
                end
            end
        end

        rand_i = randperm(length(rows),pts_per_img);      
        test_responses(sample_idx,:,i_exp) = sampled_features(rand_i,:); 
        test_xy(sample_idx,1,i_exp) = cols(rand_i);
        test_xy(sample_idx,2,i_exp) = rows(rand_i);

        [test_s(sample_idx,1,i_exp), test_phi(sample_idx,1,i_exp)] =...
            objCS1.cart2CoordSystem(rows(rand_i), cols(rand_i));
          
    end
    
    [~, fname, ~] = fileparts(case_im_names{i_case});
    s_idx = strfind(fname, 'Map_H_') + 6;
    e_idx = strfind(fname, 'mm') - 1;

    test_thickness(sample_idx) = str2double(fname(s_idx:e_idx));  
end
save([output_dir 'k_samples.mat'], 'test_responses', 'test_xy', 'test_thickness', 'test_s', 'test_phi', '-append');
%%

rf_args.prediction_type = 'rf_regression';
rf_args.n_trees = 100;
rf_args.d = [];
rf_args.w_prior = 0;
rf_args.impure_thresh = 1.0000e-004;
rf_args.split_min = 100;
rf_args.end_cut_min = 25;
rf_args.do_ubound = 0;
rf_args.quiet = 1;
rf_args.do_circular = [];
rf_args.overwrite = 0;
rf_args.minimise_size = 0;
rf_args.split_criterion = 'ssq';
rf_args.var_criterion = 'ssq';

rf_args.sampling_args.sampling_method = 'sample_saved_training_data';
rf_args.decomposition_args = [];
%
rf_dir = 'C:\isbe\density\assure\375_125\models\thickness_regression\norm_all\';
rf_args.predictor_path = [rf_dir 'rf.mat'];
rf_args.tree_dir = [rf_dir 'trees\'];
rf_args.sampling_args.y = train_thickness;
rf_args.sampling_args.X = train_responses(:,:,1);
rf_thickness_all = random_forest_reg_train(rf_args);
%
rf_dir = 'C:\isbe\density\assure\375_125\models\thickness_regression\norm_comp\';
rf_args.predictor_path = [rf_dir 'rf.mat'];
rf_args.tree_dir = [rf_dir 'trees\'];
rf_args.sampling_args.y = train_thickness;
rf_args.sampling_args.X = train_responses(:,:,2);
rf_thickness_comp = random_forest_reg_train(rf_args);
%
thickness_p_all = random_forest_reg_predict(rf_thickness_all, test_responses(:,:,1), 0);
thickness_p_comp = random_forest_reg_predict(rf_thickness_comp, test_responses(:,:,1), 0);

t_lims = [min(test_thickness) max(test_thickness)];
figure; 
subplot(1,2,1); axis equal; hold on;
plot(t_lims, t_lims, 'k');
plot(test_thickness, thickness_p_all, 'r.');
axis([t_lims t_lims]);
title('Predicting breast thickness - all pixels');
xlabel('Thickness (mm) - ground truth');
ylabel('Predicted thickness');

subplot(1,2,2); axis equal; hold on;
plot(t_lims, t_lims, 'k');
plot(test_thickness, thickness_p_comp, 'r.');
axis([t_lims t_lims]);
title('Predicting breast thickness - compressed region');
xlabel('Thickness (mm) - ground truth');
ylabel('Predicted thickness');
%%
rf_dir = 'C:\isbe\density\assure\375_125\models\position_regression\norm_y\all\';
rf_args.predictor_path = [rf_dir 'rf.mat'];
rf_args.tree_dir = [rf_dir 'trees\'];
rf_args.sampling_args.y = train_xy(:,2,1);
rf_args.sampling_args.X = train_responses(:,:,1);
rf_pos_y_all = random_forest_reg_train(rf_args);
%
rf_dir = 'C:\isbe\density\assure\375_125\models\position_regression\norm_y\comp\';
rf_args.predictor_path = [rf_dir 'rf.mat'];
rf_args.tree_dir = [rf_dir 'trees\'];
rf_args.sampling_args.y = train_xy(:,2,2);
rf_args.sampling_args.X = train_responses(:,:,2);
rf_pos_y_comp = random_forest_reg_train(rf_args);
%
rf_dir = 'C:\isbe\density\assure\375_125\models\position_regression\norm_x\all\';
rf_args.predictor_path = [rf_dir 'rf.mat'];
rf_args.tree_dir = [rf_dir 'trees\'];
rf_args.sampling_args.y = train_xy(:,1,1);
rf_args.sampling_args.X = train_responses(:,:,1);
rf_pos_x_all = random_forest_reg_train(rf_args);
%
rf_dir = 'C:\isbe\density\assure\375_125\models\position_regression\norm_x\comp\';
rf_args.predictor_path = [rf_dir 'rf.mat'];
rf_args.tree_dir = [rf_dir 'trees\'];
rf_args.sampling_args.y = train_xy(:,1,2);
rf_args.sampling_args.X = train_responses(:,:,2);
rf_pos_x_comp = random_forest_reg_train(rf_args);
%
y_p_all = random_forest_reg_predict(rf_pos_y_all, test_responses(:,:,1), 0);
y_p_comp = random_forest_reg_predict(rf_pos_y_comp, test_responses(:,:,2), 0);
%
y_lims = [min(test_xy(:,2,1)) max(test_xy(:,2,1))];
figure; 
subplot(1,2,1);axis equal; hold on;
plot(y_lims, y_lims, 'k');
plot(test_xy(:,2,1), y_p_all, 'r.');
axis([y_lims y_lims]);
title('Predicting pixel location (y) - all pixels');
xlabel('y-coordindate - ground truth');
ylabel('Predicted y');

subplot(1,2,2);axis equal; hold on;
plot(y_lims, y_lims, 'k');
plot(test_xy(:,2,2), y_p_comp, 'r.');
axis([x_lims x_lims]);
title('Predicting pixel location (y) - compressed region');
xlabel('y-coordindate - ground truth');
ylabel('Predicted y');
%
x_p_all = random_forest_reg_predict(rf_pos_x_all, test_responses(:,:,1), 0);
x_p_comp = random_forest_reg_predict(rf_pos_x_comp, test_responses(:,:,2), 0);
%
x_lims = [min(test_xy(:,1,1)) max(test_xy(:,1,1))];
figure; 
subplot(1,2,1);axis equal; hold on;
plot(x_lims, x_lims, 'k');
plot(test_xy(:,1,1), x_p_all, 'r.');
axis([x_lims x_lims]);
title('Predicting pixel location (x) - all pixels');
xlabel('x-coordindate - ground truth');
ylabel('Predicted x');

subplot(1,2,2);axis equal; hold on;
plot(x_lims, x_lims, 'k');
plot(test_xy(:,1,2), x_p_comp, 'r.');
axis([x_lims x_lims]);
title('Predicting pixel location (x) - compressed region');
xlabel('x-coordindate - ground truth');
ylabel('Predicted x');
%%
rf_dir = 'C:\isbe\density\assure\375_125\models\position_regression\norm_s\comp\';
rf_args.predictor_path = [rf_dir 'rf.mat'];
rf_args.tree_dir = [rf_dir 'trees\'];
rf_args.sampling_args.y = train_s(:,1,2);
rf_args.sampling_args.X = train_responses(:,:,2);
rf_pos_s_comp = random_forest_reg_train(rf_args);

% rf_dir = 'C:\isbe\density\assure\375_125\models\position_regression\norm_phi\comp\';
% rf_args.predictor_path = [rf_dir 'rf.mat'];
% rf_args.tree_dir = [rf_dir 'trees\'];
% rf_args.sampling_args.y = train_phi(:,1,2);
% rf_args.sampling_args.X = train_responses(:,:,2);
% rf_pos_phi_comp = random_forest_reg_train(rf_args);
%
rf_dir = 'C:\isbe\density\assure\375_125\models\position_regression\norm_phi_c\comp\';
rf_args.predictor_path = [rf_dir 'rf.mat'];
rf_args.tree_dir = [rf_dir 'trees\'];
rf_args.sampling_args.y = complex(cos(2*train_phi(:,1,2)), sin(2*train_phi(:,1,2)));
rf_args.sampling_args.X = train_responses(:,:,2);
rf_pos_phi_c_comp = random_forest_reg_train(rf_args);
%%
s_p_comp = random_forest_reg_predict(rf_pos_s_comp, test_responses(:,:,1), 0);
% phi_p_comp = random_forest_reg_predict(rf_pos_phi_comp, test_responses(:,:,2), 0);
phi_c_p_comp = random_forest_reg_predict(rf_pos_phi_c_comp, test_responses(:,:,2), 1);
phi_o_p_comp = angle(phi_c_p_comp)/2;
phi_o_p_comp(phi_o_p_comp <0) = phi_o_p_comp(phi_o_p_comp<0) + pi;
%
s_lims = [min(test_s(:,1,2)) max(test_s(:,1,2))];
phi_lims = [min(test_phi(:,1,2)) max(test_phi(:,1,2))];
figure; 
subplot(1,2,1);axis equal; hold on;
plot(s_lims, s_lims, 'k');
plot(test_s(:,1,2), s_p_comp, 'r.');
axis([s_lims s_lims]);
title('Predicting pixel location (s) - compressed region');
xlabel('Polar coordiates, distance to nipple - ground truth');
ylabel('Predicted s');

subplot(1,2,2);axis equal; hold on;
plot(phi_lims, phi_lims, 'k');
plot(test_phi(:,1,2), phi_o_p_comp, 'r.');
axis([phi_lims phi_lims]);
title('Predicting pixel location (\phi) - compressed region');
xlabel('Polar coordiates, angluar component - ground truth');
ylabel('Predicted \phi');
%
%%
decomposition_args.decomp_type = 'dt';      %Use the dual-tree
decomposition_args.win_size = 1;            %Window size about pixel to sample features
decomposition_args.levels = 1:6;            %Levels to include
decomposition_args.feature_shape = 'rect';  %Keep this as rect, although there are other things you can play with (see sample_dt_data)
decomposition_args.feature_type = 'conj';   %Again, 'conj' works best, but can use 'mag' for magnitude only, or 'phase' (see convert_complex_representation)
decomposition_args.do_max = 0;              %Keep all 6 sub-bands (0) or only keep the maximum (1)
decomposition_args.rotate = 0;              %Try and make rotation invariant (best left switched off)
decomposition_args.use_nag = 0;             %Use the NAG toolbox to interpolate (best left swicthed off)
decomposition_args.normalise = 0;

mkdir c:\isbe\density\assure\375_125\images
mkdir c:\isbe\density\assure\375_125\full_masks
mkdir c:\isbe\density\assure\375_125\compressed_masks
mkdir c:\isbe\density\assure\375_125\label_masks

for i_case = 1:num_cases
    
    mam_vdm = imread(case_im_names{i_case});
    mask = imread(case_mask_names{i_case});   
    
    [~, fname, ext] = fileparts(case_im_names{i_case});
    
%     dual_tree = compute_dual_tree(mam_vdm, 6);
% 
%     [rows cols] = find(mask == 4);
%     [responses] = compute_filter_responses(mam_vdm, decomposition_args);
%     samples_X = sample_image_features(responses, rows, cols, decomposition_args);
%     if i_case <= num_controls
%         samples_y = false(size(rows));
%     else
%         samples_y = true(size(rows));
%     end
%         
%     sample_name = [];
%     save(sample_name, 'samples_X', 'samples_y');
end
%%
im1 = imread(case_im_names{1});
mask1 = imread(case_mask_names{1});
mask1(mask1==4) = 2;

im2 = imread(case_im_names{2});
mask2 = imread(case_mask_names{2});
mask2(mask2==4) = 2;

points_mask = false(size(mask1));
points_mask(32:32:end,32:32:end) = 1;
points_mask(~(mask1==2)) = 0;
figure; imgray(points_mask);
%points_mask(381, 193) = 1;
%%
[bcs1, bcs2] = BreastCoodSystem.imageCorrespondentPoints(im1, im2, mask1, mask2, find(points_mask), 1, 0, 2, 3);
%%
    