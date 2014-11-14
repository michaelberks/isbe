apex_radius = 5;
max_apex_guess = 500;
vessel_prob_smoothing_sigma = 2;
curvature_smoothing_sigma = 2;
strong_vessel_thresh = 0.25;
curv_max = 0.5;

g = gaussian_filters_1d(vessel_prob_smoothing_sigma);
g = g / sum(g);

rsa_dir = 'rsa_study/';

cell_sz = 8;
num_cells = 8;
patch_sz = num_cells*cell_sz;
patch_sz = patch_sz + 2; %Account for padding
patch_sz2 = (patch_sz - 1)/2;

x = repmat(-patch_sz2:patch_sz2, patch_sz, 1);
y = x';
xy = [x(:) y(:)];
dist_thresh = 24^2;
base_width = 20;

apex_dir = [nailfoldroot 'data/' rsa_dir '/vessel_contours/'];
model_dir = [nailfoldroot 'data/' rsa_dir 'models/apex_templates/'];
    
%
for i_test = {'set12', 'set1g'}
    test_dir = i_test{1};
    
    image_dir = [nailfoldroot 'data/' rsa_dir test_dir '/images/'];
    fov_mask_dir = [nailfoldroot 'data/' rsa_dir test_dir '/fov_masks/'];
    vessel_mask_dir = [nailfoldroot 'data/' rsa_dir test_dir '/vessel_masks/'];
    prob_dir = [nailfoldroot 'data/' rsa_dir test_dir '/predictions/detection/rf_classification/257273/'];
    ori_dir = [nailfoldroot 'data/' rsa_dir test_dir '/predictions/orientation/rf_regression/259076/'];
    width_dir = [nailfoldroot 'data/' rsa_dir test_dir '/predictions/width/rf_regression/257847/'];
    apex_class_dir = [nailfoldroot 'data/' rsa_dir test_dir '/apex_class_data/'];
    create_folder(apex_class_dir);
    
    pred_list = dir([prob_dir '*.mat']);
    num_ims = length(pred_list);
    for i_im = 1:num_ims

        display(['Processing image ' num2str(i_im) ' of ' num2str(num_ims)]);

        vessel_im = u_load([image_dir pred_list(i_im).name(1:end-9) '.mat']);
        vessel_prob = u_load([prob_dir pred_list(i_im).name]);
        vessel_ori = u_load([ori_dir pred_list(i_im).name]);
        width_map = u_load([width_dir pred_list(i_im).name]);
        vessel_mask = u_load([vessel_mask_dir pred_list(i_im).name(1:end-9) '_v_mask.mat']);

        vessel_prob = conv2(g', g, vessel_prob, 'same');
        width_map = conv2(g', g, width_map, 'same');
        vessel_ori = conv2(g', g, vessel_ori, 'same');

        %Compute NMS centrelines
        vessel_nms = mb_non_maximal_supp(vessel_prob, angle(vessel_ori)/2);
        strong_vessels = vessel_nms > strong_vessel_thresh;
        if any(strong_vessels(:))
            [rstrong cstrong] = find(strong_vessels);
            vessel_centre_mask = bwselect(vessel_nms > 0, cstrong, rstrong, 8);
        else
            vessel_centre_mask = strong_vessels;
        end
        [vessel_centre_y vessel_centre_x] = find(vessel_centre_mask);   
        num_pts = length(vessel_centre_y);
        [n_rows n_cols] = size(vessel_im);


        %Load in the vessel apex
        if pred_list(i_im).name(1) == 'e'
            load([apex_dir 'enlarged/' pred_list(i_im).name(9:16) '_vessel_contour.mat'], 'apex_idx', 'vessel_centre');
        elseif pred_list(i_im).name(1) == 'n'
            load([apex_dir 'normal/' pred_list(i_im).name(7:14) '_vessel_contour.mat'], 'apex_idx', 'vessel_centre');
        elseif pred_list(i_im).name(1) == 'g'
            load([apex_dir 'giant/' pred_list(i_im).name(6:13) '_vessel_contour.mat'], 'apex_idx', 'vessel_centre');
        end

        %Correct apex coordinates frame
        apex_xy = [vessel_centre(apex_idx,1) vessel_centre(apex_idx,2)];

        apex_offsets = zeros(num_pts,2);
        apex_class = false(num_pts,1);
        apex_hog = zeros(num_pts, num_cells*num_cells*9);
        apex_oris = zeros(num_pts,1);
        apex_scales = zeros(num_pts,1);

        for i_pt = 1:length(vessel_centre_x)

            %Get predicted scale and orientation at this point
            vxc = vessel_centre_x(i_pt);
            vyc = vessel_centre_y(i_pt);
            ori_c = angle(vessel_ori(vyc, vxc))/2;
            width_c = width_map(vyc, vxc);

            %Get scale relative to base width a make rotation matrix
            rot = [cos(ori_c) -sin(ori_c); sin(ori_c) cos(ori_c)];
            scale = width_c / base_width;

            apex_oris(i_pt) = ori_c;
            apex_scales(i_pt) = width_c;

            %Transform points given scale and angle and translate to
            %candidate position
            xya = xy * rot * scale;
            xa = reshape(xya(:,1) + vxc, patch_sz, patch_sz);
            ya = reshape(xya(:,2) + vyc, patch_sz, patch_sz);

            apex_dist = inf;
            for i_ap = 1:length(apex_idx)
                %Transform vessel centre
                apex_xy_t = (apex_xy(i_ap,:) - [vxc vyc])*rot'/scale;
                apex_dist_i = (sum(apex_xy_t.^2));

                if apex_dist_i < apex_dist
                    apex_offsets(i_pt,:) = apex_xy_t;
                    apex_dist = apex_dist_i;
                end
            end

            %Sample vessel prob patch
            vessel_prob_patch = interp2(vessel_prob, xa, ya, '*linear', 0);
            [hog] = compute_HoG(vessel_prob_patch);       
            apex_hog(i_pt,:) = hog(:)';

            %Display patch if apex_xy_t lies within some dist
            if  vessel_mask(vyc, vxc) && (apex_dist < dist_thresh)  
                apex_class(i_pt) = 1;
            end
        end
        save([apex_class_dir pred_list(i_im).name(1:end-9) '_t.mat'], 'apex_scales', 'apex_oris', 'vessel_centre_x', 'vessel_centre_y', 'n_rows', 'n_cols');
        save([apex_class_dir pred_list(i_im).name(1:end-9) '_y.mat'], 'apex_class', 'apex_offsets');
        save([apex_class_dir pred_list(i_im).name(1:end-9) '_X.mat'], 'apex_hog');
    end
end
%
data_paths = cell(0,1);
data_counts = zeros(0,2);
im_num = 1;
for i_test = {'set12', 'set1g'}
    test_dir = i_test{1};
    apex_class_dir = [nailfoldroot 'data/' rsa_dir test_dir '/apex_class_data/'];
    
    data_list = dir([apex_class_dir '*_y.mat']);
    num_ims = length(data_list);
    
    for i_im = 1:num_ims
        data_paths{im_num,1} = [apex_class_dir data_list(i_im).name];
        
        load(data_paths{im_num,1}, 'apex_class');
        data_counts(im_num,1) = sum(apex_class);
        data_counts(im_num,2) = sum(~apex_class);
        im_num = im_num + 1;
    end
end
display(sum(data_counts));

%
    
total_ims = length(data_paths);

train_ims = rand(total_ims,1) < 0.5;
    
train_n = sum(data_counts(train_ims,1));

train_X = zeros(2*train_n, num_cells*num_cells*9);
train_c = [true(train_n,1); false(train_n, 1)];
train_y = zeros(train_n, 2);

curr_idx = 0;
for i_im = 1:total_ims
    
    if train_ims(i_im)
        load(data_paths{i_im,1}, 'apex_class', 'apex_offsets');
        load([data_paths{i_im,1}(1:end-5) 'X.mat'], 'apex_hog');
    
        idx = curr_idx + (1:data_counts(i_im,1));
        curr_idx = idx(end);
        
        n_idx = find(~apex_class);
        r_idx = randperm(data_counts(i_im,2));
        n_idx = n_idx(r_idx(1:data_counts(i_im,1)));
        
        train_X(idx,:) = apex_hog(apex_class,:);
        train_X(idx + train_n,:) = apex_hog(n_idx,:);
        
        train_y(idx,:) = apex_offsets(apex_class,:);
    end
end
%   
warning('off', 'ASYM:unexpectedArgument');
rf_args.prediction_type = 'rf_classification';
rf_args.n_trees = 100;
rf_args.d = [];
rf_args.w_prior = 0;
rf_args.impure_thresh = 1.0000e-004;
rf_args.split_min = 100;
rf_args.end_cut_min = 25;
rf_args.do_ubound = 0;
rf_args.quiet = 1;
rf_args.overwrite = 0;
rf_args.minimise_size = 0;
rf_args.split_criterion = 'gdi';
rf_args.var_criterion = 'mabs';

rf_args.sampling_args.sampling_method = 'sample_saved_training_data';
rf_args.decomposition_args = [];

rf_dir = ['C:\isbe\nailfold\models\vessel\apex_location\class' datestr(now, 30) '\'];
rf_args.tree_dir = [rf_dir 'trees/'];            

rf_args.sampling_args.y = train_c;
rf_args.sampling_args.X = train_X;
apex_class_rf = random_forest_class_train(rf_args);   
    
%  
rf_args.prediction_type = 'rf_regression';
rf_args.n_trees = 100;
rf_args.d = [];
rf_args.w_prior = 0;
rf_args.impure_thresh = 1.0000e-008;
rf_args.split_min = 100;
rf_args.end_cut_min = 0;
rf_args.do_ubound = 0;
rf_args.quiet = 1;
rf_args.do_circular = [];
rf_args.overwrite = 1;
rf_args.minimise_size = 0;
rf_args.split_criterion = 'ssq';
rf_args.var_criterion = 'ssq';

rf_args.sampling_args.sampling_method = 'sample_saved_training_data';
rf_args.decomposition_args = [];

rf_args.sampling_args.X = train_X(1:train_n,:);

rf_dir = ['C:\isbe\nailfold\models\vessel\apex_location\offset_x\' datestr(now, 30) '\'];
rf_args.tree_dir = [rf_dir 'trees/'];            
rf_args.sampling_args.y = train_y(:,1);
apex_offset_x_rf = random_forest_reg_train(rf_args); 

rf_dir = ['C:\isbe\nailfold\models\vessel\apex_location\offset_y\' datestr(now, 30) '\'];
rf_args.tree_dir = [rf_dir 'trees/'];            
rf_args.sampling_args.y = train_y(:,2);
apex_offset_y_rf = random_forest_reg_train(rf_args);

%
test_n = sum(data_counts(~train_ims,:));
all_class_pred_1 = zeros(test_n(1), 1);
all_class_pred_0 = zeros(test_n(2), 1);
all_offset_pred = zeros(test_n(1),2);

curr_idx_1 = 0;
curr_idx_0 = 0;
for i_im = 1:total_ims
    
    if ~train_ims(i_im)
        load(data_paths{i_im,1}, 'apex_class', 'apex_offsets');
        load([data_paths{i_im,1}(1:end-5) 'X.mat'], 'apex_hog');
    
        idx_1 = curr_idx_1 + (1:data_counts(i_im,1));
        curr_idx_1 = idx_1(end);
        
        idx_0 = curr_idx_0 + (1:data_counts(i_im,2));
        curr_idx_0 = idx_0(end);
        
        %do prediction
        apex_offset_x_pred = random_forest_reg_predict(apex_offset_x_rf, apex_hog);
        apex_offset_y_pred = random_forest_reg_predict(apex_offset_y_rf, apex_hog);
        [~,apex_class_pred] = random_forest_class_predict(apex_class_rf, apex_hog);
        apex_class_pred = apex_class_pred / rf_args.n_trees;
        
        %Save in main containers
        all_class_pred_1(idx_1,:) = apex_class_pred(apex_class);
        all_class_pred_0(idx_0,:) = apex_class_pred(~apex_class);
        all_offset_pred(idx_1,1) = apex_offset_x_pred(apex_class,:);
        all_offset_pred(idx_1,1) = apex_offset_y_pred(apex_class,:);
        
        %Save
        save([data_paths{i_im,1}(1:end-5) 'pred.mat'], 'apex_offset_*_pred', 'apex_class_pred');
        
    end
end
%
[roc_pts, auc, ~, ~, auc_se] = calculate_roc_curve(1-[all_class_pred_1; all_class_pred_0],[true(test_n,1); false(test_n,1)]);
        
figure;
axis([0 1 0 1]); axis equal; hold all;
plot(roc_pts(:,1), roc_pts(:,2), '-', 'linewidth', 2); 
title(['AUC = ' num2str(100*auc,3) ' \pm ' num2str(100*auc_se,3)])
%%
for i_im = [142 201]
    
    if ~train_ims(i_im)
        
        [~,apex_name] = fileparts(data_paths{i_im,1}(1:end-6));
        vessel_im = u_load([image_dir apex_name '.mat']);
        
        load([data_paths{i_im,1}(1:end-5) '_pred.mat'], 'apex_class_pred', 'apex_offset_*_pred');
        load([data_paths{i_im,1}(1:end-5) 't.mat']);
        load([data_paths{i_im,1}(1:end-5) 'y.mat']);
        
        vessel_centre_idx = sub2ind([n_rows n_cols], vessel_centre_y, vessel_centre_x);
        
        apex_class_map = zeros(n_rows, n_cols);
        apex_class_map(vessel_centre_idx) =  apex_class_pred(:,2);
        
        apex_offset_map = zeros(n_rows, n_cols);
        
        for i_pt = 1:length(vessel_centre_x)
            
            if apex_class_pred(i_pt,2) > 0%.5
                vxc = vessel_centre_x(i_pt);
                vyc = vessel_centre_y(i_pt);
                ori_c = apex_oris(i_pt);
                scale = apex_scales(i_pt) / base_width;

                %Get scale relative to base width a make rotation matrix
                rot = [cos(ori_c) -sin(ori_c); sin(ori_c) cos(ori_c)];

                apex_xy_t = [apex_offset_x_pred(i_pt) apex_offset_y_pred(i_pt)];
                apex_xy_p = scale*apex_xy_t*rot + [vxc vyc];           
               
                if all(apex_xy_p > 0.5) && apex_xy_p(1) < n_cols && apex_xy_p(2) < n_rows
                    
                    [im_r, im_c, gxy] = add_gaussian_blur(apex_xy_p(1), apex_xy_p(2), 2*scale, n_rows, n_cols);
                    
                    apex_offset_map(im_r, im_c) = apex_offset_map(im_r, im_c) + apex_class_pred(i_pt,2)*gxy;
                    
                    %apex_offset_map(round(apex_xy_p(2)), round(apex_xy_p(1))) = ...
                    %    apex_offset_map(round(apex_xy_p(2)), round(apex_xy_p(1))) + apex_class_pred(i_pt,2);
                end
            end
        end
        %apex_offset_map = conv2(g', g, apex_offset_map, 'same');
        
        figure; 
        subplot(1,3,1); imgray(vessel_im);
        subplot(1,3,2); imgray(apex_class_map); %title(num2str(i_im));
        subplot(1,3,3); imgray(apex_offset_map);
        
        
    end
end
%%
        
        