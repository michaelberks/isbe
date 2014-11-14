%--------------------------------------------------------------------------
% Random stuff to do with finding the dge of vessels I thinks
%--------------------------------------------------------------------------


        
        
        
        
        
        
    
    

%%
%%
v_files = dir([contour_dir '*vessel_contour.mat']);
segment_length = 3;

contour_mat = zeros(1e4, 2*(segment_length+1));
widths_mat = zeros(1e4, segment_length+1);
v_counter = 1;
for i_v = 1:length(v_files)
    vessel_struc = load([contour_dir v_files(i_v).name]);  
    num_pts = size(vessel_struc.vessel_centre,1);
    for i_p = 1:num_pts-segment_length
        v_segment = vessel_struc.vessel_centre(i_p:i_p+segment_length,:);
        i_segment = vessel_struc.inner_edge(i_p:i_p+segment_length,:);
        o_segment = vessel_struc.outer_edge(i_p:i_p+segment_length,:);
        contour_mat(v_counter, :) = v_segment(:)';
        widths_mat(v_counter, :) = sqrt(sum((i_segment-o_segment).^2,2))';
        v_counter = v_counter + 1;
    end    
end
contour_mat(v_counter:end,:,:) = [];
widths_mat(v_counter:end,:,:) = [];
%%

known_dims_shape = [1:segment_length segment_length+1+(1:segment_length)];
con_dims_shape = [segment_length+1 2*(segment_length+1)];
known_dims_width = 1:segment_length;
con_dims_width = segment_length+1;
known_dims = 1:3*(segment_length);
con_dims = 3*(segment_length)+(1:3);

[~, a_scales, ~, a_rots, a_trans a_origins] = align_shapes(contour_mat(:,known_dims_shape), 'area', 0, 'length', 20);
num_shapes = size(contour_mat,1);
a_shapes = zeros(num_shapes, 2*(segment_length+1));
for i_sh = 1:num_shapes
    u_shape_i = reshape(contour_mat(i_sh,:), segment_length+1, 2);
    a_shape_i = bsxfun(@plus, u_shape_i*a_rots(:,:,i_sh), a_trans(i_sh,:)) *a_scales(i_sh);
    a_shapes(i_sh,:) = a_shape_i(:)';   
end
%%

[z_data, z_mu, z_sigma] = zscore([...
    a_shapes(:,known_dims_shape)...
    widths_mat(:,known_dims_width)...
    a_shapes(:, con_dims_shape)...
    widths_mat(:,con_dims_width)]);

[P_known, B_known, L_known] = princomp(z_data(:,known_dims));

keep_modes = find(cumsum(L_known)/sum(L_known) >= 0.98, 1);
P_known(:,keep_modes+1:end) = [];
B_known(:,keep_modes+1:end) = [];

x_data = [B_known z_data(:,con_dims)];
x_mu = mean(x_data);
x_sigma = cov(bsxfun(@minus, x_data, x_mu));

%%
v_pts = vessel_struc.vessel_centre;
i_pts = vessel_struc.inner_edge;
o_pts = vessel_struc.outer_edge;

widths = sqrt(sum((vessel_struc.inner_edge-vessel_struc.outer_edge).^2,2))';
num_pts = size(v_pts,1);
scores = zeros(num_pts,11);

for i_pt = 1:num_pts;
    for i_off = -5:5
        v_pts2 = v_pts;
        v_pts2(i_pt,1) = v_pts2(i_pt,1)+i_off;
        [log_p] = gaussian_chain_prob(v_pts2, widths, x_mu, x_sigma, P_known, z_mu, z_sigma, segment_length+1);
        scores(i_pt,i_off+6) = log_p;
    end
end
figure; imgray(scores');
%%
v_pts = vessel_struc.vessel_centre;
i_pts = vessel_struc.inner_edge;
o_pts = vessel_struc.outer_edge;

widths = sqrt(sum((vessel_struc.inner_edge-vessel_struc.outer_edge).^2,2))';
num_pts = size(v_pts,1);
scores = zeros(num_pts,11);

for i_pt = 1:num_pts;
    for i_off = -5:5
        widths2 = widths;
        widths2(i_pt) = widths2(i_pt)+2*i_off;
        [log_p] = gaussian_chain_prob(v_pts, widths2, x_mu, x_sigma, P_known, z_mu, z_sigma, segment_length+1);
        scores(i_pt,i_off+6) = log_p;
    end
end
figure; imgray(scores');

%%
figure; axis equal ij; hold all;
silly_states_xy = cell(num_pts,1);
silly_states_w = cell(num_pts,1);
for i_pt = 1:num_pts
    silly_states_xy{i_pt} = bsxfun(@plus, [0.5*randn(3,1); 0.05*randn]*(o_pts(i_pt,:)-v_pts(i_pt,:)), v_pts(i_pt,:));
    silly_states_xy{i_pt} = silly_states_xy{i_pt}(randperm(4),:);
    silly_states_w{i_pt} = ones(1,4)*widths(i_pt);
    plot(silly_states_xy{i_pt}(:,1), silly_states_xy{i_pt}(:,2), 'x');
end
%%
figure; hold all;
for i_pt = 1:num_pts
    silly_states_xy{i_pt} = repmat(v_pts(i_pt,:), 4,1);
    silly_states_w{i_pt} = [5*randn(1,3) 0.00*randn] + widths(i_pt);
    plot(i_pt, silly_states_w{i_pt}, 'x');
end
%%
[best_idx, min_E, best_xy, best_w] = dynamic_gaussian_chain(silly_states_xy, silly_states_w, x_mu, x_sigma, P_known, z_mu, z_sigma, segment_length+1);
plot(best_w, 'r-', 'linewidth', 2);
%%
v_files = dir([contour_dir '*vessel_contour.mat']);
i_v = 1;
vessel_struc = load([contour_dir v_files(i_v).name]);
apex_struc = u_load([vessel_dir v_files(i_v).name(1:end-12) '.mat']);

vessel_patch = double(apex_struc.vessel_patch);
g = gaussian_filters_1d(16, 48);
g = g / sum(g);
im_edges = conv2(g', g, ones(size(vessel_patch)), 'same');
vessel_patch_smoothed = conv2(g', g, vessel_patch, 'same') ./ im_edges;
vessel_patch_equalised = vessel_patch - vessel_patch_smoothed;

figure;
imgray(vessel_patch_equalised);
plot(vessel_struc.vessel_centre(:,1), vessel_struc.vessel_centre(:,2));
%%
[profile_features displacements] = extract_vessel_profile_features(vessel_patch_equalised,  vessel_struc.vessel_centre, vessel_struc.outer_edge, vessel_struc.inner_edge, ...
    'theta_offsets', 0,...pi*[-15 0 15]/180,...
    'scale_offsets', 1,...
    'num_displacements', 2);
%%
%clear;
contour_dir = 'C:\isbe\nailfold\data\rsa_study\vessel_contours\normal\';
vessel_dir = 'C:\isbe\nailfold\data\rsa_study\apexes\normal\';
[training_x, ~, training_y] = extract_vessel_profile_features_set(contour_dir, vessel_dir, ...
    'num_profile_pts', 50,...
    'theta_offsets', pi*[-15 0 15]/180,...
    'scale_offsets', [.5 1 2],...
    'num_displacements', 10,...
    'selected_idx', 1:10,...
    'plot', 0);
%
[test_x, ~, test_y] = extract_vessel_profile_features_set(contour_dir, vessel_dir, ...
    'num_profile_pts', 50,...
    'theta_offsets', pi*[-15 0 15]/180,...
    'scale_offsets', [.5 1 2],...
    'num_displacements', 10,...
    'selected_idx', 41:50,...
    'plot', 0);
%%
warning('off', 'ASYM:unexpectedArgument');

rf_dir = 'C:\isbe\nailfold\data\rsa_study\vessel_contours\rfs\';
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
%%
rf_args.tree_dir = [rf_dir 'apex_width_prediction\trees\'];
rf_args.sampling_args.y = test_y;
rf_args.sampling_args.X = test_x;
predictor = random_forest_reg_train(rf_args);
predictions = random_forest_reg_predict(predictor, training_x, 0);

n_pts = size(training_x,1);
figure; plot(sort(abs(predictions-training_y)), (1:n_pts)/n_pts);
%%
rf_args.tree_dir = [rf_dir 'bob\trees\'];
rf_args.sampling_args.y = training_y;
rf_args.sampling_args.X = training_x(:,1:45);
predictor = random_forest_reg_train(rf_args);
predictions = random_forest_reg_predict(predictor, test_x(:,1:45), 0);
%%
%%
rf_args.tree_dir = [rf_dir 'fred\trees\'];
rf_args.sampling_args.y = training_y;
rf_args.sampling_args.X = training_x(:,46:end);
predictor = random_forest_reg_train(rf_args);
predictions = random_forest_reg_predict(predictor, test_x(:,46:end), 0);
%%
joey = zeros(3,3,3);
for ii = 1:3
    for jj = 1:3
        for kk = 1:5
            joey(ii,jj,kk) = 100*ii + 10*jj + kk;
        end
    end
end
joey(:)
%%
selected_dims = false(3,3,5);
selected_dims(2,:,:) = 1;
selected_dims = [selected_dims(:)' selected_dims(:)'];

rf_args.tree_dir = [rf_dir 'flora\trees\'];
rf_args.sampling_args.y = training_y;
rf_args.sampling_args.X = training_x(:,selected_dims);
predictor = random_forest_reg_train(rf_args);
predictions = random_forest_reg_predict(predictor, test_x(:,selected_dims), 0);
%%
old_n_pts = size(training_x,1);
new_training_x = zeros(old_n_pts*3, size(training_x,2)/3);
new_training_y = zeros(old_n_pts*3, 1);
new_training_th = zeros(old_n_pts*3, 1);

theta_offsets = pi*[-15 0 15]/180;
selected_rows = 1:old_n_pts;
for i_th = 1:3
    selected_dims = false(3,3,5);
    selected_dims(:,i_th,:) = 1;
    selected_dims = [selected_dims(:)' selected_dims(:)'];
    
    
    new_training_x(selected_rows,:) = training_x(:,selected_dims);
    new_training_y(selected_rows,:) = training_y;
    new_training_th(selected_rows,:) = theta_offsets(i_th);
    
    selected_rows = selected_rows + old_n_pts;
end
%%
old_n_pts = size(test_x,1);
new_test_x = zeros(old_n_pts*3, size(test_x,2)/3);
new_test_y = zeros(old_n_pts*3, 1);
new_test_th = zeros(old_n_pts*3, 1);

theta_offsets = pi*[-15 0 15]/180;
selected_rows = 1:old_n_pts;
for i_th = 1:3
    selected_dims = false(3,3,5);
    selected_dims(:,i_th,:) = 1;
    selected_dims = [selected_dims(:)' selected_dims(:)'];
    
    
    new_test_x(selected_rows,:) = test_x(:,selected_dims);
    new_test_y(selected_rows,:) = test_y;
    new_test_th(selected_rows,:) = theta_offsets(i_th);
    
    selected_rows = selected_rows + old_n_pts;
end
%%
training_n = size(training_x,1);

rf_args.tree_dir = [rf_dir 'tim\trees\'];
rf_args.sampling_args.y = new_training_y(training_n+1:2*training_n,:);
rf_args.sampling_args.X = new_training_x(training_n+1:2*training_n,:);
predictor = random_forest_reg_train(rf_args);
predictions = random_forest_reg_predict(predictor, new_test_x(old_n_pts+1:2*old_n_pts,:), 0);
%%
contour_dir = 'C:\isbe\nailfold\data\rsa_study\vessel_contours\normal\';
vessel_dir = 'C:\isbe\nailfold\data\rsa_study\apexes\normal\';
v_files = dir([contour_dir '*vessel_contour.mat']);
i_v = 52;
vessel_struc = load([contour_dir v_files(i_v).name]);
apex_struc = u_load([vessel_dir v_files(i_v).name(1:end-12) '.mat']);    

vessel_patch = double(apex_struc.vessel_patch);
g = gaussian_filters_1d(16, 48);
g = g / sum(g);
im_edges = conv2(g', g, ones(size(vessel_patch)), 'same');
vessel_patch_smoothed = conv2(g', g, vessel_patch, 'same') ./ im_edges;
vessel_patch_equalised = vessel_patch - vessel_patch_smoothed;

vessel_pts = spline_contour(apex_struc.v_pts, [], 4);
%%
predictor.trees = cell(100,1);
for i_tree = 1:100
    predictor.trees{i_tree} = ['rf_tree' zerostr(i_tree,4) '.mat'];
end
predictor.D = 30;
predictor.d = 5;
predictor.tree_dir = 'C:\isbe\nailfold\data\rsa_study\vessel_contours\rfs\bob\greg\trees\';
predictor.tree_root = [];
predictor.regression_method = 'rf_regression';
predictor.sampled_data_dir = [];
     
%%
profile on;
[prediction_map all_dis all_pred] = predict_vessel_profile_features_image(vessel_patch_equalised, ...
    vessel_struc.vessel_centre, 12, predictor, 'plot', 1);
profile viewer;
%%
for i_v = 41:52
    vessel_struc = load([contour_dir v_files(i_v).name]);
    apex_struc = u_load([vessel_dir v_files(i_v).name(1:end-12) '.mat']);    

    vessel_patch = double(apex_struc.vessel_patch);
    g = gaussian_filters_1d(16, 48);
    g = g / sum(g);
    im_edges = conv2(g', g, ones(size(vessel_patch)), 'same');
    vessel_patch_smoothed = conv2(g', g, vessel_patch, 'same') ./ im_edges;
    vessel_patch_equalised = vessel_patch - vessel_patch_smoothed;

    vessel_pts = spline_contour(apex_struc.v_pts, [], 4);
    
    [prediction_map] = predict_vessel_profile_features_image(vessel_patch_equalised, ...
        vessel_pts(1:2:end,:), 12, predictor, 'plot', 1);
end
%%
for i_v = 11:16%:20
    vessel_struc = load([contour_dir v_files(i_v).name]);
    apex_struc = u_load([vessel_dir v_files(i_v).name(1:end-12) '.mat']);    
    apex_struc2 = load([vessel_dir v_files(i_v).name(1:end-19) '.mat']);
    vessel_patch = double(apex_struc.vessel_patch);
    g = gaussian_filters_1d(16, 48);
    g = g / sum(g);
    im_edges = conv2(g', g, ones(size(vessel_patch)), 'same');
    vessel_patch_smoothed = conv2(g', g, vessel_patch, 'same') ./ im_edges;
    vessel_patch_equalised = vessel_patch - vessel_patch_smoothed;

    vessel_pts = spline_contour(apex_struc.v_pts, [], 4);

    apex_width = sqrt(sum(diff(apex_struc2.apex_xy).^2,2));
    [prediction_map vessel_centre_states vessel_centre_scores] = predict_vessel_profile_features_image(vessel_patch_equalised, ...
            vessel_pts, apex_width, predictor, 'plot', 1, 'debug', 0);
    num_pts = length(vessel_centre_states);

    vessel_width_states = cell(num_pts, 1);  
    figure; imgray(vessel_patch_equalised);
    for i_pt = 1:num_pts
        plot(vessel_centre_states{i_pt}(:,1), vessel_centre_states{i_pt}(:,2), 'x');
        vessel_width_states{i_pt} = 12*ones(size(vessel_centre_states{i_pt},1),1);
    end

    [best_idx, min_E, best_xy, best_w] = ...
        dynamic_gaussian_chain(vessel_centre_states, vessel_width_states,...
            x_mu, x_sigma, P_known, z_mu, z_sigma, segment_length+1);

    plot(best_xy(:,1), best_xy(:,2), 'c');
end
    


    