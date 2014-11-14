build_rf_line_detector(...
    'job_id', 'reg_test',...
    'detection_type', 'orientation',...
    'sampling_method','generate_line_training_data',...
    'num_samples', 2e3,...
    'pts_per_image', 50,...
    'use_nag', 1,...
    'width_range', [4 5],...
    'contrast_range', [4 8],...
    'decay_rate', 4,...
    'decomp_type', 'dt',...
    'num_levels', 3,...
    'win_size', 1);

%%
v1 = zeros(180,1);
v2 = zeros(180,1);

for ii = 1:180
    x = (181-ii)*rand(2e5,1);
    
    v1(ii) = var(x);
    v2(ii) = 1 - norm([mean(cosd(2*x)) mean(sind(2*x))]);
end

figure; plot(v1, v2, 'b-x');
xlabel('Linear variance');
ylabel('Circular variance');
figure; plot(181 - (1:180), v2, 'b-x');
xlabel('Range in degrees');
ylabel('Circular variance');
%%
v1 = zeros(100,1);
v2 = zeros(100,1);
for ii = 1:100
    x = (5*ii/500)*rand(2e5,1);
    
    v1(ii) = var(x);
    v2(ii) = 1 - norm([mean(cosd(2*x)) mean(sind(2*x))]);
end
figure; plot(v1, v2, 'b-x');
xlabel('Linear variance');
ylabel('Circular variance');
figure; plot(5*(1:100)/100, v2, 'b-x');
xlabel('Range in degrees');
ylabel('Circular variance');
%%
[X y] = generate_line_training_data(...
	'num_samples', 1e5,...
    'bg_dir', 'C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\',... % the mandatory arguments
    'decomp_type', 'dt',...
    'num_bgs', 200,...
    'bg_stem', 'bg',...
    'bg_zeros', 5,...
    'detection_type', 'orientation',...
    'pts_per_image', 250,...
    'bg_ratio', 1,...
    'width_range', [4 16],...
    'orientation_range', [0 360],...
    'contrast_range', [4 8],...
    'decay_rate', 4,...
    'line_type', 'curve',...
    'normalise', 0,...
    'num_levels', 3,...
    'feature_shape', 'rect',...
    'feature_type', 'all',...
    'do_max', 0,...
    'rotate', 0,...
    'win_size', 1);
save C:\isbe\asymmetry_project\data\misc\reg_data1e5.mat X y
%%
[Xt yt] = generate_line_training_data(...
	'num_samples', 1e3,...
    'bg_dir', 'C:\isbe\asymmetry_project\data\synthetic_backgrounds\real512\train\',... % the mandatory arguments
    'decomp_type', 'dt',...
    'num_bgs', 840,...
    'bg_stem', 'bg',...
    'bg_zeros', 5,...
    'detection_type', 'orientation',...
    'pts_per_image', 2,...
    'bg_ratio', 1,...
    'width_range', [4 16],...
    'orientation_range', [0 360],...
    'contrast_range', [4 8],...
    'decay_rate', 4,...
    'line_type', 'curve',...
    'normalise', 0,...
    'num_levels', 3,...
    'feature_shape', 'rect',...
    'feature_type', 'all',...
    'do_max', 0,...
    'rotate', 0,...
    'win_size', 1);
save C:\isbe\asymmetry_project\data\misc\reg_data_test.mat Xt yt
%%
[Xi yi] = generate_line_training_data(...
	'num_samples', 1e5,...
    'bg_dir', 'C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\',... % the mandatory arguments
    'decomp_type', 'dt',...
    'num_bgs', 200,...
    'bg_stem', 'bg',...
    'bg_zeros', 5,...
    'detection_type', 'orientation',...
    'pts_per_image', 250,...
    'bg_ratio', 1,...
    'width_range', [4 16],...
    'orientation_range', [0 360],...
    'contrast_range', [4 8],...
    'decay_rate', 4,...
    'line_type', 'curve',...
    'normalise', 0,...
    'num_levels', 3,...
    'feature_shape', 'rect',...
    'feature_type', 'ilp',...
    'do_max', 0,...
    'rotate', 0,...
    'win_size', 1);
save C:\isbe\asymmetry_project\data\misc\reg_data_ilp.mat Xi yi
%%
tic;
tree1 = mb_tree_rego_train(X, y);
toc;
save C:\isbe\asymmetry_project\data\misc\reg_tree1.mat tree1;
clear tree1;
%%
tic;
tree4 = mb_tree_rego_train(X, complex(cosd(2*y), sind(2*y)), 'mod', 1, 'impure_thresh', 1e-6);
toc;
save C:\isbe\asymmetry_project\data\misc\reg_tree4.mat tree4;
clear tree2;
%%
profile on;
tree3 = mb_tree_rego_train(X, y);
profile viewer;
%%
profile on;
tree4 = mb_tree_rego_train(X, complex(cosd(2*y), sind(2*y)), 'mod', 1, 'impure_thresh', 1e-4);
%tree4 = mb_tree_rego_train(X, [cosd(2*y), sind(2*y)], 'mod', 1, 'impure_thresh', 1e-4);
profile viewer;
%%
p_vals = zeros(100,1);
for ii = 1:100
    p_vals(ii) = cdf('bino', ii, 3*ii, 0.5);
end
%%
q_vals = zeros(100,1);
for ii = 1:100
    q_vals(ii) = cdf('bino', ii, 2*ii, 2/3);
end
%%
angles_count2 = zeros(180,1);

for ii = 1:200
    tree = u_load([random_forest.tree_root random_forest.tree_dir random_forest.trees{ii}]);
    
    leaves = find(~any(tree.children,2));
    for jj = 1:length(leaves)
        ori = ceil(tree.class(leaves(jj)));
        angles_count2(ori) = angles_count2(ori) + 1;
    end
end
%%
angles_count5 = zeros(180,1);
for kk = [1 2 3 6 7 8 9]
    load(['Z:\data\line_orientation_rfs\243311\random_forest' zerostr(kk,2) '.mat']);
    random_forest.tree_root = 'Z:\data\line_orientation_rfs\';
    
    for ii = 1:20
        tree = u_load([random_forest.tree_root random_forest.tree_dir random_forest.trees{ii}]);

        leaves = find(~any(tree.children,2));
        for jj = 1:length(leaves)
            ori = ceil(mod(90*angle(tree.class(leaves(jj)))/pi,180));
            angles_count5(ori) = angles_count5(ori) + 1;
        end
    end
end
%%
f1 = figure; hold all;
plot(1:180, angles_count / sum(angles_count), 'linewidth', 2);
plot(1:180, angles_count2 / sum(angles_count2), 'linewidth', 2);
plot(1:180, angles_count3 / sum(angles_count3), 'linewidth', 2);
plot(1:180, angles_count4 / sum(angles_count4), 'linewidth', 2);
plot(1:180, angles_count5 / sum(angles_count5), 'linewidth', 2);
title('Histogram of leaf node outputs for various random forests');
xlabel('Angle (degrees)');

legend({'DT-CWT, syn BGs (213209)', 'DT-CWT, real BGs (238470)', 'Linop (191959)', 'G2D (233142)', 'DT-CWT, modified regression'});
print_pdf('Z:\data\misc\RF_leaf_ouputs.pdf', f1);
%%
dt_dims_counts = zeros(432,1);
dt_dims_scores = zeros(432,1);

for jj = 1:200
    load(['C:\isbe\asymmetry_project\data\line_orientation_rfs\243311\rf_tree' zerostr(jj,4) '.mat']);
    leaves = find(~any(tree.children,2));
    for ii = 1:length(leaves)
        curr_node = leaves(ii);
        score = abs(tree.class(curr_node));

        while(curr_node ~= 1)
            curr_node = tree.parent(curr_node);
            split_dim = tree.var(curr_node);
            dt_dims_counts(split_dim) = dt_dims_counts(split_dim) + 1;
            dt_dims_scores(split_dim) = dt_dims_scores(split_dim) + score;
        end
    end
end
%
cmap = jet(9);
%%
f1 = figure; hold on; a1 = gca;
f2 = figure; hold on; a2 = gca;
f3 = figure; hold on; a3 = gca;
for pos = 1:9
    idx = (0:9:431) + pos;
    bar(a1, idx, dt_dims_counts(idx), 'facecolor', cmap(pos,:), 'barwidth', 0.1);
    bar(a2, idx, dt_dims_scores(idx), 'facecolor', cmap(pos,:), 'barwidth', 0.1);
    bar(a3, idx, dt_dims_scores(idx)./dt_dims_counts(idx), 'facecolor', cmap(pos,:), 'barwidth', 0.1);
end
%
y1 = get(a1, 'ylim');
y2 = get(a2, 'ylim');
y3 = get(a3, 'ylim');

plot(a1, [(54.5:54:432)' (54.5:54:432)'], 0.95*y1, 'k:')
plot(a2, [(54.5:54:432)' (54.5:54:432)'], 0.95*y2, 'k:')
plot(a3, [(54.5:54:432)' (54.5:54:432)'], 0.95*y3, 'k:')

plot(a1, [216.5 216.5], y1, 'r:')
plot(a2, [216.5 216.5], y2, 'r:')
plot(a3, [216.5 216.5], y3, 'r:')

set(a1, 'xtick', 27:54.5:432, 'xticklabel', {'Level 1', 'Level 2', 'Level 3', 'Level 4'});
axes(a1); text(100, y1(2), '\leftarrow Magnitude \rightarrow'); text(316, y1(2), '\leftarrow Phase \rightarrow');
title(a1, 'Number of times each dimension is used as a splitting variable');
set(a2, 'xtick', 27:54.5:432, 'xticklabel', {'Level 1', 'Level 2', 'Level 3', 'Level 4'});
axes(a2); text(100, y2(2), '\leftarrow Magnitude \rightarrow'); text(316, y2(2), '\leftarrow Phase \rightarrow');
title(a2, 'Sum of leaf dispersions for each dimension');
set(a3, 'xtick', 27:54.5:432, 'xticklabel', {'Level 1', 'Level 2', 'Level 3', 'Level 4'});
axes(a3); text(100, y3(2), '\leftarrow Magnitude \rightarrow'); text(316, y3(2), '\leftarrow Phase \rightarrow');
title(a3, 'Average leaf dispersion for each dimension');
%%
print_pdf('Z:\data\misc\dim_counts.pdf', f1);
print_pdf('Z:\data\misc\dim_scores.pdf', f2);
print_pdf('Z:\data\misc\dim_avg_scores.pdf', f3);
%%
load C:\isbe\asymmetry_project\data\misc\reg_data1e5.mat X y
[tree tree_nodes] = mb_tree_rego_train(X, complex(cosd(2*y), sind(2*y)), 'mod', 1, 'impure_thresh', 1e-6);
%%
leaves = find(~tree.var);
leaf_values = tree.class(leaves);
[leaf_oris leaves_idx] = sort(mod(90*angle(leaf_values)/pi,180), 'ascend');
%leaf_mags = abs(leaf_values(sort_idx));
num_leaves = length(leaves);


figure; hold on;
%
for ii = 1:num_leaves
    sample_oris = y(tree_nodes == leaves(leaves_idx(ii)));
    sample_oris = leaf_oris(ii) + mb_mod(sample_oris - leaf_oris(ii), 180);
    plot(ii, sample_oris, 'r.');
end
plot(1:num_leaves, leaf_oris);
plot(1:num_leaves, leaf_oris + 60, 'b--');
plot(1:num_leaves, leaf_oris - 60, 'b--');
plot(1:num_leaves, leaf_oris + 30, 'b--');
plot(1:num_leaves, leaf_oris - 30, 'b--');
plot(1:num_leaves, leaf_oris + 15, 'b--');
plot(1:num_leaves, leaf_oris - 15, 'b--');
%%
%load C:\isbe\asymmetry_project\data\misc\reg_data1e5.mat X y
[tree2 tree_nodes2] = mb_tree_reg_train(X, y, 'mod', 180);
%%
leaves = find(~tree2.var);
leaf_values = tree2.class(leaves);
[leaf_oris leaves_idx] = sort(mod(leaf_values,180), 'ascend');
%leaf_mags = abs(leaf_values(leaves_idx));
num_leaves = length(leaves);

%
figure; hold on;
%
for ii = 1:num_leaves
    sample_oris = y(tree_nodes2 == leaves(leaves_idx(ii)));
    if ~isempty(sample_oris)
        sample_oris = leaf_oris(ii) + mb_mod(sample_oris - leaf_oris(ii), 180);
        plot(ii, sample_oris, 'r.');
    end
end
plot(1:num_leaves, leaf_oris);
plot(1:num_leaves, leaf_oris + 60, 'b--');
plot(1:num_leaves, leaf_oris - 60, 'b--');
plot(1:num_leaves, leaf_oris + 30, 'b--');
plot(1:num_leaves, leaf_oris - 30, 'b--');
plot(1:num_leaves, leaf_oris + 15, 'b--');
plot(1:num_leaves, leaf_oris - 15, 'b--');
%%
dt_dims_counts = zeros(60,1);
dt_dims_scores = zeros(60,1);

for jj = 1:200
    load(['C:\isbe\asymmetry_project\data\line_detection_rfs\191905\rf_tree' zerostr(jj,4) '.mat']);
    leaves = find(~any(tree.children,2));
    for ii = 1:length(leaves)
        curr_node = leaves(ii);
        score = max(tree.classprob(curr_node,:));

        while(curr_node ~= 1)
            curr_node = tree.parent(curr_node);
            split_dim = tree.var(curr_node);
            dt_dims_counts(split_dim) = dt_dims_counts(split_dim) + 1;
            dt_dims_scores(split_dim) = dt_dims_scores(split_dim) + score;
        end
    end
end
%
cmap = jet(6);
%%
f1 = figure; hold on; a1 = gca;
f2 = figure; hold on; a2 = gca;
f3 = figure; hold on; a3 = gca;
for ori = 1:6
    idx = (0:6:59) + ori;
    bar(a1, idx, dt_dims_counts(idx), 'facecolor', cmap(ori,:), 'barwidth', 1/7);
    bar(a2, idx, dt_dims_scores(idx), 'facecolor', cmap(ori,:), 'barwidth', 1/7);
    bar(a3, idx, dt_dims_scores(idx)./dt_dims_counts(idx), 'facecolor', cmap(ori,:), 'barwidth', 1/7);
end
%%
y1 = get(a1, 'ylim');
y2 = get(a2, 'ylim');
y3 = get(a3, 'ylim');

plot(a1, [(54.5:54:432)' (54.5:54:432)'], 0.95*y1, 'k:')
plot(a2, [(54.5:54:432)' (54.5:54:432)'], 0.95*y2, 'k:')
plot(a3, [(54.5:54:432)' (54.5:54:432)'], 0.95*y3, 'k:')

plot(a1, [216.5 216.5], y1, 'r:')
plot(a2, [216.5 216.5], y2, 'r:')
plot(a3, [216.5 216.5], y3, 'r:')

set(a1, 'xtick', 27:54.5:432, 'xticklabel', {'Level 1', 'Level 2', 'Level 3', 'Level 4'});
axes(a1); text(100, y1(2), '\leftarrow Magnitude \rightarrow'); text(316, y1(2), '\leftarrow Phase \rightarrow');
title(a1, 'Number of times each dimension is used as a splitting variable');
set(a2, 'xtick', 27:54.5:432, 'xticklabel', {'Level 1', 'Level 2', 'Level 3', 'Level 4'});
axes(a2); text(100, y2(2), '\leftarrow Magnitude \rightarrow'); text(316, y2(2), '\leftarrow Phase \rightarrow');
title(a2, 'Sum of leaf dispersions for each dimension');
set(a3, 'xtick', 27:54.5:432, 'xticklabel', {'Level 1', 'Level 2', 'Level 3', 'Level 4'});
axes(a3); text(100, y3(2), '\leftarrow Magnitude \rightarrow'); text(316, y3(2), '\leftarrow Phase \rightarrow');
title(a3, 'Average leaf dispersion for each dimension');