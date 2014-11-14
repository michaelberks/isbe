
for ii = 0:15
    circ = 2^ii;
    
    N = round(1e5 / circ);
    theta = sort(2*pi*rand(N,1) / circ);
    uv = complex(cos(theta), sin(theta));

    [dummy, uv_sort_idx] = sort(abs(angle(conj(uv(1)) * uv)) + pi*rand(N,1)/5);

    n_left = (1:N-1)';
    n_right = N - n_left;

    uv_cum = cumsum(uv(uv_sort_idx));
    %uv_cum = cumsum(uv);

    uv_sum_left = uv_cum(n_left,:);
    uv_sum_right = uv_cum(end,:)-uv_sum_left;

    % Compute left/right circular dispersions the add them - note that abs is
    % *much* slower because of computing the sqrt so we use the square of the
    % dispersions as a proxy measure

    %1
    d_left1 = abs(uv_sum_left) ./ n_left;
    d_right1 = abs(uv_sum_right) ./ n_right;

    %2
    d_left2 = (real(uv_sum_left).^2 + imag(uv_sum_left).^2) ./ n_left;
    d_right2 = (real(uv_sum_right).^2 + imag(uv_sum_right).^2) ./ n_right;
    
    [dummy split1] = max(d_left1 + d_right1);
    [dummy split2] = max(d_left2 + d_right2);
    %
    figure; 
    subplot(1,2,1); hold on;
    plot(1:N-1, d_left1, 'r-x');
    plot(1:N-1, d_right1, 'g-x');
    plot(1:N-1, d_left1+d_right1, 'b-x');
    title(['Circle ' num2str(circ) ' - splitting % = ' num2str(split1/N) ' of ' num2str(N)]);

    subplot(1,2,2); hold on;
    plot(1:N-1, d_left2, 'r-x');
    plot(1:N-1, d_right2, 'g-x');
    plot(1:N-1, d_left2+d_right2, 'b-x');
    title(['Circle ' num2str(circ) ' - splitting % = ' num2str(split2/N) ' of ' num2str(N)]);
end
%%
% theta = sort(2*pi*rand(1e5,1));
% uv = complex(cos(theta), sin(theta));
load C:\isbe\asymmetry_project\data\misc\reg_data1e5.mat y
uv = complex(cosd(2*y), sind(2*y));
for ii = 0:12 
    N = length(uv);
    if N < 10
        break;
    end

    %[dummy, uv_sort_idx] = sort(abs(angle(conj(uv(1)) * uv)));% + pi*rand(N,1)/5
    %uv_sort_idx = randperm(N);
    [dummy, uv_sort_idx] = sort(abs(angle(conj(uv(1)) * uv)));%
    
    n_left = (1:N-1)';
    n_right = N - n_left;

    uv_cum = cumsum(uv(uv_sort_idx));

    uv_sum_left = uv_cum(n_left,:);
    uv_sum_right = uv_cum(end,:)-uv_sum_left;

    % Compute left/right circular dispersions the add them - note that abs is
    % *much* slower because of computing the sqrt so we use the square of the
    % dispersions as a proxy measure

    %1
    d_left1 = abs(uv_sum_left) ./ n_left;
    d_right1 = abs(uv_sum_right) ./ n_right;
    %[dummy split1] = max(d_left1 + d_right1);
    d_diff1 = -abs(d_left1 - d_right1);
    [dummy split1] = max(d_diff1);
    
    %2
    d_left2 = (real(uv_sum_left).^2 + imag(uv_sum_left).^2) ./ n_left;
    d_right2 = (real(uv_sum_right).^2 + imag(uv_sum_right).^2) ./ n_right;
    %[dummy split2] = max(d_left2 + d_right2);
    %
    d_diff2 = -abs(d_left2 - d_right2);
    [dummy split2] = max(d_diff2);
    
    figure; 
    subplot(2,2,1); hold on;
    plot(1:N-1, d_left1, 'r-x');
    plot(1:N-1, d_right1, 'g-x');
    title(['Magntiudes: Level ' num2str(ii)]);

    subplot(2,2,2); hold on;
    plot(1:N-1, d_diff1, 'b-x');
    title(['Splitting % = ' num2str(split1/N) ' of ' num2str(N)]);

    subplot(2,2,3); hold on;
    plot(1:N-1, d_left2, 'r-x');
    plot(1:N-1, d_right2, 'g-x');
    title(['Sum of squares: Level ' num2str(ii)]);
    
    subplot(2,2,4); hold on;
    plot(1:N-1, d_diff2, 'b-x');
    title(['Splitting % = ' num2str(split2/N) ' of ' num2str(N)]);
    
    uv = uv(uv_sort_idx(1:split2));
end
%%
for kk = 1:5
    %load tree
    tree = u_load([random_forest.tree_root random_forest.tree_dir random_forest.trees{kk}]);
    figure; hold on; axis equal;
    colors = 'br';
    leaves = find(~tree.var);
    for ii = 1:length(leaves)
        %Get leaf node and score at this leaf
        curr_node = leaves(ii);
        x_hat = real(tree.class(curr_node)) / abs(tree.class(curr_node));
        y_hat = imag(tree.class(curr_node)) / abs(tree.class(curr_node));

        %For each leaf, follow the branch up through the tree to the root
        branch_lr = [];
        
        while(curr_node ~= 1)
            branch_lr(end+1) = rem(curr_node,2)+1; %#ok
            curr_node = tree.parent(curr_node);    
        end
        branch_length = length(branch_lr);
        for jj = 1:branch_length
            len = [branch_length - jj branch_length - jj + 1];
            plot(x_hat*len, y_hat*len, colors(branch_lr(jj)));
        end

    end
end
%%
forest_dir = 'Z:\data\line_orientation_rfs\';
lines_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\real512';
results_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\';
figure;
a1 = subplot(1,2,1); hold on;
a2 = subplot(1,2,2); hold on;

rf_id = '243311';
rf = u_load([forest_dir rf_id '\random_forest.mat']);
rf.tree_root = forest_dir;
[ang_hist d_hist tree_means tree_medians] = forest_leaf_output(rf, 0);
errors = compute_image_orientation_errors(lines_dir, [results_dir rf_id]);
plot(a1, mean(abs(errors(:,1))), mean(tree_means), 'x');
plot(a2, mean(abs(errors(:,1))), mean(tree_medians), 'x');
figure;
subplot(1,2,1); bar(linspace(-pi/2, pi/2, 60), ang_hist); set(gca, 'xlim', [-pi/2 pi/2]);
subplot(1,2,2); bar(linspace(0, 1, 50), d_hist); set(gca, 'xlim', [0 1]);
title('Conjugate phase + magnitude, 3x3 windows');

rf_id = '243629';
rf = u_load([forest_dir rf_id '\random_forest.mat']);
rf.tree_root = forest_dir;
[ang_hist d_hist tree_means tree_medians] = forest_leaf_output(rf, 0);
errors = compute_image_orientation_errors(lines_dir, [results_dir rf_id]);
plot(a1, mean(abs(errors(:,1))), mean(tree_means), 'x');
plot(a2, mean(abs(errors(:,1))), mean(tree_medians), 'x');
figure;
subplot(1,2,1); bar(linspace(-pi/2, pi/2, 60), ang_hist); set(gca, 'xlim', [-pi/2 pi/2]);
subplot(1,2,2); bar(linspace(0, 1, 50), d_hist); set(gca, 'xlim', [0 1]);
title('Magnitude only, 3x3 windows');

rf_id = '243630';
rf = u_load([forest_dir rf_id '\random_forest.mat']);
rf.tree_root = forest_dir;
[ang_hist d_hist tree_means tree_medians] = forest_leaf_output(rf, 0);
errors = compute_image_orientation_errors(lines_dir, [results_dir rf_id]);
plot(a1, mean(abs(errors(:,1))), mean(tree_means), 'x');
plot(a2, mean(abs(errors(:,1))), mean(tree_medians), 'x');
figure;
subplot(1,2,1); bar(linspace(-pi/2, pi/2, 60), ang_hist); set(gca, 'xlim', [-pi/2 pi/2]);
subplot(1,2,2); bar(linspace(0, 1, 50), d_hist); set(gca, 'xlim', [0 1]);
title('Phase only, 3x3 windows');

rf_id = '286712';
rf = u_load([forest_dir rf_id '\random_forest.mat']);
rf.tree_root = forest_dir;
[ang_hist d_hist tree_means tree_medians] = forest_leaf_output(rf, 0);
errors = compute_image_orientation_errors(lines_dir, [results_dir rf_id]);
plot(a1, mean(abs(errors(:,1))), mean(tree_means), 'x');
plot(a2, mean(abs(errors(:,1))), mean(tree_medians), 'x');
figure;
subplot(1,2,1); bar(linspace(-pi/2, pi/2, 60), ang_hist); set(gca, 'xlim', [-pi/2 pi/2]);
subplot(1,2,2); bar(linspace(0, 1, 50), d_hist); set(gca, 'xlim', [0 1]);
title('Conjugate phase OR magnitude, 3x3 windows');

rf_id = '286713';
rf = u_load([forest_dir rf_id '\random_forest.mat']);
rf.tree_root = forest_dir;
[ang_hist d_hist tree_means tree_medians] = forest_leaf_output(rf, 0);
errors = compute_image_orientation_errors(lines_dir, [results_dir rf_id]);
plot(a1, mean(abs(errors(:,1))), mean(tree_means), 'x');
plot(a2, mean(abs(errors(:,1))), mean(tree_medians), 'x');
figure;
subplot(1,2,1); bar(linspace(-pi/2, pi/2, 60), ang_hist); set(gca, 'xlim', [-pi/2 pi/2]);
subplot(1,2,2); bar(linspace(0, 1, 50), d_hist); set(gca, 'xlim', [0 1]);
title('Magnitude only, 1x1 windows');

rf_id = '286714';
rf = u_load([forest_dir rf_id '\random_forest.mat']);
rf.tree_root = forest_dir;
[ang_hist d_hist tree_means tree_medians] = forest_leaf_output(rf, 0);
errors = compute_image_orientation_errors(lines_dir, [results_dir rf_id]);
plot(a1, mean(abs(errors(:,1))), mean(tree_means), 'x');
plot(a2, mean(abs(errors(:,1))), mean(tree_medians), 'x');
figure;
subplot(1,2,1); bar(linspace(-pi/2, pi/2, 60), ang_hist); set(gca, 'xlim', [-pi/2 pi/2]);
subplot(1,2,2); bar(linspace(0, 1, 50), d_hist); set(gca, 'xlim', [0 1]);
title('Phase only, 1x1 windows');

rf_id = '286715';
rf = u_load([forest_dir rf_id '\random_forest.mat']);
rf.tree_root = forest_dir;
[ang_hist d_hist tree_means tree_medians] = forest_leaf_output(rf, 0);
errors = compute_image_orientation_errors(lines_dir, [results_dir rf_id]);
plot(a1, mean(abs(errors(:,1))), mean(tree_means), 'x');
plot(a2, mean(abs(errors(:,1))), mean(tree_medians), 'x');
figure;
subplot(1,2,1); bar(linspace(-pi/2, pi/2, 60), ang_hist); set(gca, 'xlim', [-pi/2 pi/2]);
subplot(1,2,2); bar(linspace(0, 1, 50), d_hist); set(gca, 'xlim', [0 1]);
title('Conjugate phase + magnitude, 1x1 windows');

rf_id = '286716';
rf = u_load([forest_dir rf_id '\random_forest.mat']);
rf.tree_root = forest_dir;
[ang_hist d_hist tree_means tree_medians] = forest_leaf_output(rf, 0);
errors = compute_image_orientation_errors(lines_dir, [results_dir rf_id]);
plot(a1, mean(abs(errors(:,1))), mean(tree_means), 'x');
plot(a2, mean(abs(errors(:,1))), mean(tree_medians), 'x');
figure;
subplot(1,2,1); bar(linspace(-pi/2, pi/2, 60), ang_hist); set(gca, 'xlim', [-pi/2 pi/2]);
subplot(1,2,2); bar(linspace(0, 1, 50), d_hist); set(gca, 'xlim', [0 1]);
title('ICP phases, 3x3 windows');
%%
mkdir Z:\data\synthetic_lines\real512\masks\
for ii = 1:100
    s = load(['Z:\data\synthetic_lines\real512\labels\label' zerostr(ii,3) '.mat']);
    mask = s.label_centre & s.label == 1;
    save(['Z:\data\synthetic_lines\real512\masks\mask' zerostr(ii,3) '.mat'], 'mask');
end
%%

args.task_id = 1;
args.num_jobs = 100;
args.forest_dir = 'line_orientation_rfs';
args.mask_dir = 'synthetic_lines\real512\masks';
args.max_size = 512;
args.use_nag = 0;
tic;
classify_image_set('243311', 'synthetic_lines\real512', args);
toc;

    

