%%
% Comparison of shape alignment algorithms

% Testing the stability of the algorithm for different seed points:

% No origin shift
load C:\isbe\dev\shapes\shapes_unaligned.mat
load C:\isbe\dev\files\u_files.mat
r_idx = randsample(1:101, 20); %indices to select random seed shapes
N_pts = size(shapes_u, 2)/2;
colors = lines(20);
%
shape_means_p = zeros(20, 2*N_pts);
shape_vars_p = zeros(20, 2*N_pts);
shape_L_p = zeros(20, 100);

figure; axis equal; hold on;
for ii = 1:20
    shapes_a = align_shapes(shapes_u, 'shiftOrigin', 0, 'seed', r_idx(ii));
    shape_means_p(ii,:) = mean(shapes_a);
    shape_vars_p(ii,:) = var(shapes_a);
    [model.m model.P model.B model.L] = pca(shapes_a, 1);
    shape_L_p(ii,:) = model.L';
    clear shapes_a model
    plot(shape_means_p(ii,1:N_pts), shape_means_p(ii,N_pts+1:end), 'Color', colors(ii,:));
end
%
shape_means_o = zeros(20, 2*N_pts);
shape_vars_o = zeros(20, 2*N_pts);
shape_L_o = zeros(20, 100);

figure; axis equal; hold on;
for ii = 1:20
    shapes_a = align_shapes(shapes_u, 'shiftOrigin', 1, 'seed', r_idx(ii));
    shape_means_o(ii,:) = mean(shapes_a);
    shape_vars_o(ii,:) = var(shapes_a);
    [model.m model.P model.B model.L] = pca(shapes_a, 1);
    shape_L_o(ii,:) = model.L';
    clear shapes_a model
    plot(shape_means_o(ii,1:N_pts), shape_means_o(ii,N_pts+1:end), 'Color', colors(ii,:));
end
%
shape_means_i = zeros(20, 2*N_pts);
shape_vars_i = zeros(20, 2*N_pts);
shape_L_i = zeros(20, 100);

figure; axis equal; hold on;
for ii = 1:20
    shapes_a = align_shapes(shapes_u, 'shiftOrigin', 1, 'seed', r_idx(ii));
    shape_means_i(ii,:) = mean(shapes_a);
    shape_vars_i(ii,:) = var(shapes_a);
    [model.m model.P model.B model.L] = pca(shapes_a, 1);
    shape_L_i(ii,:) = model.L';
    clear shapes_a model
    plot(shape_means_i(ii,1:N_pts), shape_means_i(ii,N_pts+1:end), 'Color', colors(ii,:));
end
%%
shape_means_m = zeros(20, 2*N_pts);
shape_vars_m = zeros(20, 2*N_pts);
shape_L_m = zeros(20, 100);

figure; axis equal; hold on;
for ii = 1:20
    shapes_a = amb_automodelbuild (reshape_MDL(shapes_pro), 'circle', 'nIterations', 1e3,...
        'saveFrequency', 50, 'Quiet', 1,...
        'optimisePose', 1, 'optimiseOrigin', 1, 'nExamplesToOptimise', 101,...
        'initialAlign', 0, 'initialOriginOptimisation', 0,...
        'nSamplePoints', 500, 'totalMaxFevals', 1e5, 'saveDirName', 'C:\isbe\dev\shapes',...
        'optimiseParameterisation', 0,...
        'poseOptimisationsPerIteration', 1,...
        'originOptimisationsPerIteration', 1);
    shapes_a = reshape_MDL(shapes_a);
    shape_means_m(ii,:) = mean(shapes_a);
    shape_vars_m(ii,:) = var(shapes_a);
    [model.m model.P model.B model.L] = pca(shapes_a, 1);
    shape_L_m(ii,:) = model.L';
    clear shapes_a model
    plot(shape_means_m(ii,1:N_pts), shape_means_m(ii,N_pts+1:end), 'Color', colors(ii,:));
end
