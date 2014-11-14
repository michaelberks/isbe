%
% Script to compare the performance of varying MDL optimisations
%

% Load shapes and align using procrustes method (with shifted origin)
load C:\isbe\dev\files\u_files.mat
N = length(u_files1);

shapes_unaligned = get_shapes_from_masses(u_files1, 1000);
shapes_aligned = align_shapes(shapes_unaligned, 'shiftOrigin', 1);
%%
% How many sample points to use?
%
for n_pts = 100:100:1000

dirName = ['C:\isbe\dev\mdl\nSamplePoints_', num2str(n_pts, 4)];
    
amb_automodelbuild (reshape_MDL(shapes_aligned), 'circle', 'nIterations', Inf,...
    'saveFrequency', 50, 'Quiet', 1,...
    'optimisePose', 1, 'optimiseOrigin', 0, 'nExamplesToOptimise', N,...
    'initialAlign', 0, 'initialOriginOptimisation', 0,...
    'nSamplePoints', n_pts, 'totalMaxFevals', 1e5, 'saveDirName', dirName);
end

for n_ex = 1:10:101

dirName = ['C:\isbe\dev\mdl\nExamples_', num2str(n_ex, 3)];
    
amb_automodelbuild (reshape_MDL(shapes_aligned), 'circle', 'nIterations', Inf,...
    'saveFrequency', 50, 'Quiet', 1,...
    'optimisePose', 1, 'optimiseOrigin', 0, 'nExamplesToOptimise', n_ex,...
    'initialAlign', 0, 'initialOriginOptimisation', 0,...
    'nSamplePoints', 200, 'totalMaxFevals', 1e5, 'saveDirName', dirName);
end
%%
% Use best parameters with origin optimisation on
dirName = 'C:\isbe\dev\mdl\optimisedOrigin';
    
amb_automodelbuild (reshape_MDL(shapes_aligned), 'circle', 'nIterations', Inf,...
    'saveFrequency', 50, 'Quiet', 1,...
    'optimisePose', 1, 'optimiseOrigin', 1, 'nExamplesToOptimise', 101,...
    'initialAlign', 0, 'initialOriginOptimisation', 0,...
    'nSamplePoints', 1000, 'totalMaxFevals', 1e5, 'saveDirName', dirName);

%%
% Load reparameterised points and compute loo errors
er_mean_n_pts = zeros(9, 10);
er_stdv_n_pts = zeros(9, 10);
var_mean_n_pts = zeros(9, 10);

for n_pts = 100:100:1000

    dirName = ['C:\isbe\dev\mdl\nSamplePoints_', num2str(n_pts, 4)];
    mdl = load([dirName, '\ambsave_update.mat']);
    [er_mdl var_mdl] = compute_loo_shape_error(reshape_MDL(mdl.shapes),...
        [5:5:40, 1]);
    save([dirName, '\loo_error'], 'er_mdl', 'var_mdl');
    er_mean_n_pts(:,n_pts/100) = mean(er_mdl)';
    er_stdv_n_pts(:,n_pts/100) = std(er_mdl)';
    var_mean_n_pts(:,n_pts/100) = mean(var_mdl)';

    clear er_mdl var_mdl mdl
end

figure; hold on;
plot(er_mean_n_pts, '-x');
title('Mean leave-one-out errors for varying number of sample of points');
figure; hold on;
plot(er_stdv_n_pts, '-x');
title('SD of leave-one-out errors for varying number of sample of points');
figure; hold on;
plot(var_mean_n_pts, '-x');
title('Mean variance in leave-one-out model for varying number of sample of points');
%%
er_mean_n_exs = zeros(9, 11);
er_stdv_n_exs = zeros(9, 11);
var_mean_n_exs = zeros(9, 11);

n_examples = 1:10:101;
for kk = 1:11

    n_exs = n_examples(kk);
    dirName = ['C:\isbe\dev\mdl\nExamples_', num2str(n_exs, 4)];
    mdl = load([dirName, '\ambsave_update.mat']);
    [er_mdl var_mdl] = compute_loo_shape_error(reshape_MDL(mdl.shapes),...
        [5:5:40, 1]);
    save([dirName, '\loo_error'], 'er_mdl', 'var_mdl');
    er_mean_n_exs(:,kk) = mean(er_mdl)';
    er_stdv_n_exs(:,kk) = std(er_mdl)';
    var_mean_n_exs(:,kk) = mean(var_mdl)';

    clear er_mdl var_mdl mdl
end

figure; hold on;
plot(er_mean_n_exs, '-x');
title('Mean leave-one-out errors for varying number of examples');
xlabel('Number of modes retained');

figure; hold on;
plot(er_stdv_n_exs, '-x');
title('SD of leave-one-out errors for varying number of examples');
xlabel('Number of modes retained');

figure; hold on;
plot(var_mean_n_exs, '-x');
title('Mean variance in leave-one-out model for varying number of examples');
xlabel('Number of modes retained');
%%
% Since those graphs don't show much. Just show original procrustes error,
% procrustes error with shifted origin and one MDL (the origin optimised)
% error
[er_pro var_pro] = compute_loo_shape_error(shapes_pro, [1.1,2:100]);
[er_ori var_ori] = compute_loo_shape_error(shapes_aligned, [1.1,2:100]);
%%
mdl = load('C:\isbe\dev\mdl\optimisedOrigin\ambsave_update.mat');
[er_mdl var_mdl] = compute_loo_shape_error(reshape_MDL(mdl.reparameterisedPoints), [1.1,2:100]);

errors = [mean(er_pro); mean(er_ori); mean(er_mdl)]';
vars = [mean(var_pro); mean(var_ori); mean(var_mdl)]';

save C:\isbe\dev\mdl\loo_errors2.mat er* var*
save E:\mdl\loo_errors2 er* var*



