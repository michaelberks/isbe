%%
figure;
subplot(2, 2, 1); hist(model.B_ori(1,:), 30);
subplot(2, 2, 2); hist(model.B_ori(2,:), 30);
subplot(2, 2, 3); hist(model.B_ori(3,:), 30);
subplot(2, 2, 4); hist(model.B_ori(4,:), 30);

figure;
subplot(2, 2, 1); hist(model.B_pro(1,:), 30);
subplot(2, 2, 2); hist(model.B_pro(2,:), 30);
subplot(2, 2, 3); hist(model.B_pro(3,:), 30);
subplot(2, 2, 4); hist(model.B_pro(4,:), 30);
%%
[r_pts, a, a, a, a, a, shapes] = ...
    amb_automodelbuild (reshape_MDL(shapes_aligned),...
    'circle', 'nIterations', 1, 'saveFrequency', 50,...
    'Quiet', 0, 'optimisePose', 1, 'optimiseOrigin', 0,...
    'nExamplesToOptimise', 1, 'initialAlign', 0,...
    'initialOriginOptimisation', 0, 'nSamplePoints', 200,...
    'totalMaxFevals', 1e5);
%%
figure; hold on;
for ii = 1:101
    plot(r_pts(1,:,ii), r_pts(2,:,ii));
    plot(r_pts(1,1,ii), r_pts(2,1,ii), 'rx');
end
figure; hold on;
for ii = 1:101
    plot(shapes(1,:,ii), shapes(2,:,ii));
    plot(shapes(1,1,ii), shapes(2,1,ii), 'rx');
end
%%
figure; hold on;
for ii = 1:101
    plot(mdl.reparameterisedPoints(1,:,ii), mdl.reparameterisedPoints(2,:,ii));
    plot(mdl.reparameterisedPoints(1,1,ii), mdl.reparameterisedPoints(2,1,ii), 'rx');
end
figure; hold on;
for ii = 1:101
    plot(mdl.shapes(1,:,ii), mdl.shapes(2,:,ii));
    plot(mdl.shapes(1,1,ii), mdl.shapes(2,1,ii), 'rx');
end
