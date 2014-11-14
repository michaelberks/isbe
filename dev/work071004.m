

[r_pts] = amb_automodelbuild (shapes, 'circle', 'nIterations', 200, 'saveFrequency', 10,...
    'Quiet', 0, 'optimisePose', 0, 'optimiseOrigin', 0, 'nExamplesToOptimise', 25,...
    'initialAlign', 1, 'initialOriginOptimisation', 1);
shapes_mdl = zeros(50, 1000);
for ii = 1:50; shapes_mdl(ii,:) = [r_pts(1,:,ii) r_pts(2,:,ii)]; end

[mean_xp, P_xp, B_xp, L_xp] = pca(shapes_pro, 0.98);
[mean_xm, P_xm, B_xm, L_xm] = pca(shapes_mdl, length(L_xp));
%%
[P_xps, L_xps, B_xps, mean_xps] = st_pca(shapes_pro, 0.98);
[P_xms, L_xms, B_xms, mean_xms] = st_pca(shapes_mdl, 0.98, 'nmodes', [length(L_xp) length(L_xp)]);

rg_shapes_pro = repmat(mean_xps, 50, 1) + (P_xps*B_xps')';
rg_shapes_mdl = repmat(mean_xms, 50, 1) + (P_xms*B_xms')';

for ii = 1:50
%     figure; 
%     subplot(1, 2, 1)
%     hold on;
%     plot(shapes_pro(ii,1:500), shapes_pro(ii,501:end), 'r');
%     plot(rg_shapes_pro(ii,1:500), rg_shapes_pro(ii,501:end), 'b');
    dp(ii) = sum(sqrt( (shapes_pro(ii,1:500) - rg_shapes_pro(ii,1:500)).^2 + ...
                   (shapes_pro(ii,501:end) - rg_shapes_pro(ii,501:end)).^2 ));
%     subplot(1, 2, 2)
%     hold on;
%     plot(shapes_mdl(ii,1:500), shapes_mdl(ii,501:end), 'g');
%     plot(rg_shapes_mdl(ii,1:500), rg_shapes_mdl(ii,501:end), 'm');
    dm(ii) = sum(sqrt( (shapes_mdl(ii,1:500) - rg_shapes_mdl(ii,1:500)).^2 + ...
                   (shapes_mdl(ii,501:end) - rg_shapes_mdl(ii,501:end)).^2 ));
               
    display(['procrustes error = ', num2str(dp(ii)), '; MDL error = ', num2str(dm(ii))]);
end
%%
for ii = 1:20
    figure; 
    subplot(1,2,1); hold on;
    plot(shapes_t.unaligned(ii,1:500), shapes_t.unaligned(ii, 501:end), 'r');
    plot(shapes_t.pro(ii,1:500), shapes_t.pro(ii, 501:end), 'g');
    plot(shapes_rg.pro(ii,1:500), shapes_rg.pro(ii,501:end), 'b');
    
    subplot(1,2,2); hold on;
    plot(shapes_t.unaligned(ii,1:500), shapes_t.unaligned(ii, 501:end), 'r');
    plot(shapes_t.mdl(ii,1:500), shapes_t.mdl(ii, 501:end), 'g');
    plot(shapes_rg.mdl(ii,1:500), shapes_rg.mdl(ii,501:end), 'b');
    
end