%--------------------------------------------------------------------------
% Code showing how the dual-tree coefficients behave at crossing points
% - should probably extend this to look at other i2D features such as blobs
%--------------------------------------------------------------------------
% Key point: coefficients than more than just the maximal sub-band may be
% useful...

cross = ones(128);
cross(33:96, 64:65) = 2;
cross(64:65, 33:96) = 2;

dt_cross = dtwavexfm2(cross, 5, 'near_sym_b','qshift_b');
[ilp_cross icp_cross] = mb_dual_tree_transform(dt_cross);

colors = 'rgbmcy';
f1 = figure;
f2 = figure;

for level = 2:4
    dim = 2^level;
    st = (dim+1)/2;

    [max_icp_cross] = sort(ilp_cross{level}, 3, 'descend');
    [max_ilp_cross] = sort(icp_cross{level}, 3, 'descend');
    
    figure(f1);
    subplot(1,3,level - 1);
    imagesc(cross); axis image; colormap(gray(256)); hold on;
    for ori = 1:6
        quiver(st:dim:129-st, st:dim:129-st, real(ilp_cross{level}(:,:,ori)), -imag(ilp_cross{level}(:,:,ori)), colors(ori), 'LineWidth', 2.0);
    end
    figure(f2);
    subplot(1,3,level - 1);
    imagesc(cross); axis image; colormap(gray(256)); hold on;
    for ori = 1:6
        quiver(st:dim:129-st, st:dim:129-st, real(icp_cross{level}(:,:,ori)), -imag(icp_cross{level}(:,:,ori)), colors(ori), 'LineWidth', 2.0);
    end
end
%%
for level = 1:4
figure; hold on; axis equal;
    for ori = 1:6
        compass(icp_cross{level}(64 / 2^level, 64 / 2^level, ori), colors(ori));
    end
    title(['Orientations in level ', num2str(level)]);
    legend({'15 degrees', '45 degrees', '75 degrees', '105 degrees', '135 degrees', '165 degree'});
end
%%
for level = 1:4
figure; hold on; axis equal;
    for ori = 1:6
        compass(ilp_cross{level}(64 / 2^level, 64 / 2^level, ori), colors(ori));
    end
    title(['Feature type in level ', num2str(level)]);
    legend({'15 degrees', '45 degrees', '75 degrees', '105 degrees', '135 degrees', '165 degree'});
end
%%
cross = 2*ones(128);
cross(33:96, 64:65) = 1;
cross(64:65, 33:96) = 3;

dt_cross = dtwavexfm2(cross, 5, 'near_sym_b','qshift_b');
[ilp_cross icp_cross] = mb_dual_tree_transform(dt_cross);

for level = 1:4
figure; hold on; axis equal;
    for ori = 1:6
        compass(icp_cross{level}(64 / 2^level, 64 / 2^level, ori), colors(ori));
    end
    title(['Orientations in level ', num2str(level)]);
    legend({'15 degrees', '45 degrees', '75 degrees', '105 degrees', '135 degrees', '165 degree'});
end
%%
for level = 1:4
figure; hold on; axis equal;
    for ori = 1:6
        compass(ilp_cross{level}(64 / 2^level, 64 / 2^level, ori), colors(ori));
    end
    title(['Feature type in level ', num2str(level)]);
    legend({'15 degrees', '45 degrees', '75 degrees', '105 degrees', '135 degrees', '165 degree'});
end
%%
cross = zeros(128);
cross(33:96, 64:65) = 1;
cross(64:65, 33:96) = 1;
dt_cross = dtwavexfm2(cross, 5, 'near_sym_b','qshift_b');
cross2 = dtwaveifm2(dt_cross);
figure; imagesc(cross2); axis image;