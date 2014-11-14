%%
k_shape = length(mass_model.L_shape);
k_tex   = length(mass_model.L_tex);
k_ring   = length(mass_model.L_ring);

W_shape = k_shape / sum(sqrt(mass_model.L_shape));
W_tex   = k_tex / sum(sqrt(mass_model.L_tex));
W_scale = 1 / sqrt(mass_model.L_scale); %length L_scale = 1
W_ring   = k_ring / sum(sqrt(mass_model.L_ring));

mass_model.W_shape = W_shape;
mass_model.W_tex = W_tex;
mass_model.W_scale = W_scale;
mass_model.W_ring = W_ring;

combined_data = [W_shape*mass_model.B_shape; W_tex*mass_model.B_tex; W_scale*mass_model.B_scale]';
mass_model.combined_data = combined_data;

[mean_c, P_c, B_c, L_c] = pca(combined_data, 0.98);%, 0);

mass_model.mean_c = mean_c;
mass_model.P_c = P_c;
mass_model.B_c = B_c;
mass_model.L_c = L_c;

display('7: Combined model complete, function successful!');
mass_model.progress = '7: Combined model complete, function successful!';
save('C:\isbe\dev\halo_model.mat', 'mass_model');
%%
for ii = 1:1
    load(['C:\isbe\dev\annotations\', files(ii).name]);
    
    mask1 = roipoly(mass.mass_ROI,...
        mass.mass_outline(:,1), mass.mass_outline(:,2));
    mask2 = mask1;
    
    for jj = 1:100
        mask2 = imdilate(mask2, strel('disk', 1));
    end
    
    [sr sc] = ind2sub(size(mask1), find(mask1, 1));
    inner_ring = bwtraceboundary(mask1, [sr sc], 'N', 4, Inf, 'counterclockwise');
    inner_ring = fliplr(inner_ring);
    
    [sr sc] = ind2sub(size(mask2), find(mask2, 1));
    outer_ring1 = bwtraceboundary(mask2, [sr sc], 'N', 4, Inf, 'counterclockwise');
    outer_ring1 = fliplr(outer_ring1);
    
    clear mask1 mask2;
    
    outer_ring2 = repmat(mean(inner_ring), size(inner_ring,1), 1) +...
        1.25 * (inner_ring - repmat(mean(inner_ring), size(inner_ring,1), 1));
    
    figure('WindowStyle', 'docked');
    image(mass.mass_ROI); axis('image'); colormap(gray(256));
    hold on;
    plot(inner_ring(:,1), inner_ring(:,2), 'r');
    plot(outer_ring1(:,1), outer_ring1(:,2), 'g');
    plot(outer_ring2(:,1), outer_ring2(:,2), 'y');
end
    