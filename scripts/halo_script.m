%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script for experimenting with mass halos
%
load C:\isbe\dev\files\u_files.mat
r_idx = randsample(101, 20);
%%

mean_halos = zeros(20, 6, 10);
var_halos = zeros(20, 6, 10);
seg_areas = zeros(20, 6, 10);
for ii = 1:10
    idx = r_idx(ii);
    load(['C:\isbe\dev\masses\', u_files1(idx).name]);
    [mean_halos(:,:,ii) var_halos(:,:,ii) seg_areas(:,:,ii)] =...
        mass_halo(mass, mass_model, idx,...
        20, .05, .15, 'if_plot', 1, 'startingPts', 0:20:100); %:20:100
    clear mass;
end

%%
X_halos = zeros(1000, 101);
seg_areas = zeros(500, 101);
for idx = 1:101
    load(['C:\isbe\dev\masses\', u_files1(idx).name]);
    [X_halos(1:500,idx) X_halos(501:end,idx) seg_areas(:,idx)] =...
        mass_halo(mass, mass_model, idx,...
        500, .05, .15, 'if_plot', 0, 'startingPts', 1);
    clear mass;
end
X_halos = X_halos';
%%
replace_mat = repmat(nanmean(X_halos), 101, 1);
X_halos(isnan(X_halos)) = replace_mat(isnan(X_halos));

halo_model.X_halos = X_halos;
halo_model.seg_areas = seg_areas;
[halo_model.mean_h, halo_model.P_h, halo_model.B_h, halo_model.L_h] =...
    pca(X_halos, 0.98);
save C:\isbe\dev\background\halo_model halo_model
%%

mass_outline = mass_model.P_shape(:,1:10)*mass_model.B_shape(1:10,idx)...
        + mass_model.mean_shape';
mass_outline = reshape(mass_outline, [], 2);

rotation = mass_model.rotations(:,:,idx);
translation = repmat(mass_model.translations(idx,:),...
    size(mass_outline, 1), 1);
origin = mass_model.origins(idx);
scale = mass_model.X_scale(idx);
%%
mass_outline1 = (circshift(mass_outline/scale, 1-origin)...
        - translation)*inv(rotation);
mass_outline2 = (circshift(mass_outline/(scale -.1), 1-origin)...
        - translation)*inv(rotation);
mass_outline3 = (circshift(mass_outline/(scale -.15), 1-origin)...
        - translation)*inv(rotation);
    
figure; hold on; axis equal;
plot(mass_outline1(:,1), mass_outline1(:,2), 'r');
plot(mass_outline2(:,1), mass_outline2(:,2), 'g');
plot(mass_outline3(:,1), mass_outline3(:,2), 'b');