function [mass_model]...
    = generate_test_AM(mass_files, file_out, size_shape_vec, size_tex_vec, path)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% author:   Michael Berks
% date:     26/05/2006  18:26
%
% function: implementation of Caulkin's method for modelling mass
%           appearance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Create matrix of shape vectors X_shape from data
%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('1: Method started');
mass_model.progress = '1: Method started';
N = length(mass_files);
X_shape = zeros(N, 2*size_shape_vec);
mass_centroid = zeros(N, 2);

for ii = 1:N
    temp = load([path, mass_files(ii).name]);
    mass = temp.mass; clear temp;
    
    shape_vec = mass.mass_outline;
    idx = round(linspace(1, length(shape_vec(:,1)), size_shape_vec+1));
    idx = idx(1:size_shape_vec); %ensures first point is not equal to the last point!
    X_shape(ii,:) = [shape_vec(idx,1)', shape_vec(idx,2)'];
    % Get number of spicules for each mass
    mass_centroid(ii, :) = mass.mass_centroid;
end

%
% Compute shape model from shape vectors
%%%%%%%%%%%%%%%%%%%%%%%%
[mean_shape, P_shape, B_shape, L_shape,...
    mean_scale, P_scale, B_scale, L_scale] = shape_model(X_shape, 20000, 0.98);
%
% Compute model of spicule numbers
%%%%
%[mean_n, P_n, B_n, L_n] = pca(n_spicules, 0.98);

mass_model.P_shape = P_shape;
mass_model.B_shape = B_shape;
mass_model.L_shape = L_shape;
mass_model.mean_scale = mean_scale;
mass_model.P_scale = P_scale;
mass_model.B_scale = B_scale;
mass_model.L_scale = L_scale;

display('2: Shape model complete');
mass_model.progress = '2: Shape model complete';
save(file_out, 'mass_model');

%
% Compute texture model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create mean shape pixel list NEED TO CHOOSE SIZE OF TEX VEC
%%%%%%%%%%

scale_factor = sqrt(size_tex_vec / ...
    polyarea(mean_shape(1:size_shape_vec),...
    mean_shape(size_shape_vec + 1: 2*size_shape_vec)));

mean_shape = mean_shape*scale_factor;
mass_model.mean_shape = mean_shape;
mass_model.scale_factor = scale_factor;

mean_off_r = 101 - floor(min(mean_shape(size_shape_vec+1:2*size_shape_vec)));
mean_off_c = 101 - floor(min(mean_shape(1:size_shape_vec)));

mass_model.mean_off_r = mean_off_r;
mass_model.mean_off_c = mean_off_c;

mean_r = ceil(max(mean_shape(size_shape_vec+1:2*size_shape_vec)))...
    + mean_off_r + 100;
mean_c = ceil(max(mean_shape(1:size_shape_vec)))...
    + mean_off_c + 100;

mean_bw = poly2mask(mean_shape(1:size_shape_vec) + mean_off_c,...
    mean_shape(size_shape_vec+1:2*size_shape_vec) + mean_off_r, mean_r, mean_c);

mean_props = regionprops(bwlabel(mean_bw, 4), 'PixelList');

mean_shape_pl = mean_props.PixelList; clear mean_props;

mask1 = mean_bw;

for jj = 1:40
    mask1 = imdilate(mask1, strel('disk', 1));
end

mean_ring = mask1 - mean_bw; clear mean_bw mask1;

mean_props = regionprops(bwlabel(mean_ring, 4), 'PixelList'); clear mean_ring;
mean_ring_pl = mean_props.PixelList; clear mean_props;

mass_model.mean_shape_pl = mean_shape_pl;
mass_model.mean_ring_pl = mean_ring_pl;

display('4: Mean pixel list created');
mass_model.progress = '4: Mean pixel list created';
save(file_out, 'mass_model');


% Warp each shape to mean shape and extract texture vector
%%%%%%%%%%

%Define source points for TPS - as row vectors
s_x = mean_shape(1:size_shape_vec) + mean_off_c;
s_y = mean_shape(size_shape_vec+1:2*size_shape_vec) + mean_off_r;

%Define points to be interpolated by TPS - as row vectors
i_x = mean_shape_pl(:,1)';
i_y = mean_shape_pl(:,2)';

r_x = mean_ring_pl(:,1)';
r_y = mean_ring_pl(:,2)';

tps_L_inv = tps_weights(s_x, s_y);

%f1 = figure; title('Mass ROIs');
%f2 = figure; title('Subtracted and warped mass ROIs');
X_tex = zeros(N, length(i_x));
X_ring = zeros(N, length(r_x));

for ii=1:N
    display(['warping mean to shape ', int2str(ii)]);
    %Define displacement to target points
    z_x = X_shape(ii,1:size_shape_vec) + mass_centroid(ii,1);
    z_y = X_shape(ii,size_shape_vec+1:2*size_shape_vec) + mass_centroid(ii,2);
    
    %Compute displacement of interpolated points
    f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
    f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
    
    temp = load([path, mass_files(ii).name]);
    mass = temp.mass; clear temp;
    
    figure; hold on;
    plot(f_x(1:10:end), f_y(1:10:end), 'r.');
    plot(mass.mass_outline(:,1) + mass.mass_centroid(1),...
         mass.mass_outline(:,2) + mass.mass_centroid(2));
    plot(z_x, z_y, 'g');
    
    %Get texture vector
    X_tex(ii,:) = interp2(double(mass.subtract_ROI), f_x, f_y);
    
    f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, r_x, r_y);
    f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, r_x, r_y);
    
    f_x(f_x < 1) = 1;
    f_x(f_x > size(mass.mass_ROI, 2)) = size(mass.mass_ROI, 2);
    f_y(f_y < 1) = 1;
    f_y(f_y > size(mass.mass_ROI, 1)) = size(mass.mass_ROI, 1);
    
    X_ring(ii,:) = interp2(double(mass.mass_ROI), f_x, f_y);
    
%     if ii <= 20
%         temp = zeros(mean_r, mean_c);
%         temp(sub2ind([mean_r, mean_c], i_y, i_x)) = uint8(X_tex(ii,:));
%         figure(f1); subplot(4, 5, ii);
%         image(mass.mass_ROI); axis image; colormap(gray(256));
%         hold on; plot(f_x, f_y, 'r.');
%         figure(f2); subplot(4,5,ii);
%         image(temp); colormap(gray(256)); axis image; clear temp;
%         hold on; plot(r_x, r_y, 'r.');
%     end
end
clear s_x s_y z_x z_y i_x i_y f_x f_y tps_L_inv mass;

mass_model.warp_tex = X_tex;
display('5: Shapes warped to mean shape');
mass_model.progress = '5: Shapes warped to mean shape';
save(file_out, 'mass_model');
%
%Compute model of texture vectors
%%%%%%%%%%%%%%%%%%%%%
[mean_tex, P_tex, B_tex, L_tex] = pca(X_tex, 0.98); title('Modes of texture');

mass_model.mean_tex = mean_tex;
mass_model.P_tex = P_tex;
mass_model.B_tex = B_tex;
mass_model.L_tex = L_tex;

[mean_ring, P_ring, B_ring, L_ring] = pca(X_ring, 0.98); title('Modes of mass ring');

mass_model.mean_ring = mean_ring;
mass_model.P_ring = P_ring;
mass_model.B_ring = B_ring;
mass_model.L_ring = L_ring;

display('6: Texture model complete');
mass_model.progress = '6: Texture model complete';
save(file_out, 'mass_model');

%
%Calculate weights for combined model
%%%%%%%%%%%%%%%%%%%%%%
k_shape = length(L_shape);
k_tex   = length(L_tex);
k_ring  = length(L_ring);

W_shape = k_shape / sum(sqrt(L_shape));
W_tex   = k_tex / sum(sqrt(L_tex));
W_ring   = k_ring / sum(sqrt(L_ring));
W_scale = 1 / sqrt(L_scale); %length L_scale = 1
%W_n = 1 / sqrt(L_n); %length L_n = 1

%W_shape = k_shape / sum(L_shape);
%W_tex   = k_tex / sum(L_tex);
%W_scale = 1 / L_scale; %length(L_scale = 1)
%W_n = 1 / L_n; %length L_n = 1

mass_model.W_shape = W_shape;
mass_model.W_tex = W_tex;
mass_model.W_ring = W_ring;
mass_model.W_scale = W_scale;
%mass_model.W_n = W_n;

combined_data = [W_shape*B_shape; W_tex*B_tex; W_scale*B_scale;]';% W_n*B_n]';
%combined_data = [W_shape*B_shape; W_scale*B_scale;]';%
mass_model.combined_data = combined_data;

[mean_c, P_c, B_c, L_c] = pca(combined_data, 0.98); title('Combined appearance modes');

mass_model.mean_c = mean_c;
mass_model.P_c = P_c;
mass_model.B_c = B_c;
mass_model.L_c = L_c;

%C_data = [B_tex; B_shape; B_ring; B_scale];
%C_data = [B_tex; B_c; B_ring];
%C_data = [B_tex; B_shape; B_ring];
%C_com = cov([B_tex; B_shape; B_ring]');

C_com = cov([B_c; B_ring]');
%C_com = cov([B_tex; B_c; B_ring]');
mass_model.C_com = C_com; clear C_com;

display('7: Combined model complete, function successful!');
mass_model.progress = '7: Combined model complete, function successful!';
save(file_out, 'mass_model');

