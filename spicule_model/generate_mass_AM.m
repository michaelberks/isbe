function [mass_model]...
    = generate_mass_AM(mass_files, file_out, size_shape_vec, path, size_tex_vec)
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
mass_areas = zeros(N, 1);
n_spicules = zeros(N, 1);
mass_centroid = zeros(N, 2);

for ii = 1:N
    temp = load([path, mass_files(ii).name]);
    mass = temp.mass; clear temp;
    
    shape_vec = mass.mass_outline;
    idx = round(linspace(1, length(shape_vec(:,1)), size_shape_vec+1));
    idx = idx(1:size_shape_vec); %ensures first point is not equal to the last point!
    X_shape(ii,:) = [shape_vec(idx,1)', shape_vec(idx,2)'];
    mass_areas(ii) = mass.mass_area;
    % Get number of spicules for each mass
    n_spicules(ii) = 0; %mass.n_spicules;
    mass_centroid(ii, :) = mass.mass_centroid;
    
    clear mass;
end

%
% Compute shape model from shape vectors
%%%%%%%%%%%%%%%%%%%%%%%%
[mean_shape, P_shape, B_shape, L_shape,...
    mean_scale, P_scale, B_scale, L_scale, mean_target]...
    = shape_model(X_shape, mean(mass_areas),  0.98);
%
% Compute model of spicule numbers
%%%%
%[mean_n, P_n, B_n, L_n] = pca(n_spicules, 0.98);

mass_model.X_shape = X_shape;
mass_model.P_shape = P_shape;
mass_model.B_shape = B_shape;
mass_model.L_shape = L_shape;
mass_model.mean_scale = mean_scale;
mass_model.P_scale = P_scale;
mass_model.B_scale = B_scale;
mass_model.L_scale = L_scale;
mass_model.mean_target = mean_target;

% mass_model.mean_n = mean_n;
% mass_model.P_n = P_n;
% mass_model.B_n = B_n;
% mass_model.L_n = L_n;

display('2: Shape model complete');
mass_model.progress = '2: Shape model complete';
save(file_out, 'mass_model');

%
% Compute texture vectors by inerpolating grey-levels in ROI of mass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('3: Background subtraction complete');
mass_model.progress = '3: Background subtraction complete';
save(file_out, 'mass_model');

%
% Compute texture model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create mean shape pixel list
%%%%%%%%%%

if nargin < 5
    scale_factor = 1;
else
    scale_factor = sqrt(size_tex_vec / ...
        polyarea(mean_shape(1:size_shape_vec),...
        mean_shape(size_shape_vec + 1: 2*size_shape_vec)));
end

mean_shape = mean_shape*scale_factor;
mass_model.mean_shape = mean_shape;
mass_model.scale_factor = scale_factor;

mean_off_r = 101 - floor(min(mean_shape(size_shape_vec+1:2*size_shape_vec)));
mean_off_c = 101 - floor(min(mean_shape(1:size_shape_vec)));

mean_row = ceil(max(mean_shape(size_shape_vec+1:end)))...
    + mean_off_r + 100;
mean_col = ceil(max(mean_shape(1:size_shape_vec)))...
    + mean_off_c + 100;

mass_model.mean_off_r = mean_off_r;
mass_model.mean_off_c = mean_off_c;
mass_model.mean_row = mean_row;
mass_model.mean_col = mean_col;

mean_bw = poly2mask(mean_shape(1:size_shape_vec) + mean_off_c,...
    mean_shape(size_shape_vec+1:end) + mean_off_r, mean_row, mean_col);

for ii = 1:49; mean_bw = imdilate(mean_bw, strel('disk', 1)); end

[start_r start_c] = find(mean_bw, 1);
dilate_outline = bwtraceboundary(mean_bw, [start_r, start_c], 'E');
dilate_outline = [dilate_outline(:,2) dilate_outline(:,1)];

mean_props = regionprops(bwlabel(mean_bw, 4), 'PixelList'); clear mean_bw;

mean_shape_pl = mean_props.PixelList; clear mean_props;

mass_model.mean_shape_pl = mean_shape_pl;

display('4: Mean pixel list created');
mass_model.progress = '4: Mean pixel list created';
save(file_out, 'mass_model');

px = mean_shape(1) + mean_off_c;
py = mean_shape(size_shape_vec+1) + mean_off_r;
dd = sum((dilate_outline - repmat([px py], size(dilate_outline,1), 1)).^2, 2);
[mm idx] = min(dd);

dilate_outline = [dilate_outline(idx:end, :); dilate_outline(1:idx-1,:)];

idx = round(linspace(1, size(dilate_outline,1), size_shape_vec+1));
idx(end) = [];

mean_dilate = [dilate_outline(idx,1)', dilate_outline(idx,2)'];
mass_model.mean_dilate = mean_dilate;

% Warp each shape to mean shape and extract texture vector
%%%%%%%%%%

%Define source points for TPS - as row vectors
s_x = mean_dilate(1:size_shape_vec);
s_y = mean_dilate(size_shape_vec+1:end);

%Define points to be interpolated by TPS - as row vectors
i_x = mean_shape_pl(:,1)';
i_y = mean_shape_pl(:,2)';

tps_L_inv = tps_weights(s_x, s_y);

X_tex = zeros(N, length(i_x));
for ii=1:N
    
    display(['warping mean to shape ', int2str(ii)]);
    
    temp = load([path, mass_files(ii).name]);
    mass = temp.mass; clear temp;
    
    idx = round(linspace(1, size(mass.dilate_outline,1), size_shape_vec+1));
    idx(end) = [];
    x_dilate = mass.dilate_outline(idx,:);
    clear idx
    
    %Define displacement to target points
    z_x = x_dilate(:,1)';
    z_y = x_dilate(:,2)';
    clear x_dilate;    
    
    %Compute displacement of interpolated points
    f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
    f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
    
    %Get texture vector
    X_tex(ii,:) = interp2(mass.subtract_ROI, f_x, f_y);
    
%     if ~rem(ii, 20)
%         temp = zeros(mean_row, mean_col);
%         temp(sub2ind([mean_row, mean_col], i_y, i_x)) = uint8(X_tex(ii,:));
%         figure('WindowStyle', 'docked');
%         image(temp); axis image; colormap(gray(256));
%         clear temp;
%     end
        figure; hold on;
        plot(f_x(1:10:end), f_y(1:10:end), 'r.');
        plot(mass.mass_outline(:,1) + mass.mass_centroid(1),...
             mass.mass_outline(:,2) + mass.mass_centroid(2));
        plot(z_x, z_y, 'g');
    clear mass f_x f_y;
    
end
clear s_x s_y z_x z_y i_x i_y tps_L_inv;

X_tex(isnan(X_tex)) = 0;

mass_model.X_tex = X_tex;
display('5: Shapes warped to mean shape');
mass_model.progress = '5: Shapes warped to mean shape';
save(file_out, 'mass_model');
%
%Compute model of texture vectors
%%%%%%%%%%%%%%%%%%%%%
[mean_tex, P_tex, B_tex, L_tex] = pca(X_tex, 0.98);%, 0);

mass_model.mean_tex = mean_tex;
mass_model.P_tex = P_tex;
mass_model.B_tex = B_tex;
mass_model.L_tex = L_tex;

display('6: Texture model complete');
mass_model.progress = '6: Texture model complete';
save(file_out, 'mass_model');

%
%Calculate weights for combined model
%%%%%%%%%%%%%%%%%%%%%%
k_shape = length(L_shape);
k_tex   = length(L_tex);

W_shape = k_shape / sum(sqrt(L_shape));
W_tex   = k_tex / sum(sqrt(L_tex));
W_scale = 1 / sqrt(L_scale); %length L_scale = 1
% W_n = 1 / sqrt(L_n); %length L_n = 1

% W_shape = k_shape / sum(L_shape);
% W_tex   = k_tex / sum(L_tex);
% W_scale = 1 / L_scale; %length(L_scale = 1)
%W_n = 1 / L_n; %length L_n = 1

mass_model.W_shape = W_shape;
mass_model.W_tex = W_tex;
mass_model.W_scale = W_scale;
% mass_model.W_n = W_n;

%combined_data = [W_shape*B_shape; W_tex*B_tex; W_scale*B_scale; W_n*B_n]';
combined_data = [W_shape*B_shape; W_tex*B_tex; W_scale*B_scale]';
mass_model.combined_data = combined_data;

[mean_c, P_c, B_c, L_c] = pca(combined_data, 0.98);%, 0);

mass_model.mean_c = mean_c;
mass_model.P_c = P_c;
mass_model.B_c = B_c;
mass_model.L_c = L_c;

display('7: Combined model complete, function successful!');
mass_model.progress = '7: Combined model complete, function successful!';
save(file_out, 'mass_model');

