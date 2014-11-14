mass_areas = zeros(179,1);
for ii = 1:179
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    mass_areas(ii) = mass.mass_area;
end
[sl idx] = sort(mass_areas);

%%
load C:\isbe\dev\masses\m_files.mat
load C:\isbe\dev\mass_model

mass_model = combine_model(mass_model);
%%
mean_shape  = mass_model.mean_shape;
mean_shape_pl = mass_model.mean_shape_pl;
size_shape_vec = length(mean_shape) / 2;
mean_off_r  = mass_model.mean_off_r;
mean_off_c  = mass_model.mean_off_c;

mean_row = ceil(max(mean_shape(size_shape_vec+1:2*size_shape_vec)))...
    + mean_off_r + 100;
mean_col = ceil(max(mean_shape(1:size_shape_vec)))...
    + mean_off_c + 100;

mean_bw = zeros(mean_row, mean_col);
mean_bw(sub2ind(size(mean_bw), mean_shape_pl(:,1), mean_shape_pl(:,2))) = 1;

[start_r start_c] = find(mean_bw, 1);
dilate_outline = bwtraceboundary(mean_bw, [start_r, start_c], 'E');
dilate_outline = [dilate_outline(:,2) dilate_outline(:,1)];

px = mean_shape(1) + mean_off_c;
py = mean_shape(size_shape_vec+1) + mean_off_r;
dd = sum((dilate_outline - repmat([px py], size(dilate_outline,1), 1)).^2, 2);
[mm idx] = min(dd);

dilate_outline = [dilate_outline(idx:end, :); dilate_outline(1:idx-1,:)];

idx = round(linspace(1, size(dilate_outline,1), size_shape_vec+1));
idx(end) = [];

mean_dilate = [dilate_outline(idx,1)', dilate_outline(idx,2)'];
%%
%
% figure; imagesc(mean_bw); axis image; colormap(gray(256)); hold on;
% plot(mean_shape(1:500)+mean_off_c, mean_shape(501:end)+mean_off_r, 'b.')
% plot(x_dilate(1:500), x_dilate(501:1000), 'b.')
% plot(start_c, start_r, 'gx')
% plot(mean_shape(1)+mean_off_c, mean_shape(501)+mean_off_r+49, 'rx')
% plot(x_dilate(1), x_dilate(501), 'gx')
%

%Define points to be interpolated by TPS - as row vectors
i_x = mean_shape_pl(:,1)';
i_y = mean_shape_pl(:,2)';
%
%Define source points for TPS - as row vectors
% s_x = [mean_shape(1:size_shape_vec) + mean_off_c, ...
%     mean_dilate(1:size_shape_vec)];
% s_y = [mean_shape(size_shape_vec+1:end) + mean_off_r, ...
%     mean_dilate(size_shape_vec+1:end)];
%
%%
s_x = [mean_shape(1:size_shape_vec) + mean_off_c];
s_y = [mean_shape(size_shape_vec+1:end) + mean_off_r];

% s_x = [mean_dilate(1:size_shape_vec)];
% s_y = [mean_dilate(size_shape_vec+1:end)];

tps_L_inv = tps_weights(s_x, s_y);

for ii = 1:5
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    idx = round(linspace(1, length(mass.mass_outline(:,1)), size_shape_vec+1));
    idx = idx(1:size_shape_vec); %ensures first point is not equal to the last point!
    x_shape = mass.mass_outline(idx,:);
    
    shape_bw = zeros(size(mass.mass_ROI));
    shape_bw(mass.mass_list) = 1;

    [start_r start_c] = find(shape_bw, 1);
    dilate_outline = bwtraceboundary(shape_bw, [start_r, start_c], 'E');
    dilate_outline = [dilate_outline(:,2) dilate_outline(:,1)];

    px = x_shape(1,1) + mass.mass_centroid(1);
    py = x_shape(1,2) + mass.mass_centroid(2);
    dd = sum((dilate_outline - repmat([px py], size(dilate_outline,1), 1)).^2, 2);
    [mm idx] = min(dd);

    dilate_outline = [dilate_outline(idx:end, :); dilate_outline(1:idx-1,:)];

    idx = round(linspace(1, size(dilate_outline,1), size_shape_vec+1));
    idx(end) = [];

    x_dilate = dilate_outline(idx,:);
    
    %Define displacement to target points
%     z_x = [x_shape(:,1)' + mass.mass_centroid(1), x_dilate(:,1)'];
%     z_y = [x_shape(:,2)' + mass.mass_centroid(2), x_dilate(:,2)'];
    
%     %Define displacement to target points
%     z_x = [x_dilate(:,1)'];
%     z_y = [x_dilate(:,2)'];
    
    %Define displacement to target points
    z_x = [x_shape(:,1)' + mass.mass_centroid(1)];
    z_y = [x_shape(:,2)' + mass.mass_centroid(2)];

    %Compute displacement of interpolated points
    f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
    f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
    
    new_tex = interp2(mass.subtract_ROI, f_x, f_y);
    new_im = zeros(mean_row, mean_col);
    new_im(sub2ind(size(new_im), mean_shape_pl(:,1), mean_shape_pl(:,2)))...
        = new_tex;
    figure('WindowStyle', 'docked');
    imagesc(mass.subtract_ROI); axis image; colormap(gray(256)); hold on;
    plot(f_x(1:3:end), f_y(1:3:end), 'r.', 'MarkerSize', 1.0);
%     plot(mass.mass_outline(:,1) + mass.mass_centroid(1),...
%         mass.mass_outline(:,2) + mass.mass_centroid(2));
%     plot(mass.mass_outline(1,1) + mass.mass_centroid(1),...
%         mass.mass_outline(1,2) + mass.mass_centroid(2), 'gx');
    plot(x_dilate(:,1), x_dilate(:,2), 'yx');
%     plot(x_dilate(1,1), x_dilate(1,2), 'rx');
    clear mass temp
end
%%
for ii = 1:179
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    if strcmpi(m_files(ii).name(8), 'R')
        [rows cols] = size(mass.mass_ROI);
        [mlr mlc] = ind2sub([rows cols], mass.mass_list);
        mlc = cols + 1 - mlc;
        mass.mass_list = sub2ind([rows cols], mlr, mlc);
        save(['C:\isbe\dev\masses\', m_files(ii).name], 'mass');
        display(['flipped list in ', m_files(ii).name]);
    end
end
%%
for ii = 1:179
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    
    shape_bw = zeros(size(mass.mass_ROI));
    shape_bw(mass.mass_list) = 1;

    [start_r start_c] = find(shape_bw, 1);
    dilate_outline = bwtraceboundary(shape_bw, [start_r, start_c], 'E');
    clear shape_bw;
    dilate_outline = [dilate_outline(:,2) dilate_outline(:,1)];

    px = mass.mass_outline(1,1) + mass.mass_centroid(1);
    py = mass.mass_outline(1,2) + mass.mass_centroid(2);
    dd = sum((dilate_outline - repmat([px py], size(dilate_outline,1), 1)).^2, 2);
    [mm idx] = min(dd);

    dilate_outline = [dilate_outline(idx:end, :); dilate_outline(1:idx-1,:)];
    mass.dilate_outline = dilate_outline;
    save(['C:\isbe\dev\masses\', m_files(ii).name], 'mass');
end
%%
mass_model.mean_row = ceil(max(mass_model.mean_shape(201:end)))...
    + mass_model.mean_off_r + 100;
mass_model.mean_col = ceil(max(mass_model.mean_shape(1:200)))...
    + mass_model.mean_off_c + 100;

temp = zeros(mean_row, mean_col);
temp(sub2ind([mean_row mean_col],...
    mass_model.mean_shape_pl(:,1), mass_model.mean_shape_pl(:,2))) = mass_model.X_tex(5,:);
figure('WindowStyle', 'docked'); 
imagesc(temp); axis image; colormap(gray(256)); clear temp;
%%
for ii = 1:179
    if ~rem(ii+1, 20)
        temp = zeros(mean_row, mean_col);
        temp(sub2ind([mean_row mean_col],...
            mass_model.mean_shape_pl(:,1), mass_model.mean_shape_pl(:,2)))...
            = uint8(mass_model.X_tex(ii,:));
        figure('WindowStyle', 'docked');
        imagesc(temp); axis image; colormap(gray(256));
        clear temp;
    end
end
    