function [shape_error tex_error com_error] =...
    leave_one_out_model(mass_model, mass, weights)

    if nargin > 2
        mass_model = combine_model(mass_model, [weights, 1]);
    end
    
    mean_shape  = mass_model.mean_shape;
    P_shape     = mass_model.P_shape;
    mean_scale  = mass_model.mean_scale;
    P_scale     = mass_model.P_scale;
    mean_tex    = mass_model.mean_tex;
    P_tex       = mass_model.P_tex;
    mean_c      = mass_model.mean_c;
    P_c         = mass_model.P_c;
    W_shape     = mass_model.W_shape;
    W_tex       = mass_model.W_tex;
    W_scale     = mass_model.W_scale;

    mean_shape_pl = mass_model.mean_shape_pl;
    size_shape_vec = length(mean_shape) / 2;
    k_shape     = size(P_shape, 2);
    k_tex       = size(P_tex, 2);

    scale_factor = mass_model.scale_factor;
    scaled_mean = mean_shape / scale_factor;
    mean_target = mass_model.mean_target;
    mean_row  = mass_model.mean_row;
    mean_col  = mass_model.mean_col;
    mean_dilate = mass_model.mean_dilate;
    [rows cols] = size(mass.subtract_ROI);
    
    idx = round(linspace(1, length(mass.mass_outline(:,1)), size_shape_vec+1));
    idx = idx(1:size_shape_vec); %ensures first point is not equal to the last point!
    x_shape = mass.mass_outline(idx,:);
    [dd Z t] = procrustes(reshape(mean_target, size_shape_vec, 2), x_shape);

    x_pro = [Z(:,1)', Z(:,2)'];
    x_scale = t.b;

    %Define points to be interpolated by TPS - as row vectors
    i_x = mean_shape_pl(:,1)';
    i_y = mean_shape_pl(:,2)';

    %Define source points for TPS - as row vectors
    s_x = mean_dilate(1:size_shape_vec);
    s_y = mean_dilate(size_shape_vec+1:end);

    tps_L_inv = tps_weights(s_x, s_y);

    idx = round(linspace(1, size(mass.dilate_outline,1), size_shape_vec+1));
    idx(end) = [];
    x_dilate = mass.dilate_outline(idx,:);
    clear idx

    %Define displacement to target points
    z_x = x_dilate(:,1)';
    z_y = x_dilate(:,2)';

    %Compute displacement of interpolated points
    f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
    f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);

%     figure('WindowStyle', 'docked');
%     imagesc(mass.subtract_ROI); axis image; colormap(gray(256)); hold on;
%     plot(f_x(1:3:end), f_y(1:3:end), 'r.', 'MarkerSize', 0.5);

    %Get texture vector
    x_tex = interp2(mass.subtract_ROI, f_x, f_y);

    x_tex(isnan(x_tex)) = 0;

    b_shape = P_shape' * (x_pro-scaled_mean)';
    b_tex = P_tex' * (x_tex-mean_tex)';
    b_scale = P_scale' * (x_scale-mean_scale)';
    x_c = [W_shape*b_shape; W_tex*b_tex; W_scale*b_scale]';
    b_c = P_c' * (x_c - mean_c)';

    x_shape_old = P_shape*b_shape + scaled_mean';
    x_scale_old = P_scale*b_scale + mean_scale';
    x_tex_old = P_tex*b_tex + mean_tex';
    
    shape_error = x_shape_old - x_pro';
    tex_error = x_tex_old - x_tex';
    
    Q_shape = P_c(1:k_shape,:); 
    Q_tex = P_c(k_shape+1:k_shape + k_tex,:);
    Q_scale = P_c(end, :);

    x_shape_new = scaled_mean' + (P_shape*Q_shape*b_c) / W_shape;
    x_tex_new = mean_tex' + (P_tex*Q_tex*b_c) / W_tex;
    x_scale_new = mean_scale' + (P_scale*Q_scale*b_c) / W_scale;

    com_error.shape = x_shape_new - x_shape_old;
    com_error.tex = x_tex_new - x_tex_old;
    com_error.scale = x_scale_new - x_scale_old;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_shape_old = reshape(x_shape_old, size_shape_vec, 2);
    t_shape_old = t_shape_old*inv(t.T) / x_scale_old;
    t_shape_old(:,1) = t_shape_old(:,1) + mass.mass_centroid(1);
    t_shape_old(:,2) = t_shape_old(:,2) + mass.mass_centroid(2);
    old_bw = roipoly(mass.subtract_ROI, t_shape_old(:,1), t_shape_old(:,2));
    
    for jj = 1:49; old_bw = imdilate(old_bw, strel('disk', 1)); end
    rp = regionprops(bwlabel(old_bw, 4), 'PixelList');
    old_shape_pl = rp.PixelList; clear rp;
    clear rp;

    % Compute TPS warp to map from mean to old shape
    %%%%

    [start_r start_c] = find(old_bw, 1);
    dilate_outline = bwtraceboundary(old_bw, [start_r, start_c], 'E');
    dilate_outline = [dilate_outline(:,2) dilate_outline(:,1)];
    clear old_bw;

    px = t_shape_old(1,1);
    py = t_shape_old(1,2);
    dd = sum((dilate_outline - repmat([px py], size(dilate_outline,1), 1)).^2, 2);
    [mm idx] = min(dd);

    dilate_outline = [dilate_outline(idx:end, :); dilate_outline(1:idx-1,:)];
    idx = round(linspace(1, size(dilate_outline,1), size_shape_vec+1));
    idx(end) = [];
    x_dilate = dilate_outline(idx,:);
    clear idx

    %Define displacement to target points
    z_x = x_dilate(:,1)';
    z_y = x_dilate(:,2)';
    clear x_dilate;

    %Compute displacement of interpolated points        
    f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
    f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);

    display('completed warping');
    % Create old shape pixel list
    %%%%%%%%%%

    old_shape_tex = griddata(f_x, f_y, x_tex_old,...
        old_shape_pl(:,1), old_shape_pl(:,2));
    old_shape_tex(isnan(old_shape_tex)) = 0;
    old_shape_ROI = zeros(rows, cols);
    old_shape_ROI(sub2ind([rows cols], old_shape_pl(:,2), old_shape_pl(:,1)))...
        = uint8(old_shape_tex);

    display('completed grid'); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_shape_new = reshape(x_shape_new, size_shape_vec, 2);
    t_shape_new = t_shape_new*inv(t.T) / x_scale_new;
    t_shape_new(:,1) = t_shape_new(:,1) + mass.mass_centroid(1);
    t_shape_new(:,2) = t_shape_new(:,2) + mass.mass_centroid(2);
    new_bw = roipoly(mass.subtract_ROI, t_shape_new(:,1), t_shape_new(:,2));
    
    for jj = 1:49; new_bw = imdilate(new_bw, strel('disk', 1)); end
    rp = regionprops(bwlabel(new_bw, 4), 'PixelList');
    new_shape_pl = rp.PixelList; clear rp;
    clear rp;

    % Compute TPS warp to map from mean to new shape
    %%%%

    [start_r start_c] = find(new_bw, 1);
    dilate_outline = bwtraceboundary(new_bw, [start_r, start_c], 'E');
    dilate_outline = [dilate_outline(:,2) dilate_outline(:,1)];
    clear new_bw;

    px = t_shape_new(1,1);
    py = t_shape_new(1,2);
    dd = sum((dilate_outline - repmat([px py], size(dilate_outline,1), 1)).^2, 2);
    [mm idx] = min(dd);

    dilate_outline = [dilate_outline(idx:end, :); dilate_outline(1:idx-1,:)];
    idx = round(linspace(1, size(dilate_outline,1), size_shape_vec+1));
    idx(end) = [];
    x_dilate = dilate_outline(idx,:);
    clear idx

    %Define displacement to target points
    z_x = x_dilate(:,1)';
    z_y = x_dilate(:,2)';
    clear x_dilate;

    %Compute displacement of interpolated points        
    f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
    f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);

    display('completed warping');
    % Create new shape pixel list
    %%%%%%%%%%

    new_shape_tex = griddata(f_x, f_y, x_tex_new,...
        new_shape_pl(:,1), new_shape_pl(:,2));
    new_shape_tex(isnan(new_shape_tex)) = 0;
    new_shape_ROI = zeros(rows, cols);
    new_shape_ROI(sub2ind([rows cols], new_shape_pl(:,2), new_shape_pl(:,1)))...
        = uint8(new_shape_tex);

    display('completed grid'); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
%     x_shape_new = reshape(x_shape_new, size_shape_vec, 2);
%     x_shape_new = x_shape_new*inv(t.T) / x_scale_new;
% 
%     x_shape_new(:,1) = x_shape_new(:,1) + mass.mass_centroid(1);
%     x_shape_new(:,2) = x_shape_new(:,2) + mass.mass_centroid(2);

%     b_shape_c = Q_shape*b_c / W_shape;
%     b_tex_c = Q_tex*b_c / W_tex;
%     b_scale_c = Q_scale*b_c / W_scale;
    

    figure('WindowStyle', 'docked'); hold on;
    title('Effect of model on shape');
    plot(x_pro(1:size_shape_vec), x_pro(size_shape_vec+1:end), 'b');
    plot(x_shape_old(1:size_shape_vec),...
        x_shape_old(size_shape_vec+1:end), 'g:');
    plot(x_shape_new(1:size_shape_vec),...
        x_shape_new(size_shape_vec+1:end), 'r:');
    legend('Original mass','Regenerated from shape model',...
            'Regenerated from combine model');
    axis image;
    
%     figure('WindowStyle', 'docked');
%     
%     temp1 = zeros(mean_row, mean_col);
%     temp1(sub2ind(size(temp1), mean_shape_pl(:,1), mean_shape_pl(:,2))) = x_tex;
%     temp1(temp1 < 0) = 0;
%     subplot(1,3,1);
%     imagesc(temp1); axis image; colormap(gray(256)); hold on;
%     clear temp1;

    temp2 = zeros(mean_row, mean_col);
    temp2(sub2ind(size(temp2), mean_shape_pl(:,1), mean_shape_pl(:,2))) = x_tex_old;
    temp2(temp2 < 0) = 0;
    subplot(1,3,2);
    imagesc(temp2); axis image; colormap(gray(256)); hold on;
    title('Effect of model on texture');
    clear temp2;
    
    temp3 = zeros(mean_row, mean_col);
    temp3(sub2ind(size(temp3), mean_shape_pl(:,1), mean_shape_pl(:,2))) = x_tex_new;
    temp3(temp3 < 0) = 0;
    subplot(1,3,3);
    imagesc(temp3); axis image; colormap(gray(256)); hold on;
    clear temp3;
    
    figure('WindowStyle', 'docked'); hold on;
    title('Effect of model on scale');
    plot(x_shape_old(1:size_shape_vec)*x_scale_old,...
        x_shape_old(size_shape_vec+1:end)*x_scale_old, 'g:');
    plot(x_shape_new(1:size_shape_vec)*x_scale_new,...
        x_shape_new(size_shape_vec+1:end)*x_scale_new, 'r:');
    legend('Shape * orginal scale',...
        'Shape * scale from combined model');
    axis image;    
    
%     figure('WindowStyle', 'docked');
%     subplot(1,3,1);
%     imagesc(mass.subtract_ROI); axis image; colormap(gray(256));
%     hold on;
%     plot(mass.mass_outline(:,1) + mass.mass_centroid(1),...
%         mass.mass_outline(:,2) + mass.mass_centroid(2), 'r:');
%     plot(t_shape_new(:,1), t_shape_new(:,2), 'g:');
%     subplot(1,3,2);
%     imagesc(old_shape_ROI); axis image; colormap(gray(256));
%     subplot(1,3,3);
%     imagesc(new_shape_ROI); axis image; colormap(gray(256));

    figure('WindowStyle', 'docked');
    subplot(1,3,1);
    imagesc(abs(double(mass.subtract_ROI)-double(old_shape_ROI))); axis image; colormap(jet(256));
    subplot(1,3,2);
    imagesc(abs(double(new_shape_ROI)-double(old_shape_ROI))); axis image; colormap(jet(256));
    subplot(1,3,3);
    imagesc(abs(double(mass.subtract_ROI)-double(new_shape_ROI))); axis image; colormap(jet(256));
    
    the_error1 = sum(abs(double(mass.subtract_ROI(:)) - ...
        double(old_shape_ROI(:))))  / length(x_tex_new);
    display(['The error is ', num2str(the_error1)]);
    
    the_error2 = sum(abs(double(old_shape_ROI(:)) - ...
        double(new_shape_ROI(:))))  / length(x_tex_new);
    display(['The error is ', num2str(the_error2)]);
        