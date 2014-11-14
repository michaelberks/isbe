% LEAVE_ONE_OUT_MODEL
%    [shape_error tex_error com_error] = leave_one_out_model(mass_model, mass, weights)
%
%    inputs:
%       mass_model  - 
%       mass        -
%
%       optional:
%       weights     - 
%
%    outputs:
%
%    notes:
%
% Copyright: (C) 2006-2008 University of Manchester
% Author: Michael Berks

function [i_error c_error new_mass] =...
    leave_one_out_model(mass_model, mass, idx, weights)

    if nargin > 3
        mass_model = combine_model(mass_model, weights);
    end
    size_shape_vec = length(mass_model.mean_shape) / 2;
    
    t_shape = mass_model.shapes_unaligned(idx,:);
    x_shape = mass_model.X_shape(idx, :);
    x_scale = mass_model.X_scale(idx, :);
    x_tex = mass_model.X_tex(idx,:);
    rotation = mass_model.rotations(:,:,idx);
    translation = repmat(mass_model.translations(idx,:), size_shape_vec, 1);
    origin = mass_model.origins(idx);
    
    mass_model.X_shape(idx,:) = [];
    mass_model.X_scale(idx) = [];
    mass_model.X_tex(idx,:) = [];

    [mean_shape, P_shape, B_shape] = pca(mass_model.X_shape, 0.98);
    [mean_scale, P_scale, B_scale] = pca(mass_model.X_scale, 0.98);      
    [mean_tex, P_tex, B_tex] = pca(mass_model.X_tex, 0.98);%, 0);

    mean_shape_pl = mass_model.mean_shape_pl;
    scaled_mean_shape = mass_model.scaled_mean_shape;
    mean_row = mass_model.mean_row;
    mean_col = mass_model.mean_col;
    
    W_shape     = mass_model.W_shape;
    W_tex       = mass_model.W_tex;
    W_scale     = mass_model.W_scale;
    clear mass_model;

    k_shape     = size(P_shape, 2);
    k_tex       = size(P_tex, 2);

    combined_data = [W_shape*B_shape; W_tex*B_tex; W_scale*B_scale]';
    [mean_com, P_com] = pca(combined_data, 0.98);


    %Define points to be interpolated by TPS - as row vectors
    i_x = mean_shape_pl(:,1)';
    i_y = mean_shape_pl(:,2)';

    %Define source points for TPS - as row vectors
    s_x = scaled_mean_shape(1:size_shape_vec);% + mean_centre(1);
    s_y = scaled_mean_shape(size_shape_vec+1:end);% + mean_centre(2);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %masses
    [rows cols] = size(mass.subtract_ROI);
    
    b_shape = P_shape' * (x_shape-mean_shape)';
    b_tex = P_tex' * (x_tex-mean_tex)';
    b_scale = P_scale' * (x_scale-mean_scale)';
    x_c = [W_shape*b_shape; W_tex*b_tex; W_scale*b_scale]';
    b_c = P_com' * (x_c - mean_com)';

    x_shape_old = P_shape*b_shape + mean_shape';
    x_scale_old = P_scale*b_scale + mean_scale';
    x_tex_old = P_tex*b_tex + mean_tex';

    Q_shape = P_com(1:k_shape,:); 
    Q_tex = P_com(k_shape+1:k_shape + k_tex,:);
    Q_scale = P_com(end, :);

    x_shape_new = mean_shape' + (P_shape*inv(W_shape)*Q_shape*b_c);
    x_tex_new = mean_tex' + (P_tex*inv(W_tex)*Q_tex*b_c);
    x_scale_new = mean_scale' + (P_scale*Q_scale*b_c) / W_scale;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_shape_old = reshape(x_shape_old, size_shape_vec, 2);
    t_shape_old = (circshift(t_shape_old/x_scale_old, 1-origin) - translation)...
        *inv(rotation);
    old_bw = roipoly(mass.subtract_ROI, t_shape_old(:,1), t_shape_old(:,2));
    
    [old_shape_pl(:,2) old_shape_pl(:,1)]  = find(old_bw);
    clear old_bw;

    %Define displacement to target points
    z_x = t_shape_old(:,1)';
    z_y = t_shape_old(:,2)';

    T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
        'transform', 'spline');
    [pts] = geom_transformpoints([i_x; i_y], T);
    f_x = pts(1,:);
    f_y = pts(2,:);

    display('completed warping');
    % Create old shape pixel list
    %%%%%%%%%%

    old_shape_tex = griddata(f_x, f_y, x_tex_old,...
        old_shape_pl(:,1), old_shape_pl(:,2));
    old_shape_tex(isnan(old_shape_tex)) = 0;
    old_shape_ROI = zeros(rows, cols);
    old_shape_ROI(sub2ind([rows cols], old_shape_pl(:,2), old_shape_pl(:,1)))...
        = uint8(old_shape_tex);
    clear old_shape_pl

    display('completed grid'); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_shape_new = reshape(x_shape_new, size_shape_vec, 2);
    t_shape_new = (circshift(t_shape_new/x_scale_new, 1-origin) - translation)...
        *inv(rotation);
    new_bw = roipoly(mass.subtract_ROI, t_shape_new(:,1), t_shape_new(:,2));
    [new_shape_pl(:,2) new_shape_pl(:,1)]  = find(new_bw);
    clear new_bw;

    % Compute TPS warp to map from mean to new shape
    %%%%

    %Define displacement to target points
    z_x = t_shape_new(:,1)';
    z_y = t_shape_new(:,2)';

    %Compute displacement of interpolated points        
    T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
        'transform', 'spline');
    [pts] = geom_transformpoints([i_x; i_y], T);

    display('completed warping');
    % Create new shape pixel list
    %%%%%%%%%%

    new_shape_tex = griddata(pts(1,:), pts(2,:), x_tex_new,...
        new_shape_pl(:,1), new_shape_pl(:,2));
    new_shape_tex(isnan(new_shape_tex)) = 0;
    new_shape_ROI = zeros(rows, cols);
    new_shape_ROI(sub2ind([rows cols], new_shape_pl(:,2), new_shape_pl(:,1)))...
        = uint8(new_shape_tex);
    clear new_shape_pl

    display('completed grid'); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    er_shape = x_shape_old - x_shape';

    i_error.shape = mean(sqrt(er_shape(1:size_shape_vec).^2 ...
        + er_shape(size_shape_vec+1:end).^2)); clear er_shape;

    i_error.tex = mean(abs(x_tex_old - x_tex'));

    i_error.total = sum(abs(double(mass.subtract_ROI(:)) - ...
        double(old_shape_ROI(:))))  / (rows*cols);

    er_shape = x_shape_new - x_shape_old;

    c_error.shape = mean(sqrt(er_shape(1:size_shape_vec).^2 ...
        + er_shape(size_shape_vec+1:end).^2)); clear er_shape;

    c_error.tex = mean(abs(x_tex_new - x_tex_old));

    c_error.scale = mean(abs(x_scale_new - x_scale_old));

    c_error.weights = sum(abs(double(old_shape_ROI(:)) - ...
        double(new_shape_ROI(:))))  / (rows*cols);

    c_error.total = sum(abs(double(mass.subtract_ROI(:)) - ...
        double(new_shape_ROI(:))))  / (rows*cols);

    display(['Joint individual error is ', num2str(i_error.total)]);
    display(['Joint combined error is ', num2str(c_error.total)]);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    %plot effect on shapes
%     figure('WindowStyle', 'docked'); hold on;
%     title('Effect of model on shape');
%     plot(x_shape(1:size_shape_vec), x_shape(size_shape_vec+1:end), 'b');
%     plot(x_shape_old(1:size_shape_vec),...
%         x_shape_old(size_shape_vec+1:end), 'g:');
%     plot(x_shape_new(1:size_shape_vec),...
%         x_shape_new(size_shape_vec+1:end), 'r:');
%     legend('Original mass','Regenerated from shape model',...
%             'Regenerated from combined model');
%     axis image;
    
%     figure('WindowStyle', 'docked'); hold on;
%     title('Effect of model on unaligned shape');
%     plot(t_shape(1:size_shape_vec), t_shape(size_shape_vec+1:end), 'b');
%     plot(t_shape_old(1:size_shape_vec), t_shape_old(size_shape_vec+1:end), 'g:');
%     plot(t_shape_new(1:size_shape_vec), t_shape_new(size_shape_vec+1:end), 'r:');
%     legend('Original mass','Regenerated from shape model',...
%             'Regenerated from combined model');
%     axis image;
    
    %show effect on texture vectors
%     figure;
%     new_mass.orig_tex = zeros(mean_row, mean_col);
%     new_mass.orig_tex(sub2ind([mean_row, mean_col], mean_shape_pl(:,2), mean_shape_pl(:,1))) = x_tex;
%     new_mass.orig_tex(new_mass.orig_tex < 0) = 0;
%     subplot(1,3,1);
%     imagesc(new_mass.orig_tex); axis image; colormap(gray(256)); hold on;
%     clear temp;

%     new_mass.ind_tex = zeros(mean_row, mean_col);
%     new_mass.ind_tex(sub2ind([mean_row, mean_col], mean_shape_pl(:,2), mean_shape_pl(:,1))) = x_tex_old;
%     new_mass.ind_tex(new_mass.ind_tex < 0) = 0;
%     subplot(1,3,2);
%     imagesc(new_mass.ind_tex); axis image; colormap(gray(256)); hold on;
%     title('Effect of model on texture');
%     clear temp;
%     
%     new_mass.com_tex = zeros(mean_row, mean_col);
%     new_mass.com_tex(sub2ind([mean_row, mean_col], mean_shape_pl(:,2), mean_shape_pl(:,1))) = x_tex_new;
%     new_mass.com_tex(new_mass.com_tex < 0) = 0;
%     subplot(1,3,3);
%     imagesc(new_mass.com_tex); axis image; colormap(gray(256)); hold on;
%     clear temp;
    
    
%     %show effect on scale
%     figure('WindowStyle', 'docked'); hold on;
%     title('Effect of model on scale');
%     plot(x_shape(1:size_shape_vec)*x_scale_old,...
%         x_shape(size_shape_vec+1:end)*x_scale_old, 'g:');
%     plot(x_shape(1:size_shape_vec)*x_scale_old,...
%         x_shape(size_shape_vec+1:end)*x_scale_new, 'r:');
%     legend('Shape * orginal scale',...
%         'Shape * scale from combined model');
%     axis image;
    
    
    %show all three ROIs
    figure;
    subplot(1,3,1);
    imagesc(mass.subtract_ROI); axis image; colormap(gray(256));
    
    subplot(1,3,2);
    imagesc(old_shape_ROI); axis image; colormap(gray(256));
    title('Effect of model on texture');
    
    subplot(1,3,3);
    imagesc(new_shape_ROI); axis image; colormap(gray(256));

    new_mass.orig_shape = x_shape;
    new_mass.ind_shape = x_shape_old;
    new_mass.com_shape = x_shape_new;
    
    new_mass.orig_scale = x_scale_old;
    new_mass.com_scale = x_scale_new;
    
    new_mass.orig_roi = mass.subtract_ROI;
    new_mass.ind_roi = old_shape_ROI;
    new_mass.com_roi = new_shape_ROI;
    
%     % show difference maps
%     figure('WindowStyle', 'docked');
%     subplot(1,3,1);
%     imagesc(abs(double(mass.subtract_ROI)-double(old_shape_ROI))); axis image; colormap(jet(256));
%     subplot(1,3,2);
%     imagesc(abs(double(new_shape_ROI)-double(old_shape_ROI))); axis image; colormap(jet(256));
%     subplot(1,3,3);
%     imagesc(abs(double(mass.subtract_ROI)-double(new_shape_ROI))); axis image; colormap(jet(256));
           