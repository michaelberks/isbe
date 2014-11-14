% MODEL_ERRORS
%    [com_error shape_error tex_error] = model_errors(model, mass_files, mass_path, weights)
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
%    notes: % Formerly model_errors2
%
% Copyright: (C) 2006-2008 University of Manchester
% Author: Michael Berks

function [com_error ind_error] =...
    model_errors(model, mass_files, mass_path, weights)


    m = load(model);
    if nargin > 3
        m.mass_model = combine_model(m.mass_model, [weights, 1]);
        display('model combined');
    end
    if nargin < 3
        mass_path = 'C:\isbe\dev\masses\';
    end
    
    mean_shape  = m.mass_model.mean_shape;
    P_shape     = m.mass_model.P_shape;
    mean_scale  = m.mass_model.mean_scale;
    P_scale     = m.mass_model.P_scale;
    mean_tex    = m.mass_model.mean_tex;
    P_tex       = m.mass_model.P_tex;
    mean_com      = m.mass_model.mean_com;
    P_com         = m.mass_model.P_com;
    W_shape     = m.mass_model.W_shape;
    W_tex       = m.mass_model.W_tex;
    W_scale     = m.mass_model.W_scale;

    mean_shape_pl = m.mass_model.mean_shape_pl;
    scaled_mean_shape = mass_model.scaled_mean_shape;
    mean_target = m.mass_model.mean_target;
    clear m;
    
    size_shape_vec = length(mean_shape) / 2;
    k_shape     = size(P_shape, 2);
    k_tex       = size(P_tex, 2);
       
    %Define points to be interpolated by TPS - as row vectors
    i_x = mean_shape_pl(:,1)';
    i_y = mean_shape_pl(:,2)';

    %Define source points for TPS - as row vectors
    s_x = scaled_mean_shape(1:end/2);% + mean_centre(1);
    s_y = scaled_mean_shape(end/2+1:end);% + mean_centre(2);

    %tps_L_inv = tps_weights(s_x, s_y);
    
    N = length(mass_files);
    
    ind_error.shape = zeros(N,1);
    ind_error.tex = zeros(N,1);
    ind_error.total = zeros(N,1);
    
    com_error.shape = zeros(N,1);
    com_error.tex = zeros(N,1);
    com_error.scale = zeros(N,1);
    com_error.weights = zeros(N,1);
    com_error.total = zeros(N,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %individual masses
    for ii = 1:N
        
        temp = load([mass_path, mass_files(ii).name]);
        mass = temp.mass; clear temp;
        [rows cols] = size(mass.subtract_ROI);
        
        idx = round(linspace(1, length(mass.mass_outline(:,1)), size_shape_vec+1));
        idx = idx(1:size_shape_vec); %ensures first point is not equal to the last point!
        x_shape = mass.mass_outline(idx,:);
        [dd Z t] = procrustes(reshape(mean_target, size_shape_vec, 2), x_shape);

        x_pro = [Z(:,1)', Z(:,2)'];
        x_scale = t.b;

        %Define displacement to target points
        z_x = x_shape(:,1)';% + mass.mass_centroid(1);
        z_y = x_shape(:,2)';% + mass.mass_centroid(2);

        %Compute displacement of interpolated points
%         f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
%         f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
        T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
            'transform', 'spline');
        [pts] = geom_transformpoints([i_x; i_y], T);
        f_x = pts(1,:);
        f_y = pts(2,:);
        %Get texture vector
        x_tex = interp2(mass.subtract_ROI, f_x, f_y);
        x_tex(isnan(x_tex)) = 0;

        b_shape = P_shape' * (x_pro-mean_shape)';
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
        t_shape_old = t_shape_old*inv(t.T) / x_scale_old;
        t_shape_old(:,1) = t_shape_old(:,1);% + mass.mass_centroid(1);
        t_shape_old(:,2) = t_shape_old(:,2);% + mass.mass_centroid(2);
        old_bw = roipoly(mass.subtract_ROI, t_shape_old(:,1), t_shape_old(:,2));

%         for jj = 1:49; old_bw = imdilate(old_bw, strel('disk', 1)); end
        [old_shape_pl(:,2) old_shape_pl(:,1)] = find(old_bw);
        clear old_bw;
        
        % Compute TPS warp to map from mean to old shape
        %%%%

        %Define displacement to target points
        z_x = t_shape_old(:,1)';
        z_y = t_shape_old(:,2)';
        clear x_dilate;

        %Compute displacement of interpolated points        
%       f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
%       f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
        
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
        t_shape_new = t_shape_new*inv(t.T) / x_scale_new;
        t_shape_new(:,1) = t_shape_new(:,1);% + mass.mass_centroid(1);
        t_shape_new(:,2) = t_shape_new(:,2);% + mass.mass_centroid(2);
        new_bw = roipoly(mass.subtract_ROI, t_shape_new(:,1), t_shape_new(:,2));

%         for jj = 1:49; new_bw = imdilate(new_bw, strel('disk', 1)); end
        [new_shape_pl(:,2) new_shape_pl(:,1)] = find(new_bw);
        clear new_bw;

        % Compute TPS warp to map from mean to new shape
        %%%%

        %Define displacement to target points
        z_x = t_shape_new(:,1)';
        z_y = t_shape_new(:,2)';

        %Compute displacement of interpolated points        
%         f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
%         f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
        
        T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
                'transform', 'spline');
        [pts] = geom_transformpoints([i_x; i_y], T);
        f_x = pts(1,:);
        f_y = pts(2,:);
        
        display('completed warping');
        % Create new shape pixel list
        %%%%%%%%%%

        new_shape_tex = griddata(f_x, f_y, x_tex_new,...
            new_shape_pl(:,1), new_shape_pl(:,2));
        new_shape_tex(isnan(new_shape_tex)) = 0;
        new_shape_ROI = zeros(rows, cols);
        new_shape_ROI(sub2ind([rows cols], new_shape_pl(:,2), new_shape_pl(:,1)))...
            = uint8(new_shape_tex);
        clear new_shape_pl
        
        display('completed grid'); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculate errors
        
        er_shape = x_shape_old - x_pro';
        ind_error.shape(ii) = mean(sqrt(er_shape(1:size_shape_vec).^2 ...
            + er_shape(size_shape_vec+1:end).^2)); clear er_shape;
        ind_error.tex(ii) = mean(abs(x_tex_old - x_tex'));
        
        er_shape = x_shape_new - x_shape_old;
        com_error.shape(ii) = mean(sqrt(er_shape(1:size_shape_vec).^2 ...
            + er_shape(size_shape_vec+1:end).^2)); clear er_shape;
        com_error.tex(ii) = mean(abs(x_tex_new - x_tex_old));
        com_error.scale(ii) = mean(abs(x_scale_new - x_scale_old));
        
        ind_error.total(ii) = sum(abs(double(mass.subtract_ROI(:)) - ...
            double(old_shape_ROI(:))))  / (rows*cols);
        com_error.weights(ii) = sum(abs(double(old_shape_ROI(:)) - ...
            double(new_shape_ROI(:))))  / (rows*cols);
        com_error.total(ii) = sum(abs(double(mass.subtract_ROI(:)) - ...
            double(new_shape_ROI(:))))  / (rows*cols);
        
        save('com_error', 'com_error');
        display(['Joint individual error is ', num2str(ind_error.total(ii))]);
        display(['Joint combined error is ', num2str(com_error.total(ii))]);
    end
        