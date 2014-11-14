function [com_error ind_error] =...
    model_errors_loo(mass_files, weights, model_path, mass_path)
% Formerly model_errors3
    
    if nargin < 4
        mass_path = 'C:\isbe\dev\masses\';
    end
    N = length(mass_files);
    
    ind_error.shape = zeros(N,1);
    ind_error.tex = zeros(N,1);
    ind_error.total = zeros(N,1);
    
    com_error.shape = zeros(N,1);
    com_error.tex = zeros(N,1);
    com_error.scale = zeros(N,1);
    com_error.weights = zeros(N,1);
    com_error.total = zeros(N,1);
    
    for ii = 1:N
        model_name = [model_path, zerostr(ii, 3)];
        mass_name = [mass_path, mass_files(ii).name];
        [c_error i_error] = get_error(model_name, mass_name, weights);

        ind_error.shape(ii) = i_error.shape;
        ind_error.tex(ii) = i_error.tex;
        ind_error.total(ii) = i_error.total;
        
        com_error.shape(ii) = c_error.shape;
        com_error.tex(ii) = c_error.tex;
        com_error.scale(ii) = c_error.scale;
        com_error.weights(ii) = c_error.weights;
        com_error.total(ii) = c_error.total;
    end
    
    function [c_error i_error] =...
            get_error(model_name, mass_name, weights)
    
        m = load(model_name);
        if ~isempty(weights)
            m.mass_model = combine_model(m.mass_model, [weights, 1]);
        end

        mean_shape  = m.mass_model.mean_shape;
        P_shape     = m.mass_model.P_shape;
        mean_scale  = m.mass_model.mean_scale;
        P_scale     = m.mass_model.P_scale;
        mean_tex    = m.mass_model.mean_tex;
        P_tex       = m.mass_model.P_tex;
        mean_c      = m.mass_model.mean_c;
        P_c         = m.mass_model.P_c;
        W_shape     = m.mass_model.W_shape;
        W_tex       = m.mass_model.W_tex;
        W_scale     = m.mass_model.W_scale;

        mean_shape_pl = m.mass_model.mean_shape_pl;
        scale_factor = m.mass_model.scale_factor;
        mean_target = m.mass_model.mean_target;
        mean_dilate = m.mass_model.mean_dilate;
        clear m;

        size_shape_vec = length(mean_shape) / 2;
        k_shape     = size(P_shape, 2);
        k_tex       = size(P_tex, 2);

        scaled_mean = mean_shape / scale_factor;

        %Define points to be interpolated by TPS - as row vectors
        i_x = mean_shape_pl(:,1)';
        i_y = mean_shape_pl(:,2)';

        %Define source points for TPS - as row vectors
        s_x = mean_dilate(1:size_shape_vec);
        s_y = mean_dilate(size_shape_vec+1:end);

        tps_L_inv = tps_weights(s_x, s_y);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %individual masses
        
        temp = load(mass_name);
        mass = temp.mass; clear temp;
        [rows cols] = size(mass.subtract_ROI);
        
        idx = round(linspace(1, length(mass.mass_outline(:,1)), size_shape_vec+1));
        idx(end) = []; %ensures first point is not equal to the last point
        x_shape = mass.mass_outline(idx,:);
        [dd Z t] = procrustes(reshape(mean_target, size_shape_vec, 2), x_shape);

        x_pro = [Z(:,1)', Z(:,2)'];
        x_scale = t.b;

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

        Q_shape = P_c(1:k_shape,:); 
        Q_tex = P_c(k_shape+1:k_shape + k_tex,:);
        Q_scale = P_c(end, :);

        x_shape_new = scaled_mean' + (P_shape*Q_shape*b_c) / W_shape;
        x_tex_new = mean_tex' + (P_tex*Q_tex*b_c) / W_tex;
        x_scale_new = mean_scale' + (P_scale*Q_scale*b_c) / W_scale;

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
        er_shape = x_shape_old - x_pro';
        
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
        
        clear mass;
        
        save('com_error', 'com_error');
        display(['Joint individual error is ', num2str(i_error.total)]);
        display(['Joint combined error is ', num2str(c_error.total)]);
    
    end
end