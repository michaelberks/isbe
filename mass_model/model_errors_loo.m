% MODEL_ERRORS_LOO
%    [com_error ind_error] =... 
%       model_errors_loo(mass_files, weights, model_path, mass_path)
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
%    notes: % Formerly model_errors3
%
% Copyright: (C) 2006-2008 University of Manchester
% Author: Michael Berks

function [com_error ind_error] =...
    model_errors_loo(mass_files, model_path, varargin)
    
    args = u_packargs(varargin, 0,...
        'mass_path', 'C:\isbe\dev\masses\',...
        'indices', [],...
        'weights', [],...
        'shiftOrigin', 1, ...
        'plot', 1);
    
    if isempty(args.indices)
        N = length(mass_files);
        args.indices = 1:N;
    else
        N = length(args.indices);
    end
    
    ind_error.shape = zeros(N,1);
    ind_error.tex = zeros(N,1);
    ind_error.total = zeros(N,1);
    
    com_error.shape = zeros(N,1);
    com_error.tex = zeros(N,1);
    com_error.scale = zeros(N,1);
    com_error.weights = zeros(N,1);
    com_error.total = zeros(N,1);
    
    for ii = 1:N
        jj = args.indices(ii);
        model_name = [model_path, zerostr(jj, 3)];
        mass_name = [args.mass_path, mass_files(jj).name];
        [c_error i_error] =...
            get_error(model_name, mass_name, args.weights, args.shiftOrigin, args.plot);

        ind_error.shape(ii) = i_error.shape;
        ind_error.tex(ii) = i_error.tex;
        ind_error.total(ii) = i_error.total;
        
        com_error.shape(ii) = c_error.shape;
        com_error.tex(ii) = c_error.tex;
        com_error.scale(ii) = c_error.scale;
        com_error.weights(ii) = c_error.weights;
        com_error.total(ii) = c_error.total;
    end
end    
    function [c_error i_error] =...
            get_error(model_name, mass_name, weights, shiftOrigin, if_plot)
    
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
        mean_com      = m.mass_model.mean_com;
        P_com         = m.mass_model.P_com;
        W_shape     = m.mass_model.W_shape;
        W_tex       = m.mass_model.W_tex;
        W_scale     = m.mass_model.W_scale;

        mean_shape_pl = m.mass_model.mean_shape_pl;
        mean_target = m.mass_model.mean_target;
        scaled_mean_shape = m.mass_model.scaled_mean_shape;
        clear m;

        size_shape_vec = length(mean_shape) / 2;
        k_shape     = size(P_shape, 2);
        k_tex       = size(P_tex, 2);
        

        %Define points to be interpolated by TPS - as row vectors
        i_x = mean_shape_pl(:,1)';
        i_y = mean_shape_pl(:,2)';

        %Define source points for TPS - as row vectors
        s_x = scaled_mean_shape(1:size_shape_vec);% + mean_centre(1);
        s_y = scaled_mean_shape(size_shape_vec+1:end);% + mean_centre(2);

%         tps_L_inv = tps_weights(s_x, s_y);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %individual masses
        temp = load(mass_name);
        mass = temp.mass; clear temp;
        [rows cols] = size(mass.subtract_ROI);
        
        idx = round(linspace(1, length(mass.mass_outline(:,1)), size_shape_vec+1));
        idx(end) = []; %ensures first point is not equal to the last point
        shape_u = mass.mass_outline(idx,:);
        if shiftOrigin
            % cycle around the shape, shifting the origin each time
            % we choose the one with the minimal SSD to the target
            ssd = zeros(size_shape_vec, 1); 
            for origin = 1 : size_shape_vec
                shift_shape = circshift(shape_u, origin - 1); 
                ssd (origin) = procrustes(reshape(mean_target, [], 2), shift_shape);
            end
            [dummy, origin] = min (ssd);
            % apply the optimal origin shift to source shape...
            shift_shape = circshift(shape_u, origin - 1);
        end

        [dd Z t] = procrustes(reshape(mean_target, [], 2), shift_shape);

        x_shape = [Z(:,1)', Z(:,2)'];
        x_scale = t.b;

        %Define displacement to target points
        z_x = shape_u(:,1)';% + mass.mass_centroid(1);
        z_y = shape_u(:,2)';% + mass.mass_centroid(2);

        %Compute displacement of interpolated points
%         f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
%         f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);

        T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
            'transform', 'spline');
        [pts] = geom_transformpoints([i_x; i_y], T);

        %Get texture vector
        x_tex = interp2(mass.subtract_ROI, pts(1,:), pts(2,:));
        x_tex(isnan(x_tex)) = 0;

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
        t_shape_old = (circshift(t_shape_old/ x_scale_old, 1-origin) - t.c)*inv(t.T);
        old_bw = roipoly(mass.subtract_ROI, t_shape_old(:,1), t_shape_old(:,2));

%         for jj = 1:49; old_bw = imdilate(old_bw, strel('disk', 1)); end
        [old_shape_pl(:,2) old_shape_pl(:,1)]  = find(old_bw);
        clear old_bw;

        %Define displacement to target points
        z_x = t_shape_old(:,1)';
        z_y = t_shape_old(:,2)';

        %Compute displacement of interpolated points        
%         f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
%         f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);

        T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
            'transform', 'spline');
        [pts] = geom_transformpoints([i_x; i_y], T);

        display('completed warping');
        % Create old shape pixel list
        %%%%%%%%%%

        old_shape_tex = griddata(pts(1,:), pts(2,:), x_tex_old,...
            old_shape_pl(:,1), old_shape_pl(:,2));
        old_shape_tex(isnan(old_shape_tex)) = 0;
        old_shape_ROI = zeros(rows, cols);
        old_shape_ROI(sub2ind([rows cols], old_shape_pl(:,2), old_shape_pl(:,1)))...
            = uint8(old_shape_tex);
        clear old_shape_pl
        
        display('completed grid'); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t_shape_new = reshape(x_shape_new, size_shape_vec, 2);
        t_shape_new = (circshift(t_shape_new/x_scale_new, 1-origin) - t.c)*inv(t.T);
        new_bw = roipoly(mass.subtract_ROI, t_shape_new(:,1), t_shape_new(:,2));

%         for jj = 1:49; new_bw = imdilate(new_bw, strel('disk', 1)); end
        [new_shape_pl(:,2) new_shape_pl(:,1)]  = find(new_bw);
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
        
        i_error.shape = sqrt(mean(er_shape(1:size_shape_vec).^2 ...
            + er_shape(size_shape_vec+1:end).^2)); clear er_shape;
        
        i_error.tex = sqrt(mean((x_tex_old - x_tex').^2));
        
        i_error.total = sqrt(mean((double(mass.subtract_ROI(:)) - ...
            double(old_shape_ROI(:))).^2));
        
        er_shape = x_shape_new - x_shape_old;
        
        c_error.shape = sqrt(mean(er_shape(1:size_shape_vec).^2 ...
            + er_shape(size_shape_vec+1:end).^2)); clear er_shape;
        
        c_error.tex = sqrt(mean(x_tex_new - x_tex_old));
        
        c_error.scale = sqrt(mean(x_scale_new - x_scale_old));
        
        c_error.weights = sqrt(mean((double(old_shape_ROI(:)) - ...
            double(new_shape_ROI(:))).^2));
        
        c_error.total = sqrt(mean((double(mass.subtract_ROI(:)) - ...
            double(new_shape_ROI(:))).^2));
        
        clear mass;
        
        display(['Joint individual error is ', num2str(i_error.total)]);
        display(['Joint combined error is ', num2str(c_error.total)]);
        
        if if_plot
            figure;
            subplot(1,3,1);
            imagesc(mass.subtract_ROI); axis image; colormap(gray(256));

            subplot(1,3,2);
            imagesc(old_shape_ROI); axis image; colormap(gray(256));
            title('Effect of model on texture');

            subplot(1,3,3);
            imagesc(new_shape_ROI); axis image; colormap(gray(256));
        end
    
    end
