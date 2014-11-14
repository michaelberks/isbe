% GET_SHAPES_FROM_MASSES put mass outlines into shape matrix 
%    [shapes] = get_shapes_from_masses(mass_files, ...)
%
%    inputs:
%       mass_files  - Structure listing file names of input masses
%                       each DxN matrix corresponds to a particular shape
%       mass_path   - File path to mass folder
%       n_pts       - number of points to use per shape          
%                       D is the dimensionality of the data points
%       type        - {'standard'} - shape matrix is (N x 2*n_pts) 
%                      'mdl' - shape matrix is (2 x n_pts x N) 
%
%    outputs:
%       shapes      - the matrix of shapes
%
%    notes:
%
% Copyright: (C) 2006-2008 University of Manchester
% Author: Michael Berks

function [W_shape, W_scale, W_tex] = ...
    calculate_weights(model, mass_files, path, file_out)

    
    m = load(model);
    
    mean_shape  = m.mass_model.mean_shape;
    P_shape     = m.mass_model.P_shape;
    B_shape     = m.mass_model.B_shape;
    mean_scale  = m.mass_model.mean_scale;
    P_scale     = m.mass_model.P_scale;
    B_scale     = m.mass_model.B_scale;
    mean_tex    = m.mass_model.mean_tex;
    P_tex       = m.mass_model.P_tex;
    B_tex       = m.mass_model.B_tex;

    mean_target = m.mass_model.mean_target;
    mean_shape_pl = m.mass_model.mean_shape_pl;
    size_shape_vec = length(mean_shape) / 2;
    scaled_shape = mean_shape / m.mass_model.scale_factor;
    
    X_shape     = m.mass_model.X_shape;
    mean_dilate = m.mass_model.mean_dilate;
    
    clear m;
    
    
    %Define points to be interpolated by TPS - as row vectors
    i_x = mean_shape_pl(:,1)';
    i_y = mean_shape_pl(:,2)';
    clear mean_shape_pl;
    
    %Define source points for TPS - as row vectors
    s_x = mean_dilate(1:size_shape_vec);
    s_y = mean_dilate(size_shape_vec+1:end);

    clear mean_dilate;
    
    tps_L_inv = tps_weights(s_x, s_y);

    M = size(B_shape, 1);
    O = size(B_tex, 1);
        
    N = length(mass_files);
    RMS_shape = zeros(N, M);
    
    RMS_tex = zeros(N, O);
    RMS_scale = zeros(N, 2);
    
    for ii = 1:N
        
        display(['Calculate weight for mass ', num2str(ii)]);
%         f1 = figure('WindowStyle', 'docked'); hold on;
        
        temp = load([path, mass_files(ii).name]);
        [rows cols] = size(temp.mass.mass_ROI); 
        mass_centroid = temp.mass.mass_centroid; 
        subtract_ROI = temp.mass.subtract_ROI;
        clear temp;
        
        x_shape = reshape(X_shape(ii,:), [], 2);
        [dd Z trans] = procrustes(reshape(mean_target,[],2), x_shape);
        
        b_shape = B_shape(:,ii);
        x_shape_old = P_shape*b_shape + scaled_shape';
        x_shape_old = reshape(x_shape_old, [], 2);
        x_shape_old = (x_shape_old-trans.c)*inv(trans.T)/trans.b;
        
%         figure(f1);
%         plot(x_shape(:,1), x_shape(:,2));
%         plot(x_shape_old(:,1), x_shape_old(:,2), 'r:');

        shape_bw = roipoly(rows, cols,...
            x_shape_old(:,1) + mass_centroid(1),...
            x_shape_old(:,2) + mass_centroid(2));

        for kk = 1:49; 
            shape_bw = imdilate(shape_bw, strel('disk', 1)); 
        end

        [start_r start_c] = find(shape_bw, 1);
        dilate_outline = bwtraceboundary(shape_bw, [start_r, start_c], 'E');
        clear shape_bw;
        dilate_outline = [dilate_outline(:,2) dilate_outline(:,1)];

        px = x_shape_old(1,1) + mass_centroid(1);
        py = x_shape_old(1,2) + mass_centroid(2);
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
        old_tex = interp2(subtract_ROI, f_x, f_y);
        old_tex(isnan(old_tex)) = 0;
        clear f_x f_y z_x z_y
%         temp = zeros(mean_row, mean_col);
%         temp(sub2ind([mean_row mean_col],...
%             mean_shape_pl(:,1), mean_shape_pl(:,2))) = old_tex;
% 
%         figure('WindowStyle', 'docked'); 
%         imagesc(temp); axis image; colormap(gray(256)); clear temp;
        
        %
        % Compute RMS change in texture per shape mode
        %%%%%%%%%%%%%%%%%%%%%%%
        for jj = [1:M]
            b_shapei = b_shape;
            b_shapei(jj) = b_shapei(jj) + 10;
            
            x_shape_new = P_shape*b_shapei + scaled_shape';
            x_shape_new = reshape(x_shape_new, [], 2);
            x_shape_new = (x_shape_new-trans.c)*inv(trans.T)/trans.b;

%             figure(f1);
%             plot(x_shape_new(:,1), x_shape_new(:,2), 'g:');
            
            shape_bw = roipoly(rows, cols,...
                x_shape_new(:,1) + mass_centroid(1),...
                x_shape_new(:,2) + mass_centroid(2));

            for kk = 1:49; 
                shape_bw = imdilate(shape_bw, strel('disk', 1)); 
            end

            [start_r start_c] = find(shape_bw, 1);
            dilate_outline = bwtraceboundary(shape_bw, [start_r, start_c], 'E');
            clear shape_bw;
            dilate_outline = [dilate_outline(:,2) dilate_outline(:,1)];

            px = x_shape_new(1,1) + mass_centroid(1);
            py = x_shape_new(1,2) + mass_centroid(2);
            dd = sum((dilate_outline - repmat([px py], size(dilate_outline,1), 1)).^2, 2);
            [mm idx] = min(dd);

            dilate_outline = [dilate_outline(idx:end, :); dilate_outline(1:idx-1,:)];
            idx = round(linspace(1, size(dilate_outline,1), size_shape_vec+1));
            idx(end) = [];
            x_dilate = dilate_outline(idx,:);
            
%             figure(f1);
%             plot(x_dilate(:,1) - mass_centroid(1),...
%                 x_dilate(:,2) - mass_centroid(2), 'k:'); axis image;
            
            %Define displacement to target points
            z_x = x_dilate(:,1)';
            z_y = x_dilate(:,2)';
            clear x_dilate dilate_outline px py dd mm idx
            
            display(['sampling texture from shape mode ', num2str(jj)]);
            %Compute displacement of interpolated points
            f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
            f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
            
            %Get texture vector
            new_tex = interp2(subtract_ROI, f_x, f_y);
            new_tex(isnan(new_tex)) = 0;
            
%             temp = zeros(mean_row, mean_col);
%             temp(sub2ind([mean_row mean_col],...
%                 mean_shape_pl(:,1), mean_shape_pl(:,2))) = new_tex;
%             
%             figure('WindowStyle', 'docked'); 
%             imagesc(temp); axis image; colormap(gray(256)); clear temp;
            
            RMS_shape(ii, jj) = sqrt(sum((new_tex-old_tex).^2)) / 10;
            clear new_tex f_x f_y z_x z_y 
        end
%         axis image;
        
        
        %Compute RMS change in texture per scale mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        display('sampling texture from scale mode');
        x_scale_old = P_scale*B_scale(ii) + mean_scale';
        
        % + 0.1 units
        %%%%%%%%%%%
        b_scalei = B_scale(ii) + 0.1;
        
        x_scale_new = P_scale*b_scalei + mean_scale';
        scale_diff = x_scale_new / x_scale_old;
        x_shape_new = x_shape_old * scale_diff;
        
        shape_bw = roipoly(rows, cols,...
            x_shape_new(:,1) + mass_centroid(1),...
            x_shape_new(:,2) + mass_centroid(2));

        for kk = 1:49; 
            shape_bw = imdilate(shape_bw, strel('disk', 1)); 
        end

        [start_r start_c] = find(shape_bw, 1);
        dilate_outline = bwtraceboundary(shape_bw, [start_r, start_c], 'E');
        clear shape_bw;
        dilate_outline = [dilate_outline(:,2) dilate_outline(:,1)];

        px = x_shape_new(1,1) + mass_centroid(1);
        py = x_shape_new(1,2) + mass_centroid(2);
        dd = sum((dilate_outline - repmat([px py], size(dilate_outline,1), 1)).^2, 2);
        [mm idx] = min(dd);

        dilate_outline = [dilate_outline(idx:end, :); dilate_outline(1:idx-1,:)];
        idx = round(linspace(1, size(dilate_outline,1), size_shape_vec+1));
        idx(end) = [];
        x_dilate = dilate_outline(idx,:);

        %Define displacement to target points
        z_x = x_dilate(:,1)';
        z_y = x_dilate(:,2)';
        clear x_dilate dilate_outline px py dd mm idx

        %Compute displacement of interpolated points
        f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
        f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);

        %Get texture vector
        new_tex = interp2(subtract_ROI, f_x, f_y);
        new_tex(isnan(new_tex)) = 0;

        RMS_scale(ii, 1) = sqrt(sum((new_tex-old_tex).^2)) / 0.1;
        clear new_tex f_x f_y z_x z_y
        
        % - 0.1 units
        %%%%%%%%%%%
        b_scalei = B_scale(ii) - 0.1;
        
        x_scale_new = P_scale*b_scalei + mean_scale';
        scale_diff = x_scale_new / x_scale_old;
        x_shape_new = x_shape_old * scale_diff;

        shape_bw = roipoly(rows, cols,...
            x_shape_new(:,1) + mass_centroid(1),...
            x_shape_new(:,2) + mass_centroid(2));

        for kk = 1:49; 
            shape_bw = imdilate(shape_bw, strel('disk', 1)); 
        end

        [start_r start_c] = find(shape_bw, 1);
        dilate_outline = bwtraceboundary(shape_bw, [start_r, start_c], 'E');
        clear shape_bw;
        dilate_outline = [dilate_outline(:,2) dilate_outline(:,1)];

        px = x_shape_new(1,1) + mass_centroid(1);
        py = x_shape_new(1,2) + mass_centroid(2);
        dd = sum((dilate_outline - repmat([px py], size(dilate_outline,1), 1)).^2, 2);
        [mm idx] = min(dd);

        dilate_outline = [dilate_outline(idx:end, :); dilate_outline(1:idx-1,:)];
        idx = round(linspace(1, size(dilate_outline,1), size_shape_vec+1));
        idx(end) = [];
        x_dilate = dilate_outline(idx,:);

        %Define displacement to target points
        z_x = x_dilate(:,1)';
        z_y = x_dilate(:,2)';
        clear x_dilate dilate_outline px py dd mm idx

        %Compute displacement of interpolated points
        f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
        f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);

        %Get texture vector
        new_tex = interp2(subtract_ROI, f_x, f_y);
        new_tex(isnan(new_tex)) = 0;

        RMS_scale(ii, 2) = sqrt(sum((new_tex-old_tex).^2)) / 0.1;
        clear new_tex f_x f_y z_x z_y
        
        % 
        % Calculate texture weights
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        b_tex = B_tex(:,ii);
        old_tex = P_tex*b_tex + mean_tex';
        for jj = 1:O
            b_texi = b_tex;
            b_texi(jj) = b_texi(jj) + 1;
            new_tex = P_tex*b_texi + mean_tex';
            RMS_tex(ii, jj) = sqrt(sum((new_tex-old_tex).^2));
        end
        
        if nargin > 3
            save(file_out, 'RMS_scale', 'RMS_shape', 'RMS_tex');
        end
    end
    
    W_shape = mean(RMS_shape);
    W_scale = mean(RMS_scale);
    W_tex = mean(RMS_tex);       
    
end

