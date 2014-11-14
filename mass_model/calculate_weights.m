% CALCULATE_WEIGHTS calculate weights for the combined appearance model
%               using the technique described by Cootes et al (see notes)

%    [shapes] = calculate_weights(mass_files, ...)
%
%    inputs:
%       model       - 
%       mass_files  -
%       mass_path   - File path to mass folder
%       file_out    - 
%
%    outputs:
%       W_shape, W_tex, W_scale
%
%    notes:
%
% Copyright: (C) 2006-2008 University of Manchester
% Author: Michael Berks

function [W_shape, W_scale] = ...
    calculate_weights(model, mass_files, mass_path, file_out)

    warning('Use calculate_errors2 instead!!');
    
    m = load(model);
    
    mean_shape  = m.mass_model.mean_shape;
    P_shape     = m.mass_model.P_shape;
    B_shape     = m.mass_model.B_shape;
    mean_scale  = m.mass_model.mean_scale;
    P_scale     = m.mass_model.P_scale;
    B_scale     = m.mass_model.B_scale;

    mean_target = m.mass_model.mean_target;
    mean_shape_pl = m.mass_model.mean_shape_pl;
    size_shape_vec = length(mean_shape) / 2;
    scaled_mean_shape = m.mass_model.scaled_mean_shape;
%     mean_centre = m.mass_model.mean_centre;
    X_shape     = m.mass_model.X_shape;
    
    clear m;
    
    
    %Define points to be interpolated by TPS - as row vectors
    i_x = mean_shape_pl(:,1)';
    i_y = mean_shape_pl(:,2)';
    clear mean_shape_pl;
    
    %Define source points for TPS - as row vectors
    s_x = scaled_mean_shape(1:size_shape_vec);% + mean_centre(1);
    s_y = scaled_mean_shape(size_shape_vec+1:end);% + mean_centre(2);
    
    %tps_L_inv = tps_weights(s_x, s_y);

    M = size(B_shape, 1);
        
    N = length(mass_files);
    RMS_shape = zeros(N, M);
    RMS_scale = zeros(N, 2);
    
    for ii = 1:N
        
        display(['Calculate weight for mass ', num2str(ii)]);
%         f1 = figure('WindowStyle', 'docked'); hold on;
        
        temp = load([mass_path, mass_files(ii).name]); 
        mass_centroid = temp.mass.mass_centroid; 
        subtract_ROI = temp.mass.subtract_ROI;
        clear temp;
        
        x_shape = reshape(X_shape(ii,:), [], 2);
        [dd Z trans] = procrustes(reshape(mean_target,[],2), x_shape);
        
        b_shape = B_shape(:,ii);
        x_shape_old = P_shape*b_shape + mean_shape';
        x_shape_old = reshape(x_shape_old, [], 2);
        x_shape_old = (x_shape_old-trans.c)*inv(trans.T)/trans.b;
        
%         figure(f1);
%         plot(x_shape(:,1), x_shape(:,2));
%         plot(x_shape_old(:,1), x_shape_old(:,2), 'r:');

        %Define displacement to target points
        z_x = x_shape_old(1:size_shape_vec) + mass_centroid(1);
        z_y = x_shape_old(size_shape_vec+1:end) + mass_centroid(2);    

        %Compute displacement of interpolated points
%         f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
%         f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
        
        T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
                'transform', 'spline');
        [pts] = geom_transformpoints([i_x; i_y], T);
        f_x = pts(1,:);
        f_y = pts(2,:);
        
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
        for jj = 1:M
            b_shapei = b_shape;
            b_shapei(jj) = b_shapei(jj) + 10;
            
            x_shape_new = P_shape*b_shapei + mean_shape';
            x_shape_new = reshape(x_shape_new, [], 2);
            x_shape_new = (x_shape_new-trans.c)*inv(trans.T)/trans.b;

%             figure(f1);
%             plot(x_shape_new(:,1), x_shape_new(:,2), 'g:');

            
            %Define displacement to target points
            z_x = x_shape_new(1:size_shape_vec) + mass_centroid(1);
            z_y = x_shape_new(size_shape_vec+1:end) + mass_centroid(2);
            
            display(['sampling texture from shape mode ', num2str(jj)]);
            %Compute displacement of interpolated points
%           f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
%           f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
        
            T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
                    'transform', 'spline');
            [pts] = geom_transformpoints([i_x; i_y], T);
            f_x = pts(1,:);
            f_y = pts(2,:);
            
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

        %Define displacement to target points
        z_x = x_shape_new(1:size_shape_vec) + mass_centroid(1);
        z_y = x_shape_new(size_shape_vec+1:end) + mass_centroid(2);

        %Compute displacement of interpolated points
%       f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
%       f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
        
        T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
                'transform', 'spline');
        [pts] = geom_transformpoints([i_x; i_y], T);
        f_x = pts(1,:);
        f_y = pts(2,:);

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

        z_x = x_shape_new(1:size_shape_vec) + mass_centroid(1);
        z_y = x_shape_new(size_shape_vec+1:end) + mass_centroid(2);

        %Compute displacement of interpolated points
%       f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
%       f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
        
        T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
                'transform', 'spline');
        [pts] = geom_transformpoints([i_x; i_y], T);
        f_x = pts(1,:);
        f_y = pts(2,:);

        %Get texture vector
        new_tex = interp2(subtract_ROI, f_x, f_y);
        new_tex(isnan(new_tex)) = 0;

        RMS_scale(ii, 2) = sqrt(sum((new_tex-old_tex).^2)) / 0.1;
        clear new_tex f_x f_y z_x z_y
        
              
        if nargin > 3
            save(file_out, 'RMS_scale', 'RMS_shape');
        end
    end
    
    W_shape = mean(RMS_shape);
    W_scale = mean(RMS_scale);    
    
end

