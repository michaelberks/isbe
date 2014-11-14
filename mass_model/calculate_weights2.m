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

function [RMS_shape, RMS_scale] = ...
    calculate_weights2(mass_files, offsets_shape, offsets_scale, varargin)

    args = u_packargs(varargin, 0,...
        'model', 'C:\isbe\dev\models\model_i500_50k.mat',...
        'mass_path', 'C:\isbe\dev\masses\',...
        'file_out', 'C:\isbe\dev\weights\fd_weights.mat',...
        'shape_mode', 1,...
        'indices', [],...
        'if_plot', 0);

    m = load(args.model);
    
    mean_shape  = m.mass_model.mean_shape;
    P_shape     = m.mass_model.P_shape;
    B_shape     = m.mass_model.B_shape;
    mean_scale  = m.mass_model.mean_scale;
    P_scale     = m.mass_model.P_scale;
    B_scale     = m.mass_model.B_scale;
    
    mean_shape_pl = m.mass_model.mean_shape_pl;
    size_shape_vec = length(mean_shape) / 2;
    scaled_mean_shape = m.mass_model.scaled_mean_shape;
    X_scale    = m.mass_model.X_scale;
    rotations = m.mass_model.rotations;
    translations = m.mass_model.translations;
    origins = m.mass_model.origins;
    
    clear m;
    
    %Define points to be interpolated by TPS - as row vectors
    i_x = mean_shape_pl(:,1)';
    i_y = mean_shape_pl(:,2)';
    clear mean_shape_pl;
    
    %Define source points for TPS - as row vectors
    s_x = scaled_mean_shape(1:size_shape_vec);% + mean_centre(1);
    s_y = scaled_mean_shape(size_shape_vec+1:end);% + mean_centre(2);
    
%     tps_L_inv = tps_weights(s_x, s_y);
    if isempty(args.indices)
        N = length(mass_files);
        args.indices = 1:N;
    else
        N = length(args.indices);
    end
    M = length(offsets_shape);
    O = length(offsets_scale);    
    
    
    RMS_shape = zeros(N, M);
    RMS_scale = zeros(N, O);
    
    for ii = 1:N
        
        display(['Calculate weight for mass ', num2str(ii)]);
        
        kk = args.indices(ii);
        temp = load([args.mass_path, mass_files(kk).name]); 
        subtract_ROI = temp.mass.subtract_ROI;
        clear temp;
               
        x_scale = X_scale(kk, :);
        rotation = rotations(:,:,kk);
        translation = repmat(translations(kk,:), size_shape_vec, 1);
        origin = origins(kk);
        
        b_shape = B_shape(:,kk);
        x_shape_old = P_shape*b_shape + mean_shape';
        x_shape_old = reshape(x_shape_old, [], 2);
        %x_shape_old = (x_shape_old-trans.c)*inv(trans.T)/trans.b;
        
        t_shape_old = (circshift(x_shape_old/x_scale, 1-origin)...
            - translation)*inv(rotation);
        %
        if args.if_plot
            figure; hold on; axis equal;
            plot(t_shape_old(:,1), t_shape_old(:,2), 'r:');
        end

        %Define displacement to target points
        z_x = t_shape_old(1:size_shape_vec);
        z_y = t_shape_old(size_shape_vec+1:end);    

        %Compute displacement of interpolated points
        T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
                'transform', 'spline');
        [pts] = geom_transformpoints([i_x; i_y], T);
        f_x = pts(1,:);
        f_y = pts(2,:);
        
        old_tex = interp2(subtract_ROI, f_x, f_y);
        old_tex(isnan(old_tex)) = 0;
        clear f_x f_y z_x z_y
        
        %
        % Compute RMS change in texture per shape mode
        %%%%%%%%%%%%%%%%%%%%%%%
        for jj = 1:M
            b_shapei = b_shape;
            b_shapei(args.shape_mode) = b_shapei(args.shape_mode) + offsets_shape(jj);
            
            x_shape_new = P_shape*b_shapei + mean_shape';
            x_shape_new = reshape(x_shape_new, [], 2);
            %x_shape_new = (x_shape_new-trans.c)*inv(trans.T)/trans.b;

            t_shape_new = (circshift(x_shape_new/x_scale, 1-origin)...
                - translation)*inv(rotation);
            
            if args.if_plot
                plot(t_shape_new(:,1), t_shape_new(:,2), 'g:');
            end

            
            %Define displacement to target points
            z_x = t_shape_new(1:size_shape_vec);
            z_y = t_shape_new(size_shape_vec+1:end);
            
            display(['sampling texture from shape offset ', num2str(jj)]);
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
            
            RMS_shape(ii, jj) = sqrt(sum((new_tex-old_tex).^2));
            clear new_tex f_x f_y z_x z_y 
        end
%         axis image;
        
        
        %Compute RMS change in texture per scale mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for jj = 1:O
            display(['sampling texture from scale offset ', num2str(jj)]);
            %x_scale_old = P_scale*B_scale(kk) + mean_scale';

            % + scale offsets units
            %%%%%%%%%%%
            b_scalei = B_scale(kk) + offsets_scale(jj);

            x_scale_new = P_scale*b_scalei + mean_scale';
            t_shape_new = (circshift(x_shape_old/x_scale_new, 1-origin)...
                - translation)*inv(rotation);

            if args.if_plot
                plot(t_shape_new(:,1), t_shape_new(:,2), 'c:');
                legend({'model shape', 'displaced shape', 'displaced scale'});
            end
            
            %Define displacement to target points
            z_x = t_shape_new(1:size_shape_vec);
            z_y = t_shape_new(size_shape_vec+1:end);

            %Compute displacement of interpolated points        
            T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
                    'transform', 'spline');
            [pts] = geom_transformpoints([i_x; i_y], T);
            f_x = pts(1,:);
            f_y = pts(2,:);

            %Get texture vector
            new_tex = interp2(subtract_ROI, f_x, f_y);
            new_tex(isnan(new_tex)) = 0;

            RMS_scale(ii, jj) = sqrt(sum((new_tex-old_tex).^2));
            clear new_tex f_x f_y z_x z_y
        end                  
        save(args.file_out, 'RMS_scale', 'RMS_shape');
               
    end     
    
end

