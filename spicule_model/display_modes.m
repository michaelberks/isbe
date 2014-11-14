% DISPLAY_MODES - create .avi videos of variation in the principal modes of
% shape, texture and combined appearance
%
% Usage: display_modes(mass_model, path)
% 
%   Inputs:
%       mass_model  - structure of combined appearance model 
%       path        - file path of folder to save videos to
%       if_c        - {0}, 1 flag to choose whether we display the combined
%                       appearance mode variation (this takes considerably
%                       longer to compute)
% Outputs: 
%        
% Notes:
% 
%
% See also:
%
% References:
%
% Author:   Michael Berks
%           Imaging Science and Biomedical Engineering
%           University of Manchester
%

function display_modes(mass_model, path, if_c)

    
    if nargin < 3
        if_c = 0;
    end
        
    P_shape = mass_model.P_shape;
    L_shape = mass_model.L_shape;
    mean_shape = mass_model.mean_shape;
    mean_shape_pl = mass_model.mean_shape_pl;
    
    P_tex = mass_model.P_tex;
    L_tex = mass_model.L_tex;
    mean_tex = mass_model.mean_tex;
    
    P_scale = mass_model.P_scale;
%     L_scale = mass_model.L_scale;
    mean_scale = mass_model.mean_scale;
    
    P_c = mass_model.P_c;
    L_c = mass_model.L_c;
    mean_c = mass_model.mean_c;
    
    W_shape     = mass_model.W_shape;
    W_tex       = mass_model.W_tex;
    W_scale     = mass_model.W_scale;

    k_shape     = length(L_shape);
    k_tex       = length(L_tex);
    
    Q_shape = P_c(1:k_shape,:); 
    Q_tex = P_c(k_shape+1:k_shape + k_tex,:);
    Q_scale = P_c(end-1, :);
    
    scale_factor = mass_model.scale_factor;
    mean_off_r  = mass_model.mean_off_r;
    mean_off_c  = mass_model.mean_off_c;
    
    %Define source points for TPS - as row vectors
    s_x = mean_shape(1:end/2) + mean_off_c;
    s_y = mean_shape(end/2+1:end) + mean_off_r;

    %Define points to be interpolated by TPS - as row vectors
    i_x = mean_shape_pl(:,1)';
    i_y = mean_shape_pl(:,2)';

    tps_L_inv = tps_weights(s_x, s_y);
% 

    figure('Position', [100,100, 618, 758]);
    aviobj = avifile([path, 'shape_modes.avi'],'fps', 10);
    for ii = 1:20
        k = -2 + 4*ii/20;
        subplot(2,1,1);
        shape_ii = mean_shape + (k*P_shape(:,1)*sqrt(L_shape(1)))';
        plot([shape_ii(1:end/2) shape_ii(1)],...
         [shape_ii(1+end/2:end) shape_ii(1+end/2)]); axis equal;
        title('1st Mode of Shape');
        
        subplot(2,1,2);
        shape_ii = mean_shape + (k*P_shape(:,2)*sqrt(L_shape(2)))';
        plot([shape_ii(1:end/2) shape_ii(1)],...
         [shape_ii(1+end/2:end) shape_ii(1+end/2)]); axis equal;
        title('2nd Mode of Shape');
        frame = getframe(gcf);
        aviobj = addframe(aviobj,frame);
    end
    for ii = 1:20
        k = 2 - 4*ii/20;
        subplot(2,1,1);
        shape_ii = mean_shape + (k*P_shape(:,1)*sqrt(L_shape(1)))';
        plot([shape_ii(1:end/2) shape_ii(1)],...
         [shape_ii(1+end/2:end) shape_ii(1+end/2)]); axis equal;
        title('1st Mode of Shape');
        
        subplot(2,1,2);
        shape_ii = mean_shape + (k*P_shape(:,2)*sqrt(L_shape(2)))';
        plot([shape_ii(1:end/2) shape_ii(1)],...
         [shape_ii(1+end/2:end) shape_ii(1+end/2)]); axis equal;
        title('2st Mode of Shape');
        frame = getframe(gcf);
        aviobj = addframe(aviobj,frame);
    end
    aviobj = close(aviobj); clear aviobj;
    
    m_rows = ceil(max(mass_model.mean_shape(1+end/2:end)))...
         + mass_model.mean_off_r + 100;
    m_cols = ceil(max(mass_model.mean_shape(1:end/2)))...
         + mass_model.mean_off_c + 100;
    mean_idx = sub2ind([m_rows m_cols], mean_shape_pl(:,2), mean_shape_pl(:,1));
    
    figure('Position', [100,100, 618, 758]);
    aviobj = avifile([path, 'tex_modes.avi'], 'fps', 10);
    for ii = 1:20
        k = -2 + 4*ii/20;
        
        subplot(2,1,1);
        tex_ii = mean_tex + (k*P_tex(:,1)*sqrt(L_tex(1)))';
        im_ii = zeros(m_rows, m_cols); im_ii(mean_idx) = tex_ii;
        image(im_ii); axis image; colormap gray;
        title('1st Mode of Texture');
        
        subplot(2,1,2);
        tex_ii = mean_tex + (k*P_tex(:,2)*sqrt(L_tex(2)))';
        im_ii = zeros(m_rows, m_cols); im_ii(mean_idx) = tex_ii;
        image(im_ii); axis image; colormap gray;
        title('2nd Mode of Texture');
        frame = getframe(gcf);
        aviobj = addframe(aviobj,frame);
    end
    for ii = 1:20
        k = 2 - 4*ii/20;
        
        subplot(2,1,1);
        tex_ii = mean_tex + (k*P_tex(:,1)*sqrt(L_tex(1)))';
        im_ii = zeros(m_rows, m_cols); im_ii(mean_idx) = tex_ii;
        image(im_ii); axis image; colormap gray;
        title('1st Mode of Texture');
        
        subplot(2,1,2);
        tex_ii = mean_tex + (k*P_tex(:,2)*sqrt(L_tex(2)))';
        im_ii = zeros(m_rows, m_cols); im_ii(mean_idx) = tex_ii;
        image(im_ii); axis image; colormap gray;
        title('2nd Mode of Texture');
        frame = getframe(gcf);
        aviobj = addframe(aviobj,frame);
    end
    aviobj = close(aviobj); clear aviobj;
    
    if if_c
        mkdir([path, 'com_modes']);
        for ii = 1:20
            k = -2 + 4*ii/20;
            figure('Position', [100,100, 618, 758]);
            subplot(2,1,1);
            c_ii = ((mean_c + (k*P_c(:,1)*sqrt(L_c(1)))') * P_c)';
            [c_shape_ii c_ROI_ii] = make_mode_im(c_ii);
            image(c_ROI_ii); axis image; colormap gray; hold on;
            plot([c_shape_ii(1:end/2) c_shape_ii(1)],...
                  [c_shape_ii(1+end/2:end) c_shape_ii(1+end/2)], 'r');
            title('1st Combined Appearance Mode');

            subplot(2,1,2);
            c_ii = ((mean_c + (k*P_c(:,2)*sqrt(L_c(2)))') * P_c)';
            [c_shape_ii c_ROI_ii] = make_mode_im(c_ii);
            image(c_ROI_ii); axis image; colormap gray; hold on;
            plot([c_shape_ii(1:end/2) c_shape_ii(1)],...
                  [c_shape_ii(1+end/2:end) c_shape_ii(1+end/2)], 'r');
            title('2nd Combined Appearance Mode');
            saveas(gcf, [path, 'com_modes\', num2str(ii), '.fig']);
            close gcf
        end

        aviobj = avifile([path, 'com_modes.avi'], 'fps', 10);
        for ii = 1:20
            openfig([path, 'com_modes\', num2str(ii), '.fig']);
            frame = getframe(gcf);
            aviobj = addframe(aviobj,frame);
            close gcf
        end
        for ii = 1:20
            openfig([path, 'com_modes\', num2str(21-ii), '.fig']);
            frame = getframe(gcf);
            aviobj = addframe(aviobj,frame);
            close gcf
        end
        aviobj = close(aviobj); clear aviobj;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [c_shape c_ROI] = make_mode_im(c_vec)
        c_shape = mean_shape / scale_factor + (P_shape*Q_shape*c_vec)' / W_shape;
        c_tex = mean_tex + (P_tex*Q_tex*c_vec)' / W_tex;
        c_scale = mean_scale + (P_scale*Q_scale*c_vec)' / W_scale;    
        c_shape = c_shape*c_scale;
        
        off_r = 101 - round(min(c_shape(end/2+1:end)));
        off_c = 101 - round(min(c_shape(1:end/2)));

        c_shape(end/2+1:end) = c_shape(end/2+1:end) + off_r;
        c_shape(1:end/2) = c_shape(1:end/2) + off_c;

        bw_r = round(max(c_shape(end/2+1:end)))+ 100;
        bw_c = round(max(c_shape(1:end/2))) + 100;

        c_bw = roipoly(bw_r, bw_c, c_shape(1:end/2),...
            c_shape(end/2+1:end));
        
        rp = regionprops(bwlabel(c_bw, 4), 'PixelList');
        c_shape_pl = rp.PixelList; clear rp;
        
        %Define displacement to target points
        z_x = c_shape(1:end/2);
        z_y = c_shape(end/2+1:end);
 
        %Compute displacement of interpolated points        
        f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
        f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
         
        display('completed warping');
        c_shape_tex = griddata(f_x, f_y, c_tex,...
            c_shape_pl(:,1), c_shape_pl(:,2));
        
% Using griddata is much faster than original method below
%
%         N = length(c_shape_pl(:,1));
%         c_shape_tex = zeros(N, 1);
%         for jj = 1:N
%             pix_dist = (f_x - c_shape_pl(jj,1)).^2 + ...
%                 (f_y - c_shape_pl(jj,2)).^2;
%             [min_val, nearest_pix] = min(pix_dist);
%             c_shape_tex(jj) = c_tex(nearest_pix);
%         end

        clear new_tex;

        c_ROI = zeros(bw_r, bw_c);
        c_ROI(sub2ind([bw_r bw_c], c_shape_pl(:,2), c_shape_pl(:,1)))...
            = uint8(c_shape_tex);
    
    end         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end