% DISPLAY_MODES - create .avi videos of variation in the principal modes of
% shape, texture and combined appearance
%
% Usage: display_modes(mass_model, mass_path)
% 
%   Inputs:
%       mass_model  - structure of combined appearance model 
%       mass_path        - file mass_path of folder to save videos to
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

function display_modes2(mass_model, video_path, if_c, mode)

    
    if nargin < 3
        if_c = 0;
    end
    if nargin < 4
        mode = 1;
    end
        
    P_shape = mass_model.P_shape;
    L_shape = mass_model.L_shape;
    mean_shape = mass_model.mean_shape;
    mean_shape_pl = mass_model.mean_shape_pl;
    
    P_tex = mass_model.P_tex;
    L_tex = mass_model.L_tex;
    mean_tex = mass_model.mean_tex;
    
    P_scale = mass_model.P_scale;
    mean_scale = mass_model.mean_scale;
    
    P_com = mass_model.P_com;
    L_com = mass_model.L_com;
    mean_com = mass_model.mean_com;
    
    W_shape     = mass_model.W_shape;
    W_tex       = mass_model.W_tex;
    W_scale     = mass_model.W_scale;

    k_shape     = length(L_shape);
    k_tex       = length(L_tex);
    
    Q_shape = P_com(1:k_shape,:); 
    Q_tex = P_com(k_shape+1:k_shape + k_tex,:);
    Q_scale = P_com(end-1, :);
    
    scaled_mean_shape = mass_model.scaled_mean_shape;
    
    %Define source points for TPS - as row vectors
    s_x = scaled_mean_shape(1:end/2);% + mean_centre(1);
    s_y = scaled_mean_shape(end/2+1:end);% + mean_centre(2);

    %Define points to be interpolated by TPS - as row vectors
    i_x = mean_shape_pl(:,1)';
    i_y = mean_shape_pl(:,2)';

%     tps_L_inv = tps_weights(s_x, s_y);
    
    mkdir([video_path, 'shape_modes']);
    figure('Position', [100,100, 400, 400], 'WindowStyle', 'normal', 'Color', [0, 0, 0]);
    
    %aviobj = avifile([video_path, 'shape_modes.avi'],'fps', 10);
    for ii = 1:20
        k = -2 + 4*ii/20;
        shape_ii = mean_shape + (k*P_shape(:,mode)*sqrt(L_shape(mode)))';
        plot([shape_ii(1:end/2) shape_ii(1)],...
         [shape_ii(1+end/2:end) shape_ii(1+end/2)], 'g', 'LineWidth', 1.5); 
        axis equal; axis([-200 200 -200 200]); axis off;
        
        saveas(gcf, [video_path, 'shape_modes\shape', zerostr(ii, 2), '.fig']);
        shape_frames(ii) = getframe(gcf); %#ok
    end
%     for ii = 1:20
%         k = 2 - 4*ii/20;
%         shape_ii = mean_shape + (k*P_shape(:,mode)*sqrt(L_shape(mode)))';
%         plot([shape_ii(1:end/2) shape_ii(1)],...
%          [shape_ii(1+end/2:end) shape_ii(1+end/2)], 'g', 'LineWidth', 1.5);
%         axis equal; axis([-200 200 -200 200]); axis off;
%         
%         shape_frames(20+ii) = getframe(gcf); %#ok
%     end
    close gcf;
    shape_frames = [shape_frames, fliplr(shape_frames)];
    movie2avi(shape_frames,...
        [video_path, 'shape_modes', zerostr(mode, 2), '.avi'],'compression','none','fps',10)
    save([video_path, 'shape_frames', zerostr(mode, 2)], 'shape_frames');
    
    m_rows = mass_model.mean_row;
    m_cols = mass_model.mean_col;
    mean_idx = sub2ind([m_rows m_cols], mean_shape_pl(:,2), mean_shape_pl(:,1));
    
    mkdir([video_path, 'tex_modes']);
    figure('Position', [100,100, 600, 600], 'WindowStyle', 'normal', 'Color', [0, 0, 0]);
    for ii = 1:20
        k = -2 + 4*ii/20;
        
        tex_ii = mean_tex + (k*P_tex(:,mode)*sqrt(L_tex(mode)))';
        im_i1 = zeros(m_rows, m_cols); im_i1(mean_idx) = tex_ii;
        imagesc(im_i1); caxis([0 80]); axis image; colormap(gray(256)); axis off;
        
        saveas(gcf, [video_path, 'tex_modes\tex', zerostr(ii, 2), '.fig']);
        tex_frames(ii) = getframe(gcf); %#ok
    end
%     for ii = 1:20
%         k = 2 - 4*ii/20;
%        
%         tex_ii = mean_tex + (k*P_tex(:,mode)*sqrt(L_tex(mode)))';
%         im_ii = zeros(m_rows, m_cols); im_ii(mean_idx) = tex_ii;
%         imagesc(im_ii); caxis([0 80]); axis image; colormap(gray(256)); axis off;
% 
%         tex_frames(20+ii) = getframe(gcf); %#ok
%     end
    close gcf;
    tex_frames = [tex_frames, fliplr(tex_frames)];
    movie2avi(tex_frames,[video_path, 'tex_modes', zerostr(mode, 2), '.avi'],'compression','none','fps',10)
    save([video_path, 'tex_frames', zerostr(mode, 2)], 'tex_frames'); 
    if if_c
        mkdir([video_path, 'com_modes']);
        figure('Position', [100,100, 600, 600], 'WindowStyle', 'normal', 'Color', [0, 0, 0]);
        for ii = 1:20
            k = -2 + 4*ii/20;
            
            c_ii = ((mean_com + (k*P_com(:,mode)*sqrt(L_com(mode)))') * P_com)';
            [c_shape_ii c_ROI_i1] = make_mode_im(c_ii);
            imagesc(c_ROI_i1); caxis([0 80]); caxis([0 80]); axis image; colormap(gray(256)); hold on;
            plot([c_shape_ii(1:end/2) c_shape_ii(1)],...
                  [c_shape_ii(1+end/2:end) c_shape_ii(1+end/2)], 'r');
            axis off;
            
            saveas(gcf, [video_path, 'com_modes\com', zerostr(ii, 2), '.fig']);
            com_frames(ii) = getframe(gcf); %#ok
            %save([video_path, 'com_modes\com_ROI', zerostr(ii, 2)], 'c_ROI_i1');
            hold off;
        end
        close gcf
%         for ii = 1:20
%             openfig([video_path, 'com_modes\com', zerostr(ii, 2), '.fig']);
%             com_frames(ii) = getframe(gcf); %#ok
%             close gcf
%         end
%         for ii = 1:20
%             openfig([video_path, 'com_modes\com', zerostr(21-ii, 2), '.fig']);
%             com_frames(20+ii) = getframe(gcf); %#ok
%             close gcf
%         end
        com_frames = [com_frames, fliplr(com_frames)];
        movie2avi(com_frames,[video_path, 'com_modes', zerostr(mode, 2), '.avi'],'compression','none','fps',10)
        save([video_path, 'com_frames', zerostr(mode, 2)], 'com_frames');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [c_shape c_ROI] = make_mode_im(c_vec)
        c_shape = mean_shape + (P_shape*inv(W_shape)*Q_shape*c_vec)';
        c_tex = mean_tex + (P_tex*inv(W_tex)*Q_tex*c_vec)';
        c_scale = mean_scale + (P_scale*Q_scale*c_vec)' / W_scale;    
        c_shape = c_shape*c_scale;
        
        off_r = 51 - round(min(c_shape(end/2+1:end)));
        off_c = 51 - round(min(c_shape(1:end/2)));

        c_shape(end/2+1:end) = c_shape(end/2+1:end) + off_r;
        c_shape(1:end/2) = c_shape(1:end/2) + off_c;

        bw_r = round(max(c_shape(end/2+1:end)))+ 50;
        bw_c = round(max(c_shape(1:end/2))) + 50;

        c_bw = roipoly(bw_r, bw_c, c_shape(1:end/2),...
            c_shape(end/2+1:end));
        
        rp = regionprops(bwlabel(c_bw, 4), 'PixelList');
        c_shape_pl = rp.PixelList; clear rp;
        
        %Define displacement to target points
        z_x = c_shape(1:end/2);
        z_y = c_shape(end/2+1:end);
 
        %Compute displacement of interpolated points        
%         f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
%         f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
        T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
                'transform', 'spline');
        [pts] = geom_transformpoints([i_x; i_y], T);
        f_x = pts(1,:);
        f_y = pts(2,:);
         
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