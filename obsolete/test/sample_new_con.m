function [new_masses] = sample_new_con(mass_model, no_of_masses, if_plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% author:   Michael Berks
% date:     08/06/2006  09:52
%
% function: generate new masses using mass appearance model
%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Get args from mass_model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_shape  = mass_model.mean_shape;
P_shape     = mass_model.P_shape;
L_shape     = mass_model.L_shape;
mean_scale  = mass_model.mean_scale;
P_scale     = mass_model.P_scale;
%L_scale     = mass_model.L_scale;
mean_tex    = mass_model.mean_tex;
P_tex       = mass_model.P_tex;
L_tex       = mass_model.L_tex;
P_ring      = mass_model.P_ring;
%L_ring      = mass_model.L_ring;
mean_ring   = mass_model.mean_ring;
%mean_c      = mass_model.mean_c;
P_c         = mass_model.P_c;
L_c         = mass_model.L_c;

W_shape     = mass_model.W_shape;
%W_tex       = mass_model.W_tex;
W_scale     = mass_model.W_scale;

C_com       = mass_model.C_com;

k_shape     = length(L_shape);
k_tex       = length(L_tex);
%k_ring      = length(L_ring);

Q_shape = P_c(1:k_shape,:); 
Q_scale = P_c(end, :);

%
% get C_mm, C_nm
C_com_inv = pinv(C_com);
C_mm = C_com_inv(1:k_tex, 1:k_tex);
C_nm = C_com_inv(1:k_tex, k_tex+1:end); clear C_com C_com_inv;
C_mm_inv = pinv(C_mm); clear c_mm;

mean_shape_pl = mass_model.mean_shape_pl;
mean_ring_pl = mass_model.mean_ring_pl;

size_shape_vec = length(mean_shape) / 2;

scale_factor = mass_model.scale_factor;
mean_off_r  = mass_model.mean_off_r;
mean_off_c  = mass_model.mean_off_c;

bad_shape = 0;

if if_plot,
    f1 = figure;
    f2 = figure;
end

%Define source points for TPS - as row vectors
s_x = mean_shape(1:size_shape_vec) + mean_off_c;
s_y = mean_shape(size_shape_vec+1:2*size_shape_vec) + mean_off_r;

%Define points to be interpolated by TPS - as row vectors
i_x = mean_shape_pl(:,1)';
i_y = mean_shape_pl(:,2)';

i_xr = mean_ring_pl(:,1)';
i_yr = mean_ring_pl(:,2)';
        
tps_L_inv = tps_weights(s_x, s_y);

%
% Generate specified number of masses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_masses(no_of_masses) = struct('mass_outline',[],'mass_ROI',[],...
    'B_tex',[],'B_shape',[], 'B_ring',[]);
ii = 1;
while ii <= no_of_masses

    % sample new combined appearance vectors - assume normal distribution
    % of modes
    %%%
    
    %
    % Sample new shape and vector - check mean?
    B_c = (randn(length(L_c), 1) .* sqrt(L_c));
    B_shape = Q_shape*B_c / W_shape;
    B_scale = Q_scale*B_c / W_scale;
    
    new_shape = mean_shape / scale_factor + ...
        (P_shape*B_shape)';
    new_scale = mean_scale + (P_scale*B_scale)';
    
    %
    % Reconstruct new shape including offsets
    new_shape = new_shape*new_scale;
    
    %off_r = 101 - round(min(new_shape(size_shape_vec+1:2*size_shape_vec)));
    %off_c = 101 - round(min(new_shape(1:size_shape_vec)));
    off_r = 300; off_c = 300;
    
    new_shape(size_shape_vec+1:2*size_shape_vec) = ...
        new_shape(size_shape_vec+1:2*size_shape_vec) + off_r;
    new_shape(1:size_shape_vec) = ...
        new_shape(1:size_shape_vec) + off_c;
    
    %x = new_shape(1:size_shape_vec);
    %y = new_shape(size_shape_vec+1:end);
    %new_masses(ii).dif = (max(x) - min(x) - max(y) + min(y))/2;
    
    % Compute pixel list for new shape
    %%%%
    %bw_r = round(max(new_shape(size_shape_vec+1:2*size_shape_vec)))+ 100;
    %bw_c = round(max(new_shape(1:size_shape_vec))) + 100;
    bw_r = 600; bw_c = 600;
    
    new_bw = roipoly(bw_r, bw_c, new_shape(1:size_shape_vec),...
        new_shape(size_shape_vec+1:2*size_shape_vec));
        
    rp = regionprops(bwlabel(new_bw, 4), 'PixelList');
    rp_l = size(rp, 1);
    
    if rp_l == 1, 
        new_shape_pl = rp.PixelList; clear rp;
        
        % Compute pixel list for ROI ring
        %%%%
        mask1 = new_bw;
        for jj = 1:40
            mask1 = imdilate(mask1, strel('disk', 1));
        end

        %new_ring = mask1 - new_bw; clear new_bw mask1;
        %rp = regionprops(bwlabel(new_ring, 4), 'PixelList'); clear new_ring;
        %new_ring_pl = rp.PixelList; clear rp;
        
        %Define displacement to target points
        z_x = new_shape(1:size_shape_vec);
        z_y = new_shape(size_shape_vec+1:2*size_shape_vec);
        
        %
        % Compute TPS warp from mean to new ring
        f_xr = tps_warp(tps_L_inv, s_x, s_y, z_x, i_xr, i_yr);
        f_yr = tps_warp(tps_L_inv, s_x, s_y, z_y, i_xr, i_yr);
        
        % Make new background texture
        %%%%%%%%%%
        new_shape_ROI = 128 + 32*randn(1)+ 4*randn(bw_r, bw_c);
        %new_masses(ii).ring_par = mean(new_shape_ROI(:));
        % Sample new ring text vector
        %%%%%
        new_ring_tex = interp2(double(new_shape_ROI), f_xr, f_yr);
        clear f_xr f_yr;
        B_ring = P_ring'*(new_ring_tex - mean_ring)';
        %
        % Sample new texture vector conditioned on ring and shape
        %
        %%%%
        mu = -C_mm_inv*C_nm*[B_c; B_ring];
        %mu = -C_mm_inv*C_nm*[B_shape; B_ring];
        %mu = -C_mm_inv*C_nm*B_ring
        [e_vec, e_val] = eig(C_mm_inv);
        
        B_tex = e_vec*sqrtm(e_val)*randn(k_tex, 1) + mu;
        %display(['b_tex = ', num2str(B_tex), '; b_ring = ', num2str(B_ring)]);
        
        new_tex = mean_tex + (P_tex*B_tex)';
        %new_masses(ii).tex_par = mean(new_tex);
        % Compute TPS warp to map from mean to new shape
        %%%%       
        f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
        f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
        

        % For each pixel in new_shape find nearest warped pixel
        display('find nearest');

        N = length(new_shape_pl(:,1));
        new_shape_tex = zeros(N, 1);
        for jj = 1:N
            pix_dist = (f_x - new_shape_pl(jj,1)).^2 + ...
                (f_y - new_shape_pl(jj,2)).^2;
            [min_val, nearest_pix] = min(pix_dist);
            new_shape_tex(jj) = new_tex(nearest_pix);
        end
        clear new_tex;

        
        new_shape_ROI(sub2ind([bw_r bw_c], new_shape_pl(:,2), new_shape_pl(:,1)))...
            = uint8(new_shape_tex);

        if if_plot
            figure(f1);
            subplot(2,2, rem(ii-1,4)+1);       
            imagesc(new_shape_ROI); axis image; hold on;
            %plot(new_shape(1:size_shape_vec), new_shape(size_shape_vec+1:2*size_shape_vec), 'y','LineWidth', .5);
            %plot(f_x, f_y, 'rx');
            figure(f2);
            subplot(2,2,rem(ii-1,4)+1);       
            image(new_shape_ROI); axis image; colormap(gray(256)); hold on;
            plot(new_shape(1:size_shape_vec), new_shape(size_shape_vec+1:2*size_shape_vec), 'r','LineWidth', .5);
            if ~rem(ii,4)
                f1 = figure;
                f2 = figure;
            end
        end
        %new_masses(ii).er = ((new_masses(ii).dif > 0) == (new_masses(ii).ring_par < new_masses(ii).tex_par))...
        %    | isnan(new_masses(ii).tex_par);
        new_masses(ii).mass_outline = new_shape;
        new_masses(ii).mass_ROI = new_shape_ROI;
        new_masses(ii).B_tex = B_tex;
        new_masses(ii).B_shape = B_shape;
        new_masses(ii).B_ring = B_ring;
        hold on;
        plot3(B_tex, B_ring, B_shape, 'bo');
        clear new_shape_tex new_shape_pl new_shape_ROI;
         ii = ii + 1;
    else
        bad_shape = bad_shape + 1;
    end
end
display(['Number of bad shapes = ', num2str(bad_shape)]);
%new_masses = 0;