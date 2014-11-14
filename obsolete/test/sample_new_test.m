function [new_masses] = sample_new_test(mass_model, no_of_masses, if_plot)

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
% mean_n      = mass_model.mean_n;
% P_n         = mass_model.P_n;
% L_n         = mass_model.L_n;

%mean_c      = mass_model.mean_c;
P_c         = mass_model.P_c;
L_c         = mass_model.L_c;

W_shape     = mass_model.W_shape;
W_tex       = mass_model.W_tex;
W_scale     = mass_model.W_scale;
%W_n         = mass_model.W_n;

mean_shape_pl = mass_model.mean_shape_pl;

size_shape_vec = length(mean_shape) / 2;

k_shape     = length(L_shape);
k_tex       = length(L_tex);

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
        
tps_L_inv = tps_weights(s_x, s_y);

Q_shape = P_c(1:k_shape,:); 
Q_tex = P_c(k_shape+1:k_shape + k_tex,:);
Q_scale = P_c(end, :);
%Q_n = P_c(end, :);

%
% Generate specified number of masses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_masses(no_of_masses) = struct('mass_outline', [], 'mass_ROI', []);
ii = 1;
while ii <= no_of_masses

    % sample new combined appearance vectors - assume normal distribution
    % of modes
    %%%

    
    % Compute new shape vector, texture vector and scale
    new_c = (randn(length(L_c), 1) .* sqrt(L_c));
    
    new_shape = mean_shape / scale_factor + (P_shape*Q_shape*new_c)' / W_shape;
    new_tex = mean_tex + (P_tex*Q_tex*new_c)' / W_tex;
    new_scale = mean_scale + (P_scale*Q_scale*new_c)' / W_scale;
    %new_masses(ii).n_spicules = round(mean_n + (P_n*Q_n*new_c)' / W_n);
    %if new_masses(ii).n_spicules < 0, new_masses(ii).n_spicules = 0; end;
    
    new_shape = new_shape*new_scale;
    
    off_r = 101 - round(min(new_shape(size_shape_vec+1:2*size_shape_vec)));
    off_c = 101 - round(min(new_shape(1:size_shape_vec)));
    
    new_shape(size_shape_vec+1:2*size_shape_vec) = ...
        new_shape(size_shape_vec+1:2*size_shape_vec) + off_r;
    new_shape(1:size_shape_vec) = ...
        new_shape(1:size_shape_vec) + off_c;
    
    bw_r = round(max(new_shape(size_shape_vec+1:2*size_shape_vec)))+ 100;
    bw_c = round(max(new_shape(1:size_shape_vec))) + 100;
        
    new_bw = roipoly(bw_r, bw_c, new_shape(1:size_shape_vec),...
        new_shape(size_shape_vec+1:2*size_shape_vec));
        
    rp = regionprops(bwlabel(new_bw, 4), 'PixelList');
    clear new_bw;
    rp_l = size(rp, 1);
    
    if rp_l == 1, 
        new_shape_pl = rp.PixelList; clear rp;

        % Compute TPS warp to map from mean to new shape
        %%%%

        %Define displacement to target points
        z_x = new_shape(1:size_shape_vec);% - mean_shape(1:size_shape_vec) + mean_off_c;
        z_y = new_shape(size_shape_vec+1:2*size_shape_vec);%...
            %- mean_shape(size_shape_vec+1:2*size_shape_vec);

        %Compute displacement of interpolated points        
        f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
        f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
        
        display('completed warping');
        % Create new shape pixel list
        %%%%%%%%%%

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

        new_shape_ROI = reshape(uint8(255*rand(1) + 4*rand(bw_r*bw_c, 1)), bw_r, bw_c);
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

        new_masses(ii).mass_outline = new_shape;
        new_masses(ii).mass_ROI = new_shape_ROI;
        clear new_shape_tex new_shape_pl new_shape_ROI;
        ii = ii + 1;
    else
        bad_shape = bad_shape + 1;
    end
end
display(['Number of bad shapes = ', num2str(bad_shape)]);
%new_masses = 0;