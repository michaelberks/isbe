function [new_masses] = sample_new_masses(mass_model, no_of_masses, if_plot)

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
% B_shape     = mass_model.B_shape;
L_shape     = mass_model.L_shape;
mean_scale  = mass_model.mean_scale;
P_scale     = mass_model.P_scale;
% B_scale     = mass_model.B_scale;
% L_scale     = mass_model.L_scale;
mean_tex    = mass_model.mean_tex;
P_tex       = mass_model.P_tex;
% B_tex       = mass_model.B_tex;
L_tex       = mass_model.L_tex;
% mean_n      = mass_model.mean_n;
% P_n         = mass_model.P_n;
% B_n         = mass_model.B_n;
% L_n         = mass_model.L_n;

% mean_c      = mass_model.mean_c;
P_c         = mass_model.P_c;
% B_c         = mass_model.B_c;
L_c         = mass_model.L_c;

W_shape     = mass_model.W_shape;
W_tex       = mass_model.W_tex;
W_scale     = mass_model.W_scale;
% W_n         = mass_model.W_n;

mean_shape_pl = mass_model.mean_shape_pl;
mean_dilate = mass_model.mean_dilate;
size_shape_vec = length(mean_shape) / 2;

k_shape     = length(L_shape);
k_tex       = length(L_tex);

scale_factor = mass_model.scale_factor;

bad_shape = 0;
  
%Define source points for TPS - as row vectors
s_x = mean_dilate(1:size_shape_vec);
s_y = mean_dilate(size_shape_vec+1:end);

%Define points to be interpolated by TPS - as row vectors
i_x = mean_shape_pl(:,1)';
i_y = mean_shape_pl(:,2)';
        
tps_L_inv = tps_weights(s_x, s_y);

%
% Generate specified number of masses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

new_masses(no_of_masses).mass_outline = [];
new_masses(no_of_masses).mass_ROI = [];
ii = 1;

while ii <= no_of_masses

    % sample new combined appearance vectors - assume normal distribution
    % of modes
    %%%

    
    % Compute new shape vector, texture vector and scale
    new_c = (randn(length(L_c), 1) .* sqrt(L_c));
    
    Q_shape = P_c(1:k_shape,:); 
    Q_tex = P_c(k_shape+1:k_shape + k_tex,:);
    %Q_scale = P_c(end-1, :);
    %Q_n = P_c(end, :);
    Q_scale = P_c(end, :);
    
    new_shape = mean_shape / scale_factor + (P_shape*Q_shape*new_c)' / W_shape;
    new_tex = mean_tex + (P_tex*Q_tex*new_c)' / W_tex;
    new_scale = mean_scale + (P_scale*Q_scale*new_c)' / W_scale;
    %new_masses(ii).n_spicules = round(mean_n + (P_n*Q_n*new_c)' / W_n);
    %if new_masses(ii).n_spicules < 0, new_masses(ii).n_spicules = 0; end;
    
    new_shape = new_shape*new_scale;
    new_shape = reshape(new_shape, [], 2);
    
    off_r = 101 - round(min(new_shape(:,2)));
    off_c = 101 - round(min(new_shape(:,1)));
    
    new_shape(:,2) = new_shape(:,2) + off_r;
    new_shape(:,1) = new_shape(:,1) + off_c;
    
    bw_r = round(max(new_shape(:,2)))+ 100;
    bw_c = round(max(new_shape(:,1))) + 100;
        
    new_bw = roipoly(bw_r, bw_c, new_shape(:,1), new_shape(:,2));
    
    for jj = 1:49; new_bw = imdilate(new_bw, strel('disk', 1)); end
    [new_label, num_obj] = bwlabel(new_bw, 4);    

    if num_obj == 1
        rp = regionprops(new_label, 'PixelList'); clear new_label;
        new_shape_pl = rp.PixelList; clear rp;
        clear rp;
        
        % Compute TPS warp to map from mean to new shape
        %%%%
    
        [start_r start_c] = find(new_bw, 1);
        dilate_outline = bwtraceboundary(new_bw, [start_r, start_c], 'E');
        dilate_outline = [dilate_outline(:,2) dilate_outline(:,1)];
        clear new_bw;
        
        px = new_shape(1,1);
        py = new_shape(1,2);
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
        
        new_shape_tex = griddata(f_x, f_y, new_tex,...
            new_shape_pl(:,1), new_shape_pl(:,2));
        new_shape_tex(isnan(new_shape_tex)) = 0;
        new_shape_ROI = zeros(bw_r, bw_c);
        new_shape_ROI(sub2ind([bw_r bw_c], new_shape_pl(:,2), new_shape_pl(:,1)))...
            = new_shape_tex;
        
        display('completed grid');

        if if_plot
            
            %figure('WindowStyle', 'docked');
            %imagesc(new_shape_ROI); axis image; colormap(gray(256));
            %hold on
            %plot(new_shape(:,1),new_shape(:,2), 'y','LineWidth', .5);
            %plot(f_x(1:3:end), f_y(1:3:end), 'rx', 'MarkerSize', .5);
            figure('WindowStyle', 'docked');
            imagesc(uint8(new_shape_ROI)); axis image; colormap(gray(256));
        end

        new_masses(ii).mass_outline = new_shape;
        new_masses(ii).mass_ROI = new_shape_ROI;
        clear new_shape_tex new_shape_pl new_shape_ROI new_tex;
        ii = ii + 1;
    else
        bad_shape = bad_shape + 1;
    end
end
display(['Number of bad shapes = ', num2str(bad_shape)]);
%new_masses = 0;