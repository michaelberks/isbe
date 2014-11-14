% SAMPLE_NEW_SPICULES - sample new spicules from spicule appearance model 
%
% Usage: [new_spicules] = sample_new_spicules(spicule_AM, n_spicules,...
%                                               width, if_plot)
% 
% Arguments:
% spicule_AM - spicule appearance model
% n_spicules - scalar integer, number of new spicules to sample
% width - scalar integer, the half width of the new spicule normal profiles
% if_plot - [0 1] turn on/off displayinf new spciules
%
% Returns: 
% new_spicules - n_spicules-length structure of new spicules
%
% Notes:
% 
%
% See also: GENERATE_SPICULE_AM SAMPLE_NEW_MASSES
%
% References:
%
% Author:   Michael Berks
%           Imaging Science and Biomedical Engineering
%           University of Manchester
%
function [new_spicules] = sample_new_spicules(spicule_AM, n_spicules, width, if_plot);

    mean_s  = spicule_AM.mean_s;
    P_s     = spicule_AM.P_s;
    B_s     = spicule_AM.B_s;
    L_s     = spicule_AM.L_s;
    mean_w  = spicule_AM.mean_w;
    P_w     = spicule_AM.P_w;
    B_w     = spicule_AM.B_w;
    L_w     = spicule_AM.L_w;
    mean_b  = spicule_AM.mean_b;
    P_b     = spicule_AM.P_b;
    B_b     = spicule_AM.B_b;
    L_b     = spicule_AM.L_b;
    mean_p  = spicule_AM.mean_p;
    P_p     = spicule_AM.P_p;
    B_p     = spicule_AM.B_p;
    L_p     = spicule_AM.L_p;
    mean_l  = spicule_AM.mean_l;
    P_l     = spicule_AM.P_l;
    B_l     = spicule_AM.B_l;
    L_l     = spicule_AM.L_l;

    mean_c  = spicule_AM.mean_c;
    P_c     = spicule_AM.P_c;
    B_c     = spicule_AM.B_c;
    L_c     = spicule_AM.L_c;

    W_s     = spicule_AM.W_s;
    W_w     = spicule_AM.W_w;
    W_b     = spicule_AM.W_b;
    W_p     = spicule_AM.W_p;
    W_l     = spicule_AM.W_l;

    m = length(mean_s) / 2;

    k_s = length(L_s);
    k_w = length(L_w);
    k_b = length(L_b);
    k_p = length(L_p);
    
    bad_shape = 0;

%     if if_plot,
%         f1 = figure;
%         f2 = figure;
%     end


    %
    % Generate specified number of masses
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ii = 1;
    while ii <= n_spicules

        % sample new combined appearance vectors - assume normal distribution
        % of modes
        %%%


        % Compute new shape vector, texture vector and scale
        new_c = (randn(length(L_c), 1) .* sqrt(L_c));

        Q_s = P_c(1:k_s,:); 
        Q_w = P_c(k_s+1:k_s + k_w,:);
        Q_b = P_c(k_s+k_w+1:k_s+k_w+k_b,:);
        Q_p = P_p(k_s+k_w+k_b+1:end-1,:);
        Q_l = P_c(end, :);

        new_s = mean_s + (P_s*Q_s*new_c)' / W_s;
        new_w = mean_w + (P_w*Q_w*new_c)' / W_w;
        new_b = mean_b + (P_b*Q_b*new_c)' / W_b;
        new_p = mean_p + (P_p*Q_b*new_p)' / W_p;
        new_l = mean_l + (P_l*Q_l*new_c)' / W_l;
        
        new_b(new_b < 0) = 0;
        new_w(new_w < 1) = 1;

        new_s = new_s*new_l / m;
        new_shape(:,1) = new_s(1:m);
        new_shape(:,2) = new_s(m+1:2*m);
        clear new_s;

%         new_shape = spicule_AM.spicule;
%         new_b = spicule_AM.b;
%         new_w = spicule_AM.w;
        
        %calculate length of new spicule
        V = diff(new_shape);
        N = [0; cumsum(sqrt(V(:,1).^2 + V(:,2).^2))];
        
        %create landmark points spaced along pixel length for new shape and
        %profile
        n_pts = round(N(end) / 5);%
        landmarks = interp1(N, new_shape, linspace(0, N(end), n_pts));
        
        % Reconstruct profile appearance vectors
        profile_x = [1:N(end)+1];
        prows = 2*width + 1; pcols = length(profile_x);
        
        landmarks_xp = linspace(1, floor(N(end) + 1), n_pts);
        landmarks_yp = repmat(width+1, 1, n_pts);
        landmarks_w = ...
            interp1(linspace(1, pcols, length(new_w)), new_w, landmarks_xp);    
        upper_yp = landmarks_yp + landmarks_w;
        lower_yp = landmarks_yp - landmarks_w;
        
        %calculate widths at landmarks along profile from new_w
        w_vector = interp1(linspace(1, pcols, length(new_w)), new_w, 1:pcols);
        w_matrix = repmat(w_vector, prows, 1);
        
        %calculate brightness along profile from new_b
        b_vector = interp1(linspace(1, pcols, length(new_b)), new_b, 1:pcols);
        b_matrix = repmat(b_vector, prows, 1);

        offsets = repmat([-width:width]', 1, pcols);
        offsets(abs(offsets) >= w_matrix) = w_matrix(abs(offsets) >= w_matrix);
        
        new_profile = b_matrix .* sqrt(1 - (offsets.^2 ./ w_matrix.^2));
        
        %create upper and lower border landmarks for spicule shape
        [fx, fy] = gradient(landmarks);
        fy = fy ./ [sqrt(sum((fy.^2)')'), sqrt(sum((fy.^2)')')];
        
        landmarks_x = landmarks(:,1)'; 
        landmarks_y = landmarks(:,2)'; clear landmarks;
        
        
        upper_x = landmarks_x - landmarks_w.*fy(:,2)';
        upper_y = landmarks_y + landmarks_w.*fy(:,1)';
        lower_x = landmarks_x + landmarks_w.*fy(:,2)';
        lower_y = landmarks_y - landmarks_w.*fy(:,1)';
        
        %translate xy-coords to fit in box
        xbuff = 10 - min([upper_x, lower_x]);
        upper_x = upper_x + xbuff;
        lower_x = lower_x + xbuff;
        landmarks_x = landmarks_x + xbuff;
        new_shape(:,1) = new_shape(:,1) + xbuff;
        
        ybuff = 10 - min([upper_y, lower_y]); %NB: if spicule "hooks" then min could be in lower or upper
        upper_y = upper_y + ybuff;
        lower_y = lower_y + ybuff;
        landmarks_y = landmarks_y + ybuff;
        new_shape(:,2) = new_shape(:,2) + ybuff;
            
        %create box big enough for new spicule to fit into
        new_r = round(max([upper_y, lower_y])) + 10;
        new_c = round(max([upper_x, lower_x])) + 10;
        new_spicule_ROI = zeros(new_r, new_c);
            
        %create pixel list for new spicule, mask outline is upper border
        %landmarks from 1 to end, then lower landmarks from end to 1
        new_bw = roipoly(new_r, new_c, [upper_x fliplr(lower_x)],...
            [upper_y fliplr(lower_y)]);
        
        rp = regionprops(bwlabel(new_bw, 4), 'PixelList');
        %clear new_bw;
        
        %if new spicule shape is valid, warp intensity values from profile
        if size(rp, 1) == 1, 
            new_shape_pl = rp.PixelList; clear rp;
    
            % Compute TPS warp to map from mean to new shape
            %%%%
            %Define source points for TPS - as row vectors
            s_x = [upper_x landmarks_x lower_x];
            s_y = [upper_y landmarks_y lower_y];
    
            %Define points to be interpolated by TPS - as row vectors
            i_x = new_shape_pl(:,1)';
            i_y = new_shape_pl(:,2)';
    
            tps_L_inv = tps_weights(s_x, s_y);
    
            %Define displacement to target points
            z_x = [landmarks_xp landmarks_xp landmarks_xp];
            z_y = [upper_yp landmarks_yp lower_yp];
    
            %Compute displacement of interpolated points        
            f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y);
            f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y);
            
            new_tex = interp2(double(new_profile), f_x, f_y);            
            new_spicule_ROI(sub2ind([new_r, new_c], i_y, i_x)) =...
                uint8(new_tex);
    
            new_spicules(ii).new_shape = new_shape;
            new_spicules(ii).new_spicule_ROI = new_spicule_ROI;
            
            new_spicules(ii).s_orientation = newp(1);
            new_spicules(ii).s_location = newp(2);
            new_spicules(ii).s_distance = newp(3);
            
            if if_plot
                %subplot(4,5,ii);
                figure; 
                subplot(2,1,1); hold on;
                imagesc(new_profile); axis image; colormap gray;
                plot(landmarks_xp, upper_yp);
                plot(landmarks_xp, landmarks_yp);
                plot(landmarks_xp, lower_yp);
                plot(landmarks_xp, upper_yp, 'rx');
                plot(landmarks_xp, landmarks_yp, 'yx');
                plot(landmarks_xp, lower_yp, 'gx');
                %plot(f_x, f_y, 'c.');

                subplot(2,1,2); hold on;
                imagesc(new_spicule_ROI); axis image; colormap gray;
                plot(upper_x, upper_y);
                plot(landmarks_x, landmarks_y);
                plot(lower_x, lower_y);
                plot(upper_x, upper_y, 'rx');
                plot(landmarks_x, landmarks_y, 'yx');
                plot(lower_x, lower_y, 'gx');
                %plot(i_x, i_y, 'c.');
            end
            clear new_shape_pl new_shape new_spicule_ROI;
            ii = ii + 1;
        else
            bad_shape = bad_shape + 1;
        end
    end
    display(['Number of bad shapes = ', num2str(bad_shape)]);
end