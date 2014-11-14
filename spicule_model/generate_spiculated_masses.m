% generate_spiculated_masses - sample new spicules from spicule appearance model 
%
% Usage: [new_spicules] = generate_spiculated_masses(spicule_AM, mass_AM
%                                           n_masses, width, if_plot)
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
function [s_masses] = generate_spiculated_masses(spicule_AM,...
    mass_AM, n_masses, width, if_plot)
    
    f1 = figure;
    f2 = figure;
    for ii = 1:n_masses
        m = sample_new_masses(mass_AM, 1, 0);
        mass_outline = [m.mass_outline(1:end/2)'...
            m.mass_outline(end/2 +1:end)'];
        mass_ROI = m.mass_ROI;
        [r c] = size(mass_ROI);
        
        n_spicules = m.n_spicules, clear m;
        centroid = mean(mass_outline);
        
        jj = 1;
        while jj <= n_spicules
            s = sample_new_spicules1(spicule_AM, 1, 40, 0);
            s_shape = s.landmarks;
            s_location = s.s_location;
            s_distance = s.s_distance;
            s_orientation = s.s_orientation;
            landmarks_w = s.landmarks_w;
            profile = s.profile;
            
            m_lengths = [0; cumsum(sqrt(diff(mass_outline(:,1)).^2 ...
                + diff(mass_outline(:,2)).^2))];
            m_lengths = m_lengths / m_lengths(end);

            [mm idx] = min(abs(m_lengths - s_location));

            intersect = mass_outline(idx, :);
            c_to_m = intersect - centroid;

            theta = atan2(c_to_m(2), c_to_m(1)) + s_orientation;

            s_shape = ([cos(theta) -sin(theta); sin(theta) cos(theta)]*s_shape')'; 

            s_start = centroid + s_distance*c_to_m;
            s_shape = s_shape + repmat(s_start, size(s_shape, 1), 1);
            
            for kk = 1:jj-1
                if_cross(kk) = line_cross(s_shape, spicules(kk).s_shape);
            end
            if jj > 1, if sum(if_cross); continue; end; end;
            
            %try
                
                %create upper and lower border landmarks for spicule shape
            [fx, fy] = gradient(s_shape);
            fy = fy ./ [sqrt(sum((fy.^2)')'), sqrt(sum((fy.^2)')')];

            upper_x = s_shape(:,1)' - landmarks_w.*fy(:,2)';
            upper_y = s_shape(:,2)' + landmarks_w.*fy(:,1)';
            lower_x = s_shape(:,1)' + landmarks_w.*fy(:,2)';
            lower_y = s_shape(:,2)' - landmarks_w.*fy(:,1)';

            %translate ROI to fit spicule
            xbuff = floor(min(min([upper_x, lower_x]) - 10, 0));
            if xbuff
                upper_x = upper_x - xbuff;
                lower_x = lower_x - xbuff;
                s_shape(:,1) = s_shape(:,1) - xbuff;
                for kk = 1:jj - 1;
                    spicules(kk).s_shape(:,1) = spicules(kk).s_shape(:,1) - xbuff;
                end    
                centroid(1) = centroid(1) - xbuff;
                mass_outline(:,1) = mass_outline(:,1) - xbuff;
                mass_ROI = [zeros(r, -xbuff) mass_ROI];
                [r c] = size(mass_ROI);
            end

            xbuff = ceil(max(max([upper_x, lower_x]) - c + 10, 0));
            if xbuff
                mass_ROI = [mass_ROI zeros(r, xbuff)];
                [r c] = size(mass_ROI);
            end

            ybuff = floor(min(min([upper_y, lower_y]) - 10, 0));
            if ybuff
                upper_y = upper_y - ybuff;
                lower_y = lower_y - ybuff;
                s_shape(:,2) = s_shape(:,2) - ybuff;
                for kk = 1:jj - 1;
                    spicules(kk).s_shape(:,2) = spicules(kk).s_shape(:,2) - ybuff;
                end 
                centroid(2) = centroid(2) - ybuff;
                mass_outline(:,2) = mass_outline(:,2) - ybuff;
                mass_ROI = [zeros(-ybuff, c); mass_ROI];
                [r c] = size(mass_ROI);
            end

            ybuff = ceil(max(max([upper_y, lower_y]) - r + 10, 0));
            if ybuff
                mass_ROI = [mass_ROI; zeros(ybuff, c)];
                [r c] = size(mass_ROI);
            end

            %create pixel list for new spicule, mask outline is upper border
            %landmarks from 1 to end, then lower landmarks from end to 1
            new_bw = roipoly(r, c, [upper_x fliplr(lower_x)],...
                [upper_y fliplr(lower_y)]);

            rp = regionprops(bwlabel(new_bw, 4), 'PixelList');
            %clear new_bw;

            %if new spicule shape is valid, warp intensity values from profile
            if size(rp, 1) > 1, continue; end;

            new_shape_pl = rp.PixelList; clear rp;

            n_pts = size(s_shape, 1);
            [prows, pcols] = size(profile);
            landmarks_xp = linspace(1, pcols, n_pts);
            landmarks_yp = repmat((prows+1)/2, 1, n_pts);    
            upper_yp = landmarks_yp + landmarks_w;
            lower_yp = landmarks_yp - landmarks_w;

            % Compute TPS warp to map from mean to new shape
            %%%%
            %Define source points for TPS - as row vectors
            s_x = [upper_x s_shape(:,1)' lower_x];
            s_y = [upper_y s_shape(:,2)' lower_y];

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

            new_tex = interp2(double(profile), f_x, f_y);            
            mass_ROI(sub2ind([r, c], i_y, i_x)) =...
                mass_ROI(sub2ind([r, c], i_y, i_x)) + new_tex;

            spicules(jj).s_shape = s_shape;
            jj = jj + 1;
            %end
        end
        figure(f1); subplot(4,5,ii);
        imagesc(mass_ROI); colormap gray; axis image;
        
        figure(f2); subplot(4,5,ii); hold on;
        plot([mass_outline(:,1); mass_outline(1,1)],...
             [mass_outline(:,2); mass_outline(1,2)], 'LineWidth', 1.0);
        
        for jj = 1:n_spicules
            
            plot(spicules(jj).s_shape(:,1), spicules(jj).s_shape(:,2),...
                'LineWidth', 1.0);
        end
        axis equal ij;
        s_masses(ii).mass_outline = mass_outline;
        s_masses(ii).mass_ROI = mass_ROI;
        if n_spicules; s_masses(ii).spicules = spicules; clear spicules; end
    end
    %s_masses = 0;
end