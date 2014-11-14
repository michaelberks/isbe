% Work experimenting on the relationship between dual-tree coefficients of
% mammographic texture within subbands of the same levels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First lets look at just the coefficients

for lev = 1:3
    
    %Get a load of data from our texture regions
    [dt_data ilp_data icp_data] = mb_get_dual_tree_data('C:\isbe\dev\background\dual_tree\normal512', lev, 32, 1, 1);
    
    %Find out which subband is maximal
    [max_dt dt_idx] = max(dt_data, [], 2);
    [max_ilp ilp_idx] = max(ilp_data, [], 2);
    [max_icp icp_idx] = max(icp_data, [], 2);

    %Circular shift each row of data so maximal sub-band is in column 1
    rot_dt = dt_data;
    rot_ilp = ilp_data;
    rot_icp = icp_data;
    for band = 1:6
        rot_dt(dt_idx == band,:) = circshift(dt_data(dt_idx == band, :), [0 2 - band]);
        rot_ilp(ilp_idx == band,:) = circshift(ilp_data(ilp_idx == band, :), [0 2 - band]);
        rot_icp(icp_idx == band,:) = circshift(icp_data(icp_idx == band, :), [0 2 - band]);
    end
    
    %Now produce histograms of:
    
    % 1) absolute magnitudes
%     figure('name', 'DT magnitudes');
%     for band = 1:6
%         dt_lev_abs = abs(rot_dt(:,band));
%         subplot(2,3,band); hist(dt_lev_abs);
%         title(['Level ', num2str(lev), ' Sub-band ', num2str(band)]);
%     end
    
    % 2) percentage of magnitudes for dt, ilp and icp
    figure('name', 'DT percentage of maximum magnitudes');
    for band = 1:6
        dt_lev_abs_percent = abs(rot_dt(:,band)) ./ abs(rot_dt(:,2));
        subplot(2,3,band); hist(dt_lev_abs_percent, 100);
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','b','EdgeColor','b');
        title(['Level ', num2str(lev), ' Sub-band ', num2str(band)]);
    end
    saveas(gcf, ['C:\isbe\dev\background\dt_subbands\figures\dt_subband_mag_percent_level', num2str(lev), '.eps']);
    saveas(gcf, ['C:\isbe\dev\background\dt_subbands\figures\dt_subband_mag_percent_level', num2str(lev), '.bmp']);
    
    figure('name', 'ILP percentage of maximum magnitudes');
    for band = 1:6
        ilp_lev_abs_percent = abs(rot_ilp(:,band)) ./ abs(rot_ilp(:,2));
        subplot(2,3,band); hist(ilp_lev_abs_percent, 100);
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','b','EdgeColor','b');
        title(['Level ', num2str(lev), ' Sub-band ', num2str(band)]);
    end
    saveas(gcf, ['C:\isbe\dev\background\dt_subbands\figures\ilp_subband_mag_percent_level', num2str(lev), '.eps']);
    saveas(gcf, ['C:\isbe\dev\background\dt_subbands\figures\ilp_subband_mag_percent_level', num2str(lev), '.bmp']);

    figure('name', 'ICP percentage of maximum magnitudes');
    for band = 1:6
        icp_lev_abs_percent = abs(rot_icp(:,band)) ./ abs(rot_icp(:,2));
        subplot(2,3,band); hist(icp_lev_abs_percent, 100);
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','b','EdgeColor','b');
        title(['Level ', num2str(lev), ' Sub-band ', num2str(band)]);
    end
    saveas(gcf, ['C:\isbe\dev\background\dt_subbands\figures\icp_subband_mag_percent_level', num2str(lev), '.eps']);
    saveas(gcf, ['C:\isbe\dev\background\dt_subbands\figures\icp_subband_mag_percent_level', num2str(lev), '.bmp']);
    
    % 3) phases
%     figure('name', 'DT angles');
%     for band = 1:6
%         dt_lev_angle = angle(rot_dt(:,band));
%         subplot(2,3,band); hist(dt_lev_angle);
%         title(['Level ', num2str(lev), ' Sub-band ', num2str(band)]);
%     end
    
    % 4) difference of phases
    figure('name', 'DT difference angles from maximum magnitudes');
    for band = 1:6
        dt_lev_angle_diff = angle(rot_dt(:,band) .* conj(rot_dt(:,2)));
        subplot(2,3,band); hist(dt_lev_angle_diff, 100);
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','b','EdgeColor','b');
        title(['Level ', num2str(lev), ' Sub-band ', num2str(band)]);
    end  
    saveas(gcf, ['C:\isbe\dev\background\dt_subbands\figures\subband_phase_diff_level', num2str(lev), '.eps']);
    saveas(gcf, ['C:\isbe\dev\background\dt_subbands\figures\subband_phase_diff_level', num2str(lev), '.bmp']);
    
    figure('name', 'ILP difference angles from maximum magnitudes');
    for band = 1:6
        ilp_lev_angle_diff = angle(rot_ilp(:,band) .* conj(rot_ilp(:,2)));
        subplot(2,3,band); hist(ilp_lev_angle_diff, 100);
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','b','EdgeColor','b');
        title(['Level ', num2str(lev), ' Sub-band ', num2str(band)]);
    end  
    saveas(gcf, ['C:\isbe\dev\background\dt_subbands\figures\ilp_subband_phase_diff_level', num2str(lev), '.eps']);
    saveas(gcf, ['C:\isbe\dev\background\dt_subbands\figures\ilp_subband_phase_diff_level', num2str(lev), '.bmp']);
    
    figure('name', 'ICP difference angles from maximum magnitudes');
    for band = 1:6
        icp_lev_angle_diff = angle(rot_icp(:,band) .* conj(rot_icp(:,2)));
        subplot(2,3,band); hist(icp_lev_angle_diff, 100);
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','b','EdgeColor','b');
        title(['Level ', num2str(lev), ' Sub-band ', num2str(band)]);
    end  
    saveas(gcf, ['C:\isbe\dev\background\dt_subbands\figures\icp_subband_phase_diff_level', num2str(lev), '.eps']);
    saveas(gcf, ['C:\isbe\dev\background\dt_subbands\figures\icp_subband_phase_diff_level', num2str(lev), '.bmp']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
ilp_lev_angle_diff = angle(rot_ilp .* repmat(conj(rot_ilp(:,2)),1,6));
figure; plot3(ilp_lev_angle_diff(1:10000,3), ilp_lev_angle_diff(1:10000,4), ilp_lev_angle_diff(1:10000,5), 'r.');
figure; plot(ilp_lev_angle_diff(1:10000,3), ilp_lev_angle_diff(1:10000,1), 'r.');
figure; plot(ilp_lev_angle_diff(1:10000,4), ilp_lev_angle_diff(1:10000,6), 'r.');
%%
icp_lev_angle_diff = angle(rot_icp .* repmat(conj(rot_icp(:,2)),1,6));
figure; plot3(icp_lev_angle_diff(1:10000,3), icp_lev_angle_diff(1:10000,4), icp_lev_angle_diff(1:10000,5), 'r.');
figure; plot(icp_lev_angle_diff(1:10000,3), icp_lev_angle_diff(1:10000,1), 'r.');
figure; plot(icp_lev_angle_diff(1:10000,4), icp_lev_angle_diff(1:10000,6), 'r.');
%%
clear
angle_wrap = [...
    -pi/2 -pi -pi/4 -pi -pi -3*pi/4 ;...
    -3*pi/4 -pi -pi/4 -pi/4 -pi -3*pi/4 ; ...
    -3*pi/4 -pi -pi/2 -pi/4 -pi -pi ;
    -pi/2 -pi -pi/4 -pi -pi -3*pi/4 ;...
    -3*pi/4 -pi -pi/4 -pi/4 -pi -3*pi/4 ; ...
    -3*pi/4 -pi -pi/2 -pi/4 -pi -pi ;...
    ];

band_centres(:,:,1) = [0 pi; 0 0; pi/4 5*pi/4; -3*pi/4 pi/4; -pi/2 pi/2; -pi/4 3*pi/4];
band_centres(:,:,2) = [-pi/4 3*pi/4; 0 0; pi/4 5*pi/4; pi/2 3*pi/2; -pi/2 pi/2; -pi/4 3*pi/4];
band_centres(:,:,3) = [-pi/4 3*pi/4; 0 0; 0 pi; pi/4 5*pi/4; -pi/2 pi/2; -pi/4 3*pi/4];
band_centres(:,:,4) = [0 pi; 0 0; pi/4 5*pi/4; -3*pi/4 pi/4; -pi/2 pi/2; -pi/4 3*pi/4];
band_centres(:,:,5) = [-pi/4 3*pi/4; 0 0; pi/4 5*pi/4; pi/2 3*pi/2; -pi/2 pi/2; -pi/4 3*pi/4];
band_centres(:,:,6) = [-pi/4 3*pi/4; 0 0; 0 pi; pi/4 5*pi/4; -pi/2 pi/2; -pi/4 3*pi/4];

bin_centres = linspace(-2*pi, 2*pi, 201);
bin_width = bin_centres(2) - bin_centres(1);

%    
% Come again sunshine - but let's look at each band individually
for lev = 1:2
    
    %Get a load of data from our texture regions
    [dt_data ilp_data icp_data] = mb_get_dual_tree_data('C:\isbe\dev\background\dual_tree\normal512', lev, 32, 1, 1);
    
    %Find out which subband is maximal
    [max_dt dt_idx] = max(dt_data, [], 2);
    [max_ilp ilp_idx] = max(ilp_data, [], 2);
    [max_icp icp_idx] = max(icp_data, [], 2);
%
    %Circular shift each row of data so maximal sub-band is in column 1
    for sb = 1:6
        band_dt_data = circshift(dt_data(dt_idx == sb, :), [0 2 - sb]);
        band_ilp_data = circshift(ilp_data(ilp_idx == sb, :), [0 2 - sb]);
        band_icp_data = circshift(icp_data(icp_idx == sb, :), [0 2 - sb]);
        n_pts = size(band_icp_data, 1);
        
%         figure('name', ['DT percentage of maximum magnitudes - Band ', num2str(sb)] );
%         for band = 1:6
%             dt_lev_abs_percent = abs(band_dt_data(:,band)) ./ abs(band_dt_data(:,2));
%             subplot(2,3,band); hist(dt_lev_abs_percent, 100);
%             h = findobj(gca,'Type','patch');
%             set(h,'FaceColor','b','EdgeColor','b');
%             title(['Level ', num2str(lev), ' Sub-band ', num2str(band)]);
%         end
%         
%         figure('name', ['DT difference angles from maximum magnitudes - Band ', num2str(sb)]);
%         for band = 1:6
%             dt_lev_angle_diff = angle(band_dt_data(:,band) .* conj(band_dt_data(:,2)));
%             subplot(2,3,band); hist(dt_lev_angle_diff, 100);
%             h = findobj(gca,'Type','patch');
%             set(h,'FaceColor','b','EdgeColor','b');
%             title(['Level ', num2str(lev), ' Sub-band ', num2str(band)]);
%         end
        
%         figure('name', ['ILP difference angles from maximum magnitudes - Band ', num2str(sb)]);
%         for band = 1:6
%             ilp_lev_angle_diff = angle(band_ilp_data(:,band) .* conj(band_ilp_data(:,2)));
%             subplot(2,3,band); hist(ilp_lev_angle_diff, 100);
%             h = findobj(gca,'Type','patch');
%             set(h,'FaceColor','b','EdgeColor','b');
%             title(['Level ', num2str(lev), ' Sub-band ', num2str(band)]);
%         end
%         
        figure('name', ['ICP difference angles from maximum magnitudes - Band ', num2str(sb)]);
        for band = 1:6
            icp_lev_angle_diff = angle(band_icp_data(:,band) .* conj(band_icp_data(:,2)));
            swap_idx = icp_lev_angle_diff < angle_wrap(sb, band);
            icp_lev_angle_diff(swap_idx) = icp_lev_angle_diff(swap_idx) + pi;%2*pi;
            swap_idx = icp_lev_angle_diff > (angle_wrap(sb, band) + pi);
            icp_lev_angle_diff(swap_idx) = icp_lev_angle_diff(swap_idx) - pi;
            
            subplot(2,3,band); hist(icp_lev_angle_diff, 100);
            h = findobj(gca,'Type','patch');
            set(h,'FaceColor','b','EdgeColor','b');
            title(['Level ', num2str(lev), ' Sub-band ', num2str(band)]);
        end

%         figure('name', ['ICP difference angles from maximum magnitudes - Band ', num2str(sb)]);
%         for band = [1 3:6];
%             %icp_lev_angle_diff = angle(band_icp_data(:,band) .* conj(band_icp_data(:,2)));          
%             %swap_idx = icp_lev_angle_diff < angle_wrap(sb, band);
%             %icp_lev_angle_diff(swap_idx) = icp_lev_angle_diff(swap_idx) + 2*pi;
%              
%             icp_lev_angle_diff = band_icp_data(:,band) .* conj(band_icp_data(:,2));
%             icp_lev_angle_diff_v = [real(icp_lev_angle_diff) imag(icp_lev_angle_diff)];
%             
%             [cluster_idx, cc, SUMD, D] = ...
%                 kmeans(icp_lev_angle_diff_v, 2, 'Distance', 'cosine', 'start', [0 1; 0 -1]);
%                 %kmeans(icp_lev_angle_diff, 2, 'start', band_centres(band,:,sb)');            
%             
%             icp_lev_angle_diff = angle(icp_lev_angle_diff);
%             swap_idx = icp_lev_angle_diff < angle_wrap(sb, band);
%             icp_lev_angle_diff(swap_idx) = icp_lev_angle_diff(swap_idx) + 2*pi;
%             
%             cluster_centres = angle(complex(cc(:,1), cc(:,2)));
%             
%             c_means(1) = mean(icp_lev_angle_diff(cluster_idx == 1));
%             c_means(2) = mean(icp_lev_angle_diff(cluster_idx == 2));
%             c_var(1) = var(icp_lev_angle_diff(cluster_idx == 1));
%             c_var(2) = var(icp_lev_angle_diff(cluster_idx == 2));
%             
%             subplot(2,3,band); hold on;
%             
% %             [c1 x1] = hist(icp_lev_angle_diff(cluster_idx == 1), bin_centres); 
% %             h1 = bar(x1, c1/(bin_width*n_pts)); set(h1,'FaceColor','b','EdgeColor','b');
% %             plot([cluster_centres(1), cluster_centres(1)], [0 max(c1/(bin_width*n_pts))], 'y');
% %             [c2 x2] = hist(icp_lev_angle_diff(cluster_idx == 2), bin_centres); 
% %             h2 = bar(x2, c2/(bin_width*n_pts)); set(h2,'FaceColor','r','EdgeColor','r');
% %             plot([cluster_centres(2), cluster_centres(2)], [0 max(c2/(bin_width*n_pts))], 'y');
% %             
% %             plot(bin_centres, 0.5*exp(-(bin_centres - c_means(1)).^2 ./ (2*c_var(1)))./sqrt(2*pi*c_var(1)), 'c');
% %             plot(bin_centres, 0.5*exp(-(bin_centres - c_means(2)).^2 ./ (2*c_var(2)))./sqrt(2*pi*c_var(2)), 'c');
%             
%             if c_means(1) < c_means(2)
%                 icp_lev_angle_diff(cluster_idx == 1) = icp_lev_angle_diff(cluster_idx == 1) + pi;
%             else
%                 icp_lev_angle_diff(cluster_idx == 2) = icp_lev_angle_diff(cluster_idx == 2) + pi;
%             end
%             
%             [c1 x1] = hist(icp_lev_angle_diff(cluster_idx == 1), bin_centres);
%             [c2 x2] = hist(icp_lev_angle_diff(cluster_idx == 2), bin_centres); 
%             h = bar(x2, [c1' c2']/(bin_width*n_pts), 'stacked');
%             set(h(1),'FaceColor','b','EdgeColor','b');
%             set(h(2),'FaceColor','r','EdgeColor','r');
%             title(['Level ', num2str(lev), ' Sub-band ', num2str(band)]);
%         end

%         figure('name', ['ICP difference angles from maximum magnitudes - Band ', num2str(sb)]);
%         
%         icp_lev_angle_diff1 = angle(band_icp_data(1:1e4,5) .* conj(band_icp_data(1:1e4,2)));
%         swap_idx = icp_lev_angle_diff1 < angle_wrap(sb, 5);
%         icp_lev_angle_diff1(swap_idx) = icp_lev_angle_diff1(swap_idx) + 2*pi;
%         icp_lev_angle_diff3 = angle(band_icp_data(1:1e4,6) .* conj(band_icp_data(1:1e4,2)));
%         swap_idx = icp_lev_angle_diff3 < angle_wrap(sb, 6);
%         icp_lev_angle_diff3(swap_idx) = icp_lev_angle_diff3(swap_idx) + 2*pi;
%         icp_lev_angle_diff4 = angle(band_icp_data(1:1e4,1) .* conj(band_icp_data(1:1e4,2)));
%         swap_idx = icp_lev_angle_diff4 < angle_wrap(sb, 1);
%         icp_lev_angle_diff4(swap_idx) = icp_lev_angle_diff4(swap_idx) + 2*pi;
%         
%         plot3(icp_lev_angle_diff1, icp_lev_angle_diff3, icp_lev_angle_diff4, 'rx');
    end
%

end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I wanna try something else...
% For each subband in turn, look at the angle difference compared to its
% neighbouring band relative to its magnitude
%%

[dt_data] = mb_get_dual_tree_data('C:\isbe\dev\background\dual_tree\normal512', 2, 64);

% A1 magnitude of 1 against angle difference between 1 and 2
angle12_abs1 = zeros(size(dt_data, 1), 2);
angle12_abs1(:,1) = angle(dt_data(:,1) .* conj(dt_data(:,2)));
angle12_abs1(:,2) = abs(dt_data(:,1));
%
max_abs = max(angle12_abs1(:,2));
%
figure;

idx = angle12_abs1(:,2) < max_abs/8;
subplot(2,2,1); hist(angle12_abs1(idx, 1), 20);

idx = angle12_abs1(:,2) >= max_abs/8 & angle12_abs1(:,2) < 2*max_abs/8;
subplot(2,2,2); hist(angle12_abs1(idx, 1), 20);

idx = angle12_abs1(:,2) >= 2*max_abs/8 & angle12_abs1(:,2) < 3*max_abs/8;
subplot(2,2,3); hist(angle12_abs1(idx, 1), 20);

idx = angle12_abs1(:,2) >= 3*max_abs/8 & angle12_abs1(:,2);
subplot(2,2,4); hist(angle12_abs1(idx, 1), 10);