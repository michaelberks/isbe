%--------------------------------------------------------------------------
% Use function mb_dt_orientation sum to measure orientation across scales
%--------------------------------------------------------------------------

bg_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');

%load in a background region - lets measure the orientations for different
%rotations of the region
bg_0 = double(imread(['C:\isbe\dev\background\images\normal1024\', bg_list(4).name]));
%%
%Compute the mean orientation at each scale
[orientation_0, orientation_levels_0] = mb_dt_orientation_mean('Image', bg_0);

%%
colors = 'rgbymc';
figure; hold on; axis equal
for lev = 1:5
    plot([0 orientation_levels_0(lev)], colors(lev), 'Linewidth', 2, 'DisplayName', ['Level ', num2str(lev)]);
end
legend('show');
%%
figure;
for a = 1:6
    %Extend the bg, rotate, then take take the central region
    [r c] = size(bg_0);
    pad_vec = ceil(0.5*(sqrt(r^2 + c^2) - [r c]));
    
    bg_0_pad = padarray(bg_0, pad_vec, 'symmetric');
    bg_angle_pad = imrotate(bg_0_pad, angles(a), 'nearest', 'crop');
    bg_angle = bg_angle_pad(pad_vec(1)+(1:r), pad_vec(2)+(1:c));
    
    [o, o_levels] = mb_dt_orientation_mean('Image', bg_angle, 'Levels', 6);
    
    for lev = 1:6      
        subplot(2,3,lev);
        hold on; axis equal
        plot([o_levels(lev,1)/abs(o_levels(lev,1)) 0 o_levels(lev,2)/abs(o_levels(lev,2))],...
            colors(a), 'Linewidth', 2, 'DisplayName', [num2str(angles(a)), ' degrees']);
        if a == 6
            legend('show');
        end
    end
    
    %imagesc(bg_angle); axis image; colormap(gray(256));
end
    
%%
% Plot the rose diagram of ICP angle for each level of the bg
dt_bg_0 = dtwavexfm2(bg_0, 7);
[ilp icp] = mb_dual_tree_transform(dt_bg_0);
%
figure('Name', 'ICP orientation histograms at each scale - original image');
for lev = 1:6
    max_icp = max(icp{lev}, [], 3);
    subplot(2,3,lev); title(['Level ', num2str(lev)]);
    rose(angle(max_icp(:)), 100);
end
%%
%Now do this for each rotated copy of the image
for lev = 1:6
    f(lev) = figure('Name', ['ICP orientation histograms at level - ', num2str(lev)]); %#ok
end
%
angles = 15:30:165;
for a = 1:6
    %Extend the bg, rotate, then take take the central region
    [r c] = size(bg_0);
    pad_vec = ceil(0.5*(sqrt(r^2 + c^2) - [r c]));
    
    bg_0_pad = padarray(bg_0, pad_vec, 'symmetric');
    bg_angle_pad = imrotate(bg_0_pad, angles(a), 'nearest', 'crop');
    bg_angle = bg_angle_pad(pad_vec(1)+(1:r), pad_vec(2)+(1:c));
    
    dt_bg_angle = dtwavexfm2(bg_angle, 7);
    [ilp icp] = mb_dual_tree_transform(dt_bg_angle);
    %
    for lev = 1:6
        max_icp = max(icp{lev}, [], 3);
        %idx = imag(max_icp) < 0;
        %max_icp(idx) = -max_icp(idx); 
        
        figure(f(lev));
        subplot(2,3,a);
        rose(angle(max_icp(:)), 100);
        title([num2str(angles(a)), ' degrees']);
    end
    
    %imagesc(bg_angle); axis image; colormap(gray(256));
end
%%
%--------------------------------------------------------------------------
%A quick test on max icp angles...
% Make an image of concentric circles
for width = 1:10
    circles = ones(256);
    [xx yy] = find(circles);
    for rads = 0:(2*width):128
        idx = ((xx-128).^2 + (yy-128).^2 > rads^2) & ((xx-128).^2 + (yy-128).^2 <= (rads+width)^2);
        circles(idx) = 0;
    end
    figure; imagesc(circles); axis image;
    %
    dt_circles = dtwavexfm2(circles, 7);
    [ilp icp] = mb_dual_tree_transform(dt_circles);
    %%
    figure('Name', ['ICP orientation histograms at each scale - circles image, width =', num2str(width)]);
    for lev = 1:6
        max_icp = max(icp{lev}, [], 3);

        centre = 2^(7-lev);
        [xx yy] = find(ones(size(max_icp)));

        idx = (xx-centre).^2 + (yy-centre).^2 < centre^2;

        subplot(2,3,lev);
        rose(angle(max_icp(idx)), 100);
        title(['Level ', num2str(lev)]);
    end
end
%%
figure('Name', 'ICP orientations at each scale - circles image');
for lev = 1:6
    max_icp = max(icp{lev}, [], 3);
    
    centre = 2^(7-lev);
    [xx yy] = find(ones(size(max_icp)));
    
    idx = (xx-centre).^2 + (yy-centre).^2 > centre^2;
    
    max_icp(idx) = 0;
    subplot(2,3,lev); title(['Level ', num2str(lev)]);
    imagesc(2*angle(max_icp)); axis image; colormap([hsv(128); hsv(128)]);

end
%%
%--------------------------------------------------------------------------
%This really needs more work: are all orientation calculation biased
%towards +/-15,75 degree bands (and against the +/-45)

%Load in a load of data from each level across all our mammos and rose plot
%the orientation histograms. Across the data we would not expect any
%particular orientation to be more prevalent than any other, so should see
%approximately uniform histograms

%May as well also look at dt and ilp phase
f_dt = figure('Name', 'Max dual-tree phase across all mammo data');
f_ilp = figure('Name', 'Max ILP phase across all mammo data');
f_icp = figure('Name', 'Max ICP phase across all mammo data');
f_icp2 = figure('Name', 'Max ICP phase across all mammo data');
for lev = 2:5
    
    %Get a load of data from our texture regions
    [dt_data ilp_data icp_data] = mb_get_dual_tree_data('C:\isbe\dev\background\dual_tree\normal1024', lev, 32, 1, 1);
    
    %Find out which subband is maximal
    [max_dt dt_idx{lev}] = max(dt_data, [], 2); %#ok
    [max_ilp ilp_idx] = max(ilp_data, [], 2);
    [max_icp icp_idx] = max(icp_data, [], 2);
    
    figure(f_dt); subplot(2,2,lev-1);
    rose(angle(max_dt(:)), 100);
    title(['Level ', num2str(lev)]);
    
    figure(f_ilp); subplot(2,2,lev-1);
    rose(angle(max_ilp(:)), 100);
    title(['Level ', num2str(lev)]);
    
    figure(f_icp); subplot(2,2,lev-1);
    rose(angle(max_icp(:)), 100);
    title(['Level ', num2str(lev)]);
    
    %Produce magnitude weighted ICP bins
    figure(f_icp2); subplot(2,2,lev-1);
    weighted_complex_rose(max_icp, 100);
    title(['Level ', num2str(lev)]);
end
%
saveas(f_dt, 'C:\isbe\dev\background\orientations\max_dt_phase.fig');
saveas(f_ilp, 'C:\isbe\dev\background\orientations\max_ilp_phase.fig');
saveas(f_icp, 'C:\isbe\dev\background\orientations\max_icp_phase.fig');
saveas(f_icp2, 'C:\isbe\dev\background\orientations\max_icp_phase2.fig');
%%
% Does it make a difference if we use the ILP magnitudes instead? No - it
% is the underlying DT magnitudes causing the bias
icp_phase_ilp_mag = zeros(size(icp_data,1), 1);
for band = 1:6
    idx = dt_idx == band;
    icp_phase_ilp_mag(idx) = icp_data(idx, band);
end
figure;
rose(angle(icp_phase_ilp_mag), 100);
title(['Level ', num2str(lev)]);

%%
% So what is the bias, what percentage of pixels are maximal within in each
% subband across all the mammo data?
for lev = 2:5
    for band = 1:6
        display(['Band ', num2str(band), ': ', num2str(100*(sum(dt_idx{lev}==band)) / length(dt_idx{lev})), '% of points']);
    end
    display('Next level');
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Let's try and compute some orientations with gabor filters
% instead:
[response orientation] = gabor_filter('ImageIn', bg_0);
%%
figure; imagesc(response); axis image; colormap(gray(256));
figure; imagesc(bg_0); axis image; colormap(gray(256));
hold on;
quiver(1:4:1024, 1:4:1024,...
    response(1:4:end, 1:4:end).*cos(orientation(1:4:end, 1:4:end)),...
    -response(1:4:end, 1:4:end).*sin(orientation(1:4:end, 1:4:end)));
%%
[responses orientations] = mb_get_gabor_data('C:\isbe\dev\background\images\normal512', 32);

%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Let's try and compute some orientations with the monogenic signal
% instead:
[dummy, dummy, dummy, A, theta, psi] = monofilt(bg_0, 3, 4, 2, 0.65, 0);
clear dummy;
pack