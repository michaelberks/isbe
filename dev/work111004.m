[line_strength, line_orientation] = ...
    gaussian_clover_line(nailfold, [2 4]);%gaussian_2nd_derivative_line2

%Apply non-maximal suppression to skeletonise the line strength map
line_nms = mb_non_maximal_supp(line_strength, line_orientation);

%Apply 
[line_mask] = bwareaopen(line_nms > 0, 80);
[y_vessels x_vessels] = find(line_mask);
%%
nailfold = nailfold(:,:,1);
nailfold_mask = imerode(nailfold > 10 & nailfold < 250, strel('disk', 20));
n_min = min(nailfold(nailfold_mask));
n_max = max(nailfold(nailfold_mask));
%%
figure; imagesc(nailfold); axis image; colormap(gray(256)); hold on; caxis([n_min n_max]);
plot(maxima_pos12(:,1), maxima_pos12(:,2), 'r+', 'markersize', 2);
%%
num_pts = length(maxima_vals12);
for ii = 4%:num_pts
    x = maxima_pos12(ii,1);
    y = maxima_pos12(ii,2);
    neighbours = (abs(x_vessels - x) < 10) & (abs(y_vessels - y) < 10);
            
    if ~any(neighbours); continue; end
    
    xn = x_vessels(neighbours);
    yn = y_vessels(neighbours);

    dists = (xn - x).^2 + (yn - y).^2;
    [dummy min_idx] = min(dists);
    x = xn(min_idx);
    y = yn(min_idx);
    
    try_again = false;
    try
        vessel_east = bwtraceboundary(line_mask, [y x], 'E', 8, inf, 'clockwise');
    catch
        try_again = 1;
    end
    if try_again
        try
            vessel_east = bwtraceboundary(line_mask, [y x], 'NE', 8, inf, 'clockwise');
            try_again = false;
        catch
            try_again = 1;
        end
    end
    if try_again
        try
            vessel_east = bwtraceboundary(line_mask, [y x], 'SE', 8, inf, 'clockwise');
        catch
            vessel_east = [x y];
            display('total fail on east');
        end
    end
%     try_again = false;
%     try
%         vessel_west = bwtraceboundary(line_mask, [y x], 'W', 8, inf, 'counterclockwise');
%     catch
%         try_again = 1;
%     end
%     if try_again
%         try
%             vessel_west = bwtraceboundary(line_mask, [y x], 'NW', 8, inf, 'counterclockwise');
%             try_again = false;
%         catch
%             try_again = 1;
%         end
%     end
%     if try_again
%         try
%             vessel_west = bwtraceboundary(line_mask, [y x], 'SW', 8, inf, 'counterclockwise');
%         catch
%             vessel_west = [x y];
%             display('total fail on west');
%         end
%     end    
    
    
    plot(vessel_east(:,2), vessel_east(:,1), 'm-');
%     plot(vessel_west(:,2), vessel_west(:,1), 'g-');
    
end   

num_v = size(vessel_east,1);
keep = true(num_v,1);
for jj = 2:num_v
    keep(jj) = ~ismember(vessel_east(jj,:), vessel_east(1:jj-1,:), 'rows');
end
%%
figure; hold on;
plot(vessel_east(:,2), vessel_east(:,1), 'k.');  axis equal ij;

line_pts_l = track_line(line_mask, [y x], 'dir_order', {'W', 'S', 'E', 'N', 'SW', 'SE', 'NE', 'NW'});
line_pts_r = track_line(line_mask, [y x], 'dir_order', {'E', 'S', 'W', 'N', 'SE', 'SW', 'NW', 'NE'});
plot(line_pts_l(:,2), line_pts_l(:,1), 'r-');
plot(line_pts_r(:,2), line_pts_r(:,1), 'g-');
%%
% figure; imagesc(nailfold); axis image; colormap(gray(256)); hold on; caxis([n_min n_max]);
figure; imagesc(line_mask); axis image; colormap(gray(256)); hold on;
plot(maxima_pos12(:,1), maxima_pos12(:,2), 'r+', 'markersize', 2);
num_pts = length(maxima_vals12);
for ii = 1:num_pts
    x = maxima_pos12(ii,1);
    y = maxima_pos12(ii,2);
    neighbours = (abs(x_vessels - x) < 10) & (abs(y_vessels - y) < 10);
            
    if ~any(neighbours); continue; end
    
    xn = x_vessels(neighbours);
    yn = y_vessels(neighbours);

    dists = (xn - x).^2 + (yn - y).^2;
    [dummy min_idx] = min(dists);
    x = xn(min_idx);
    y = yn(min_idx);
    
    line_pts = track_line(line_mask, [y x], 'dir_order', {'W', 'S', 'E', 'N', 'SW', 'SE', 'NE', 'NW'});
    line_pts = track_line(line_mask, flipud(line_pts), 'dir_order', {'E', 'S', 'W', 'N', 'SE', 'SW', 'NW', 'NE'});
    plot(line_pts(:,2), line_pts(:,1), 'm-');
    
end   

