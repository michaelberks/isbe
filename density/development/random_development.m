load goatfiend

%Display step wedge
figure; imagesc(sw_image); axis image; colormap(gray(256));

%Define some points for 5 steps
u(1,:) = [370, 63];
for ii = 2:6
    u(ii,:) = u(1,:) + (ii-1)*[12 0]; 
end
l = u;
l(:,2) = l(:,2) + 40;

%Create a tranformable poly
hp = impoly(gca, [u(1,:); l(1,:); l(2,:); u(2,:); u(3,:); l(3,:); l(4,:); u(4,:); u(5,:); l(5,:); l(6,:); u(6,:)], 'closed', 0);
api = iptgetapi(hp);

%Nudge horizontally (we can just drag though)
pos = api.getPosition(); pos(:,1) = pos(:,1) + 5; api.setPosition(pos);

%Rotate
theta = 5*pi/180; rot_mat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
pos = api.getPosition(); pos = (pos - repmat(mean(pos), size(pos,1), 1))*rot_mat + repmat(mean(pos), size(pos,1), 1); api.setPosition(pos);

%Stetch horizontally
pos = api.getPosition(); pos(:,1) = (pos(:,1) - pos(1,1))*1.1 + pos(1,1); api.setPosition(pos);
%%
%Go back again to the step wedge - can we find the edges of it?
canny_sw = edge(sw_image,'canny');

figure; imagesc(canny_sw); axis image; colormap(gray(256));
%%
[ll, no_objects] = bwlabel(canny_sw, 8); %#ok
[nn] = hist(ll(:), no_objects + 1);
[mm ind] = max(nn(2:end));
canny_mask = ll == ind;

figure; imagesc(canny_mask); axis image;
%%
[H,T,R] = hough(canny_sw);

H(:,[1:80 110:180]) = 0;
figure; imagesc(H, 'XData',T,'YData',R);
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;

P  = houghpeaks(H,10,'threshold',ceil(0.3*max(H(:))), 'NHoodSize', [11 11]);
x = T(P(:,2)); y = R(P(:,1));
plot(x,y,'s','color','white');
%
% Find lines and plot them
lines = houghlines(canny_sw,T,R,P,'FillGap',5,'MinLength',35);
figure; imagesc(sw_image); hold on; axis image; colormap(gray(256));
max_len = 0;
mid_xy = zeros(length(lines), 2);
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   mid_xy(k,:) = mean(xy);
   
   % Plot beginnings and ends of lines
   plot(mid_xy(k,1), mid_xy(k,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end

% highlight the longest line segment
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','cyan');

[p,S] = polyfit(mid_xy(:,1),mid_xy(:,2),1);
sx = 1; sy = sum(p);
ex = size(sw_image,2); ey = p(1)*ex + p(2);

plot([sx ex], [sy ey], 'r:');
%%
%--------------------------------------------------------------------------
% Convert the mammograms stored on C: to jpegs to make simple viewing
% easier
mam_list = dir('C:\isbe\density\mammograms\*.tif');
for ii = 1:length(mam_list)
    mam = imread(['C:\isbe\density\mammograms\', mam_list(ii).name]);
    mam = imresize(mam, [1024 NaN], 'bilinear');
    mam = uint8(mam / 2^8);
    if ~isempty(strfind(mam_list(ii).name, '1824'))
        mam = rot90(mam);
    end
    imwrite(mam, ['C:\isbe\density\mammograms\jpeg\', mam_list(ii).name(1:end-3), 'jpg']);
end
%%
mam_list = dir('C:\isbe\density\mammograms\*.tif');
for ii = 1:40
    mam = imread(['C:\isbe\density\mammograms\', mam_list(ii).name]);
    if ~isempty(strfind(mam_list(ii).name, '1824'))
        filmsizes = 1;
        max_pairs = 3;
    else
        filmsizes = 2;
        max_pairs = 4;
    end

    if filmsizes == 1
        % this is only needed for the 1824 film sizes
        mam = rot90(mam);
    end
    
    mam = medfilt2(mam);
    
    mam = imresize(mam,0.176);
    mam = mam / 16;
    mam = 4095 - mam;
    
    sw_image = rot90(mam(200:end-50, 75:175));    
    clear mam;
    
    canny_sw = edge(sw_image,'canny', 0.05, 2);
    figure; 
    subplot(2,1,1); imagesc(sw_image); axis image; colormap(gray(256)); hold on;
    subplot(2,1,2); imagesc(canny_sw); axis image; colormap(gray(256)); hold on;
%     [H,T,R] = hough(canny_sw, 'ThetaResolution', 2, 'RhoResolution', 2);
% 
%     H(:,[1:40 51:90]) = 0;
%     P  = houghpeaks(H,10);
%     lines = houghlines(canny_sw,T,R,P,'FillGap',5,'MinLength',10);
%     if length(lines) > 1
%         max_len = 0;
%         for k = 1:length(lines)
%            xy = [lines(k).point1; lines(k).point2];
%            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
% 
%            mid_xy(k,:) = mean(xy);
% 
%            % Plot beginnings and ends of lines
%            plot(mid_xy(k,1), mid_xy(k,2),'x','LineWidth',2,'Color','red');
% 
%            % Determine the endpoints of the longest line segment
%            len = norm(lines(k).point1 - lines(k).point2);
%            if ( len > max_len)
%               max_len = len;
%               xy_long = xy;
%            end
%         end
%         % highlight the longest line segment
%         plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','cyan');
% 
%         [p,S] = polyfit(mid_xy(:,1),mid_xy(:,2),1);
%         sx = 1; sy = sum(p);
%         ex = size(sw_image,2); ey = p(1)*ex + p(2);
% 
%         plot([sx ex], [sy ey], 'r:');
%     end
    
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
%Some work on Hough circles

%load an image
mam_list = dir('C:\isbe\density\mammograms\*.tif');
for ii = 1:10
    mam = rot90(imread(['C:\isbe\density\mammograms\', mam_list(ii).name]));
    figure; imagesc(mam); axis image; colormap(gray(256));
end

mam_region = mam([1:500 end-499:end],:);
%display it and impixel it
figure; imagesc(mam_region); axis image; colormap(gray(256));
[x_pts y_pts p] = impixel;

%
tic
for ii = 1:length(x_pts);
    
    [xc, yc, rad] = find_marker(mam_region, x_pts(ii), y_pts(ii), 0, gcf);
end
toc
%%
for ii = 1:length(x_pts);
    %for each point selected
    xi = x_pts(ii); yi = y_pts(ii);
    [r c] = size(mam_region);

    %specify size of box in which to search for marker
    box_size = 100;

    % get boundaries for box, ensuring we don't sample outside the image      
    x_start = max(1, xi-box_size);
    x_end = min(c, xi+box_size);
    y_start = max(1, yi-box_size);
    y_end = min(r, yi+box_size);

    %adjust the co-ordinates of the original point
    xii = xi - x_start + 1;
    yii = yi - y_start + 1;

    %Extract small box
    small = mam_region(y_start:y_end, x_start:x_end);
    small_med = medfilt2(small, [7 7]);
    [rs cs] = size(small);
    
    %Work out valid centres
    xx = repmat(1:cs, rs, 1);
    yy = repmat((1:rs)', 1, cs);
    cc = (xx-xii).^2 + (yy-yii).^2 < 1600; %max_r = 40^2
    [centres_y centres_x] = find(cc);
    
    %Do an initial canny edge detector to work out strong straight lines so
    %we can discard these
    cannysmall = canny_edge(small_med,0.5,2,[],0.5);
    
    %Compute the initial ignore_map, discarding outside the large central
    %circle
    ignore_map = (xx-xii).^2 + (yy-yii).^2 > 6400;
    
    %separate edges into 8-connected structures
    [edge_labels n] = bwlabel(cannysmall,8);
    
    %for each edge, add straight line segments to the ignore map
    for kk = 1:n
        these_edges = edge_labels == kk;
        
        %if it's of reasonable size
        if sum(these_edges(:)) > 20
            
            %work out the gradient
            [ey ex] = find(these_edges);
            [dummy fxy] = gradient([ex ey]);
            std_fxy = std(fxy);
            
            %if this is fairly constant (i.e. has a low SD) then we
            %probably have a straight line
            if sum(std_fxy) < 1;
                ignore_map(these_edges) = 1;
            end
        end
    end

    %Now repeat the cann edge detector using the new ignore map
    circular_edges = canny_edge(small_med,[],2,0.995,0.5, ignore_map, 1);  
    
    figure; 
    subplot(1,2,1); imagesc(small); axis image; colormap(gray(256));
    subplot(1,2,2); imagesc(circular_edges); axis image; colormap(gray(256));
    
    if any(circular_edges(:))
        %
        rad_scores = hough_circle1(circular_edges, [centres_x centres_y], 20, 75);
        rad_scores_pos = sort(rad_scores(rad_scores > 0), 'descend');
        [pos_circ_cent pos_circ_rad] = find(rad_scores >= rad_scores_pos(3));
        xc = zeros(length(pos_circ_cent), 1);
        yc = zeros(length(pos_circ_cent), 1);
        rad = zeros(length(pos_circ_cent), 1);
        for jj = 1:length(pos_circ_cent)
            xc(jj) = centres_x(pos_circ_cent(jj));
            yc(jj) = centres_y(pos_circ_cent(jj));
            rad(jj) = pos_circ_rad(jj) + 19;
            drawcircle(xc(jj),yc(jj),rad(jj),'g');
        end

        xc = mean(xc);
        yc = mean(yc);
        rad = mean(rad);
        drawcircle(xc,yc,rad,'r');
        subplot(1,2,1);
        drawcircle(xc,yc,rad,'r');
    end
end
%%
[max_rows cent_idx] = max(rad_scores);
[dummy max_circ_rad] = max(max_rows);
max_circ_cent = cent_idx(max_circ_rad);
rad = max_circ_rad + 19;
xc = centres_x(max_circ_cent);
yc = centres_y(max_circ_cent);
figure; 
subplot(1,2,1); imagesc(small_med); axis image; colormap(gray(256));
drawcircle(xc,yc,rad,'r');