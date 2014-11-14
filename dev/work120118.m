nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');
num_nf = length(nf_files);

%%
ori_bins = linspace(0, pi, 90);
for nn = 1:num_nf
    image_path = ['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp'];
    nailfold = imread(image_path);
    figure; imgray(nailfold);
    
    %Load in the annotated vessels
    v = read_vessels_from(['C:\isbe\nailfold\images\anonymous_oct\annotations_qt\' nf_files(nn).name]);
    vessels = [];
    ori_counts = zeros(1,90);
    for ii = 1:length(v)
        if size(v{ii},1) > 1;       
            %Discard duplicate points
            keep = [true; any(diff(v{ii}),2)];
            vessel = [v{ii}(keep,1) v{ii}(keep,2)];
            vessel = [vessel; vessel(1,:)];

            %Sample points evenly along the vessel
            dists = cumsum([0; sqrt(sum(diff(vessel).^2, 2))]);
            vessels{end+1,1} = interp1(dists(1:end-1), vessel(1:end-1,:), linspace(0, dists(end-1), floor(dists(end-1))), 'spline'); %#ok
            vessel_circ = interp1(dists, vessel, linspace(0, dists(end), floor(dists(end))), 'linear');
            plot(vessel_circ(:,1),vessel_circ(:,2), 'm');
            plot(vessels{end,1}(:,1),vessels{end,1}(:,2), 'r');
            
%             num_pts = size(vessels{end,1},1);
%             i1 = ceil(num_pts/3);
%             i2 = floor(2*num_pts/3);
%             
%             v1 = vessels{end,1}(1:i1,:);
%             v2 = vessels{end,1}(i2:end,:);
%             
%             p1 = polyfit(v1(:,1), -v1(:,2), 1);
%             v1_ori = atan(p1(1));
%             %text(v1(1,1), v1(1,2), num2str(round(180*v1_ori/pi)));
%             
%             cxy = mean(vessels{end,1});
%             pc = polyfit(vessels{end,1}(:,1), vessels{end,1}(:,2), 1);
%             c_ori = -atan(pc(1));
%             
%             plot(cxy(1) + [-10 10], cxy(2) + pc(1)*[-10 10]);
%             text(v1(1,1), v1(1,2), num2str(round(180*c_ori/pi)));
            
             
            [r_x, r_y, x0 , y0, x_axis] = fit_ellipse(vessel_circ);
            [ex ey] = ellipse(r_x,r_y,x0,y0,x_axis);
            plot(ex, ey, 'g');
            
            smooth_diff_vec = conv2(diff(vessels{end,1}), ones(5,1)/5, 'valid');
            ori_counts_ii = hist(mod(atan(smooth_diff_vec(:,2)./smooth_diff_vec(:,1)), pi), ori_bins);
            
            %figure; bar(ori_bins, ori_counts_ii);
            
            ori_counts = ori_counts + ori_counts_ii;
        end
    end
    %figure; bar(ori_bins, ori_counts);
    num_vessels = length(vessels);
end
%%
for nn = 1:12
    image_path = ['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp'];
    nailfold = imread(image_path);
    bw = nailfold < 250;
    bw_edge = bw & ~imerode(bw, strel('disk', 1));
    bw_edge([1 end], :) = 0;
    bw_edge(:,[1 end]) = 0;
    [H,T,R] = hough(bw_edge,'RhoResolution',2,'Theta',-90:1:89);
    peaks = houghpeaks(H, 10, 'thresh', 0);
    lines = houghlines(bw_edge,T,R,peaks(1:10,:));
    
    oris = mb_mod(1 - peaks(:,2), 180);
    oris(abs(oris) > 45) = mb_mod(oris(abs(oris) > 45) - 90, 180);
    
    rot_angle = mode(oris);
    
    % Display the Hough matrix.
    figure; 
    subplot(2,1,1); imgray(imrotate(nailfold, -rot_angle));
    subplot(2,1,2); imgray(bw_edge);

    max_len = 0;
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];
       plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

       % Plot beginnings and ends of lines
       plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

       % Determine the endpoints of the longest line segment
       len = norm(lines(k).point1 - lines(k).point2);
       if ( len > max_len)
          max_len = len;
          xy_long = xy;
       end
    end

    % highlight the longest line segment
    plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','blue');
    
    
    %figure; hist(oris, -45:45);
    xlabel(['Angle = ' num2str(rot_angle) '^o']);
end
%%
for nn = 1:num_nf
    image_path = ['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp'];
    nailfold = imread(image_path);
    
    
    bw = nailfold < 250;
    bw_edge = bw & ~imerode(bw, strel('disk', 1));
    bw_edge([1 end], :) = 0;
    bw_edge(:,[1 end]) = 0;
    [H,T,R] = hough(bw_edge,'RhoResolution',2,'Theta',-90:1:89);
    peaks = houghpeaks(H, 10, 'thresh', 0);
    lines = houghlines(bw_edge,T,R,peaks(1:10,:));
    
    oris = mb_mod(1 - peaks(:,2), 180);
    oris(abs(oris) > 45) = mb_mod(oris(abs(oris) > 45) - 90, 180);
    
    rot_angle = mode(oris);
    rot_mat = [cosd(rot_angle) sind(rot_angle); -sind(rot_angle) cosd(rot_angle)];
    
    nailfold_rot = imrotate(nailfold, -rot_angle, 'nearest');
    
    [r1 c1] = size(nailfold);
    [r2 c2] = size(nailfold_rot);
    
    figure;
    h1 = subplot(2,1,1); imgray(nailfold);
    h2 = subplot(2,1,2); imgray(nailfold_rot);
    
    %Load in the annotated vessels
    v = read_vessels_from(['C:\isbe\nailfold\images\anonymous_oct\annotations_qt\' nf_files(nn).name]);
    vessels = [];
    vessels_rot = [];
    
    for ii = 1:length(v)
        if size(v{ii},1) > 1;       
            %Discard duplicate points
            keep = [true; any(diff(v{ii}),2)];
            vessel = [v{ii}(keep,1) v{ii}(keep,2)];

            %Sample points evenly along the vessel
            dists = cumsum([0; sqrt(sum(diff(vessel).^2, 2))]);
            vessel = interp1(dists, vessel, linspace(0, dists(end), floor(dists(end))), 'spline');
            vessels{end+1,1} = vessel; %#ok          
            vessel_rot = bsxfun(@plus, bsxfun(@minus, vessel, [c1 r1]/2)*rot_mat, [c2 r2]/2);
            plot(h1, vessel(:,1),vessel(:,2), 'r');
            plot(h2, vessel_rot(:,1),vessel_rot(:,2), 'r');
        end
    end
    num_vessels = length(vessels);
    
    %Now choose each vessel apex as the point with highest y-val (I know this
    %may not strictly be accurate given local deviations, rotations etc - but
    %it appears to work more robustly than more complex methods such as
    %searching for points of maximum curvature etc)
    vessel_apex = zeros(num_vessels, 2);
    vessel_tops = cell(num_vessels, 1);
    
    for ii = 1:num_vessels
        [y_min min_idx] = min(vessels{ii}(:,2));
        vessel_apex(ii,:) = [vessels{ii}(min_idx,1) y_min];

        vessel_l = vessels{ii}(min_idx:-1:1,:);
        dists = cumsum([0; sqrt(sum(diff(vessel_l).^2, 2))]);
        vessel_l = interp1(dists, vessel_l, 2:2:30, 'spline');

        vessel_r = vessels{ii}(min_idx:end,:);
        dists = cumsum([0; sqrt(sum(diff(vessel_r).^2, 2))]);
        vessel_r = interp1(dists, vessel_r, 2:2:30, 'spline');

        vessel_tops{ii} = [flipud(vessel_l); vessels{ii}(min_idx,:); vessel_r];
        if vessel_tops{ii}(1,1) > vessel_tops{ii}(end,1)
            vessel_tops{ii} = flipud(vessel_tops{ii});
        end


        plot(h1, vessels{ii}(:,1), vessels{ii}(:,2)); 
        plot(h1, vessel_apex(ii,1), vessel_apex(ii,2), 'rx');
        %plot(vessel_l(:,1), vessel_l(:,2), 'g.');
        %plot(vessel_r(:,1), vessel_r(:,2), 'y.');
        plot(h1, vessel_tops{ii}(:,1), vessel_tops{ii}(:,2), 'm-.');
        plot(h1, vessel_tops{ii}(1,1), vessel_tops{ii}(1,2), 'gx');
        plot(h1, vessel_tops{ii}(end,1), vessel_tops{ii}(end,2), 'yx');
    end
end