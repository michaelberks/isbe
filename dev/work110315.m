seg_list = dir('C:\isbe\asymmetry_project\data\segmentations\2004_screening\abnormals\*LML*.mat');

for ii = 1:10
    mammo = imresize(...
        u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\' seg_list(ii).name(1:6) '.mat']),...
        [1024 nan], 'bilinear');
    [rows cols] = size(mammo);
    seg = u_load(['C:\isbe\asymmetry_project\data\segmentations\2004_screening\abnormals\' seg_list(ii).name]);
    
    breast_border = seg.breast_border;
    breast_air = seg.breast_air;
    breast_mask = roipoly(mammo, breast_border(:,1), breast_border(:,2));
    
    [yy xx] = find(breast_mask);
    breast_centre = [mean(xx) mean(yy)];
    
    %Extract pectoral region as the box bounded by the start of the breast
    %air border and the centorid of the breast
    p_rows = round(breast_centre(2));
    p_cols = round(breast_border(breast_air(1),1));
    pectoral = double(mammo(1:p_rows, 1:p_cols));
    
    %Filter the pectoral region and discard high/low pixels
    pectoral_med = medfilt2(pectoral, [9 9]);
    pectoral_mask = pectoral_med > 240 | pectoral_med < 2;
    pectoral_med(pectoral_mask) = 0;
    
    %Perform canny edge detection on pectoral region
    [pectoral_canny ori] = canny_edge(pectoral_med, [], 1, 0.90, .6);
    
    %Discard edges that aren't approximatelt aligned on an upperleft
    %diagonal
    ignore_ori = ori > pi/4 | ori < pi/24;
    pectoral_ori = pectoral_canny;
    pectoral_ori(ignore_ori) = 0;
    
    %Take hough transformation of edges to locate most likely staright line
    [line_scores rhos thetas] = hough_line(pectoral_ori);
    [max_rhos rho_idx] = max(line_scores);
    [dummy theta_idx] = max(max_rhos);
    
    %Convert from rho,theta parametisation to y = mx + c;
    theta = thetas(theta_idx);
    rho = rhos(rho_idx(theta_idx));
    m = tan(theta); 
    c = rho / cos(theta);
    
    %Work out points on pectoral edge at top and bottom of image
    y1 = 1; x1 = (y1 - c)/m;
    y2 = rows; x2 = (y2 - c)/m;
    pectoral_edge = [x1 y1; x2 y2];
    
    figure; 
    subplot(1,2,1); imagesc(pectoral_ori); axis image; hold on; colormap(gray(256));
    title(seg_list(ii).name(1:6));
    
    %--------------------------
    [yy xx] = find(pectoral_ori);
    x_norm = -sin(theta);
    y_norm = cos(theta);

    edge_dists = abs((xx - x1)*x_norm + (yy - y1)*y_norm);
    edge_oris = pi/2 - ori(pectoral_ori);

    valid_dists = edge_dists < 20;
    valid_ori = abs(mb_mod(edge_oris + theta, pi)) < pi/12;
    valid_edges = valid_dists & valid_ori;

    xx = xx(valid_edges);
    yy = yy(valid_edges);
    plot(xx, yy, 'm.');

    N = length(xx);
    n = 5;
    x_plot = 1:p_cols;
    
    max_score = 0;
    kk = 0;
    while kk < 100
        r_idx = randperm(N);
        x_fit = xx(r_idx(1:n));
        if all(diff(sort(x_fit)))
            kk = kk+1;
            y_fit = yy(r_idx(1:n));


            %p = polyfit(x_fit,y_fit,n-1);
            %y_plot = polyval(p, x_plot);
            
            y_plot = interp1(x_fit, y_fit, x_plot, 'cubic');


            curr_score = 0;
            for jj = 1:N
                r2 = (xx(jj) - x_plot).^2 + (yy(jj) - y_plot).^2;
                curr_score = curr_score + (min(r2) < 10);
            end
            if curr_score > max_score;
                x_max = x_plot;
                y_max = y_plot;
                max_score = curr_score;
            end
        end
    end
    
    %--------------------------
    
    
    plot([x1 x2], [y1 y2], 'r');
    plot(x_max, y_max, 'g');
    subplot(1,2,2); imagesc(pectoral); axis image; hold on; colormap(gray(256));
    title(seg_list(ii).name(1:6));
    plot([x1 x2], [y1 y2], 'r');
    plot(x_max, y_max, 'g');
end