function [r_x, r_y, x0 , y0, x_axis] = fit_ellipse(outline)

% Calculate the best fit ellipse to an irregular polygonal boundary
    
    %to complete boundary need to append first point to the end, if not
    %done already
    %if ~sum(outline(1,:) - outline(end,:))
    %    outline(end+1,:) = outline(1,:);
    %end
    
    diffs = diff(outline([1:end,1],:));
    dists = sqrt(sum(diffs.^2, 2));
    dists = dists / sum(dists(:));
    theta_diff = 2*pi*dists;
    theta = 2*pi*cumsum(dists);
    theta = theta([end,1:end-1]);
    
    x0_int = outline(:,1);
    x0 = int_calc(x0_int, theta_diff) / (2*pi);
    y0_int = outline(:,2);
    y0 = int_calc(y0_int, theta_diff) / (2*pi);
    a1_int = outline(:,1) .* cos(theta);
    a1 = int_calc(a1_int, theta_diff) / pi;
    b1_int = outline(:,1) .* sin(theta);
    b1 = int_calc(b1_int, theta_diff) / pi;
    c1_int = outline(:,2) .* cos(theta);
    c1 = int_calc(c1_int, theta_diff) / pi;
    d1_int = outline(:,2) .* sin(theta);
    d1 = int_calc(d1_int, theta_diff) / pi;
    
    alpha = a1^2 + b1^2 + c1^2 + d1^2;
    beta = a1*d1 - c1*b1;
    r_x = sqrt((alpha + sqrt(alpha^2 - 4*beta^2)) / 2);
    r_y = sqrt(2*beta^2 / (alpha + sqrt(alpha^2 - 4*beta^2)));
    
    x_axis(1) = r_x*a1 - r_y*d1;
    x_axis(2) = r_x*c1 + r_y*b1;
    x_axis = x_axis / sqrt(sum(x_axis.^2));
end

function int_sum = int_calc(ints, thetas)
    int_sum = sum((ints + ints([2:end, 1])) .* thetas / 2);
end
    