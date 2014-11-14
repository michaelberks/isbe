function [boundary_out] = active_contour_ferrari(boundary_in, image_in)
%im_in is the input image
%boundary is an an Nx2 array of (x,y) coordinates of boundary points

%make sure input image is a double
if ~isa(image_in, 'double')
    image_in = double(im_in);
end

%Get number of points on the boundary
no_pts = length(boundary_in(:,1));

%Get current version of the boundary
cur_b = boundary_in;

%Get size of image
[y_max, x_max] = size(i1);

%Set constants for method
alpha = 0.2; beta = 0.2; gamma = 1.0;
ro = 1/(2*cos(2*pi/no_pts));

%Compute image gradient and maximum of gradient across image
gh = imfilter(image_in, fspecial('sobel') /8,'replicate');
gv = imfilter(image_in, fspecial('sobel')'/8,'replicate');
g_max = sqrt(max(gh(:).^2 + gv(:).^2));

do = 1; e_new = Inf; count = 0;
while (do)
    avg_dist = sum(sum(diff(cur_b).^2)) / (no_pts - 1);
    Energy_tot = 0;
    for jj = 2:no_pts-1
        v = cur_b(jj,:); %current point
        v_p = cur_b(jj - 1, :); %previous point
        v_n = cur_b(jj + 1, :); %next point
        
        for xi = -3:3
            for yi = -3:3
                x = cur_b(jj, 1) + xi;
                y = cur_b(jj, 2) + yi;
                p = [x y];
                if (x >= 1 && y >= 1 && x <= x_max && y <= y_max)
                
                    t = (v - v_p) / dist(v_p,v) + (v_n - v) / dist(v_n, v);
                    n(1) = -t(2); n(2) = t(1);
                    
                    E_con(yi+4, xi+4) = (dist(p, ro*(v_p + v_n))^2) / avg_dist;
                    E_bal(yi+4, xi+4) = n * (v - p)';
                    E_gra(yi+4, xi+4) = 1 - sqrt((gh(round(y), round(x))^2 + gv(round(y), round(x))^2)) / g_max;
                    E_ext(yi+4, xi+4) = -n*[gh(round(y), round(x)) gv(round(y), round(x))]';
                     
                else
                    E_con(yi+4, xi+4) = NaN;
                    E_bal(yi+4, xi+4) = NaN;
                    E_ext(yi+4, xi+4) = NaN;
                end                    
            end
        end
        E_con = (E_con - min(E_con(:))) / (max(E_con(:)) - min(E_con(:)));
        E_bal = (E_bal - min(E_bal(:))) / (max(E_bal(:)) - min(E_bal(:)));
        E_bal = E_bal .* E_gra;
        E_ext = (E_ext - min(E_ext(:))) / max((max(E_ext(:)) - min(E_ext(:))), g_max);
        Energy = alpha*E_con + beta*E_bal + gamma*E_ext;
        
        [E_min r_min] = min(Energy); [E_min c_min] = min(E_min); r_min = r_min(c_min);
        Energy_tot = Energy_tot + E_min;
        cur_b(jj, 1) = cur_b(jj,1) + c_min - 4;
        cur_b(jj, 2) = cur_b(jj,2) + r_min - 4;
    end
    count = count+1;
    e_old = e_new;
    e_new = Energy_tot
    do = e_new < e_old;
end

count
boundary_out = cur_b;