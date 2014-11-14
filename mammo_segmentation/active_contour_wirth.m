function [boundary_out] = active_contour_wirth(boundary_in, bin_im, no_its)
% implementation of the active contour algorithm

no_pts = length(boundary_in(:,1));
cur_b = boundary_in;
[y_max, x_max] = size(bin_im);

alpha = 2.2; beta = 2.5; gamma = 2.0;

for ii = 1:no_its
    avg_dist = sqrt( mean(diff(cur_b))*mean(diff(cur_b))' );
    for jj = 2:no_pts-1
        if jj > 1, v_p = cur_b(jj - 1, :); end %previous point
        if jj < no_pts, v_n = cur_b(jj + 1, :); end %next point
        %v_mid = 0.5*(v_p + v_n); %mid point
        
        Energy(7,7) = 0;
        for xi = -3:3
            for yi = -3:3
                x = cur_b(jj, 1) + xi;
                y = cur_b(jj, 2) + yi;
                if (x >= 1 && y >= 1 && x <= x_max && y <= y_max)
                
                    %E_int = ( alpha * (v_p-[x y])*(v_p - [x y])' + beta * (v_mid-[x y])*(v_mid-[x y])' ) / 2;

                    %E_line = 1000*bin_im(round(y), round(x)); %check x and y

                    %E_edge = 0;

                    %Energy(xi+4, yi+4) =  E_int + E_line + E_edge;
                    if jj == 1,
                        E_cont = abs(avg_dist - dist(v_n,[x y]));
                        E_curvX = 0;
                        E_curvY = 0;
                    elseif jj == no_pts
                        E_cont = abs(avg_dist - dist(v_p,[x y]));
                        E_curvX = 0;
                        E_curvY = 0;
                    else
                        E_cont = abs(avg_dist - dist(v_p,[x y])) + abs(avg_dist - dist(v_n,[x y]));
                        E_curvX = (v_p(1) - x) / dist(v_p,[x y]) + (v_n(1) - x) / dist(v_n,[x y]);
                        E_curvY = (v_p(2) - y) / dist(v_p,[x y]) + (v_n(2) - y) / dist(v_n,[x y]);
                    end
                    E_image = 0;
                    
                    Energy(yi+4, xi+4) = alpha*E_cont;% + beta*(E_curvX^2 + E_curvY^2) + gamma*E_image; 
                else
                    Energy(yi+4, xi+4) = 1000;
                end                    
            end
        end
        Energy;
        [E_min r_min] = min(Energy); [E_min c_min] = min(E_min); r_min = r_min(c_min);
        cur_b(jj, 1) = cur_b(jj,1) + c_min - 4;
        cur_b(jj, 2) = cur_b(jj,2) + r_min - 4;
    end
end

boundary_out = cur_b;
        