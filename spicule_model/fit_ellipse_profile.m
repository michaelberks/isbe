function [profile_data] = fit_ellipse_profile...
    (s_profile, offset, width, jump, n_pts, E, B, solutions)

%local ellipse profile fitting
    if_plot = 1;
    
    e_profile = enhance_spicules(s_profile, 20, 20);
    [r c] = size(e_profile);
    if offset + width < r/2;
        if if_plot
            f1 = figure;
            imagesc(e_profile); axis image; colormap gray; hold on;
        end
        
        mid = ceil((r+1) / 2);
        
        if nargin < 6
            [E B] = calculate_error;            
        end
        if nargin < 8
            solutions = get_solutions(E, B);
        end
        
        [x0_local w_local b_local] = get_local_profile(solutions, jump);
        w_params = ...
            interp1(1:c, w_local, linspace(1, c, n_pts))';
        
        w_global = ...
            interp1(linspace(1, c, n_pts), w_params, 1:c)';
        
        b_params = ...
            interp1(1:c, b_local, linspace(1, c, n_pts))';
        
        b_global = ...
            interp1(linspace(1, c, n_pts), b_params, 1:c)';
        
        profile_data.x0_local = x0_local;
        profile_data.w_local = w_local;
        profile_data.w_params = w_params;
        profile_data.b_local = b_local;
        profile_data.b_params = b_params;
        profile_data.errors = E;
        profile_data.solutions = solutions;
        
        if if_plot
            %figure(f1);
            figure
            imagesc(e_profile); axis image; colormap gray; hold on;
            %Plot shortest path on spicule and solution widths
            plot(x0_local(:,1), x0_local(:,2)+mid, 'y');
            plot(x0_local(:,1), x0_local(:,2)+mid-w_global, 'g');
            plot(x0_local(:,1), x0_local(:,2)+mid+w_global, 'r');
            
            w_matrix = repmat(w_local', r, 1);
            b_matrix = repmat(b_local', r, 1);
            offsets = repmat([-40:40]', 1, c);
            offsets(abs(offsets) >= w_matrix) = w_matrix(abs(offsets) >= w_matrix);        
            nprofile = b_matrix .* sqrt(1 - (offsets.^2 ./ w_matrix.^2));
            figure;
            imagesc([e_profile;nprofile]); axis image; colormap gray; hold on;    
        end
        
        
        
    else
        display('offset and width values are too large for profile');
        solutions = 0;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Auxillary functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [E B] = calculate_error
   
        stretch_e = imresize(e_profile, [5*r c], 'bilinear');
        
        if if_plot
            plot(1:c, mid - offset, 'w');
            plot(1:c, mid + offset, 'w');
        end
        
        for ii = mid - offset:mid + offset;
            for jj = 1:width;
                
                w = 5*jj;
                x0 = 5*ii;
                
                x = repmat([-w:w]', 1, c);
                N = 2*w + 1;
                
                p = stretch_e(x0-w:x0+w, :);
                b = calculate_b(p, x, w); %calculate_b
                b(~b) = 0.01; %supress divide by 0;
                b_mat = repmat(b, N, 1);
                
                E(ii+offset-mid+1, jj, :) = ...
                    sqrt(sum((p - g(x, 0, b_mat, w)).^2)) ./ (N*b);
                B(ii+offset-mid+1, jj, :) = b;
            end
        end
    end     
        
    
    function ellipse = g(x, x0, b, w)
        ellipse = b .* sqrt(1 - (x - x0).^2 / w.^2);
    end

    function b = calculate_b(p, x, w)
        b = sum(sqrt(1 - (x / w).^2) .* p) ./ sum(1 - (x / w).^2);
    end

    function solutions = get_solutions(E, B)
        
        [x0, w, len] = size(E);
        
        %find first solution
        %[min1 ind1] = min(reshape(E, [x0*w, len]));
        %[i_x0 i_w] = ind2sub([x0 w], ind1);
        %solutions(1).profile = [i_x0'+20, i_w'];
        
        for ii = 1:len
            er = E(:,:,ii);
            
            %find first solution
            [min1 ind1] = min(er(:));
            [i_x0 i_w] = ind2sub(size(er), ind1);
            i_b = B(ind1)';
            sol(1).profile(ii,:) = [i_x0-offset-1 i_w i_b];
            sol(1).error(ii) = min1;
            if if_plot, plot(ii, i_x0+mid-offset-1, 'r.'); end;
            
            %find second solution
            er(max(i_x0-i_w, 1):min(i_x0+i_w,end), :) = NaN;
            [min1 ind1] = min(er(:));
            [i_x0 i_w] = ind2sub(size(er), ind1);
            i_b = B(ind1)';
            sol(2).profile(ii,:) = [i_x0-offset-1 i_w i_b];
            sol(2).error(ii) = min1;
            if if_plot, plot(ii, i_x0+mid-offset-1, 'g.'); end
            
            %find third solution
            er(max(i_x0-i_w, 1):min(i_x0+i_w,end), :) = NaN;
            [min1 ind1] = min(er(:));
            [i_x0 i_w] = ind2sub(size(er), ind1);
            i_b = B(ind1)';
            sol(3).profile(ii,:) = [i_x0-offset-1 i_w i_b];
            sol(3).error(ii) = min1;
            if if_plot, plot(ii, i_x0+mid-offset-1, 'b.'); end
            
            %find fourth solution
            er(max(i_x0-i_w, 1):min(i_x0+i_w,end), :) = NaN;
            [min1 ind1] = min(er(:));
            [i_x0 i_w] = ind2sub(size(er), ind1);
            i_b = B(ind1)';
            sol(4).profile(ii,:) = [i_x0-offset-1 i_w i_b];
            sol(4).error(ii) = min1;
            if if_plot, plot(ii, i_x0+mid-offset-1, 'y.'); end
        end
        
        solutions = sol;
        
    end

    function [x0_local w_local b_local] = get_local_profile(solutions, jump, n)

        %figure, imagesc(e_profile); colormap gray; hold on; axis image;
        hold on;

        if nargin < 3
            n = c;
        end
        N = 4*n;

        %Set node x-co-ordinates
        nodes_x = repmat(1:n, 4, 1);
        nodes_x = nodes_x(:);

        %Set node y-co-ordinates
        nodes_y(1,:) = solutions(1).profile(1:n,1)';
        nodes_y(2,:) = solutions(2).profile(1:n,1)';
        nodes_y(3,:) = solutions(3).profile(1:n,1)';
        nodes_y(4,:) = solutions(4).profile(1:n,1)';  
        nodes_y = nodes_y(:);

        %Set node indices for first level of links - idx1(i) -> idx2(i)
        idx1 = repmat(1:N-4, 4, 1);
        idx1 = idx1(:);

        idx2 = repmat(reshape(5:N, 4, []), 4, 1);
        idx2 = idx2(:);

        %Set node indices for additional levels of links
        for ii = 2:jump
            idx1 = [idx1; idx1(1:end-16)];
            idx2 = [idx2; idx2(17:end)];
        end

        % calculate weights as distance between linked nodes
        weights = (nodes_x(idx1)-nodes_x(idx2)).^2 + ...
                        (nodes_y(idx1)-nodes_y(idx2)).^2;

        %Set adjacency matrix
        adj = sparse(zeros(N+2, N+2));
        adj(sub2ind([N+2, N+2], idx1+1, idx2+1)) = weights;

        %update nodes and adjacency matrix to include start and end nodes 
        nodes_x = [0; nodes_x; n+1];
        nodes_y = [0; nodes_y; 0];

        adj(1,2:5) = 1;
        adj(N-2:N+1, N+2) = 1;

        %Can plot links between vertices though looks very messy
        %figure, gplot(adj, [nodes_x nodes_y]); hold on;
        
        %Use dijkstra algorithm to calculate shortest paths between start node
        %and all other nodes
        [dp pred] = dijkstra_sp(adj,1);

        %Obtain shortest path from start to end node by following predecessor
        %nodes
        shortest_path(1) = pred(N+2);
        ii = 1;
        while (shortest_path(ii) > 1)
            shortest_path(ii+1) = pred(shortest_path(ii));
            ii = ii +1;
        end
        shortest_path = fliplr(shortest_path);
        shortest_path(1) = [];

        %Obtain x0 and width values of shortest path nodes :- the local
        %solutions to the ellipse model
        x0_x = nodes_x(shortest_path);
        x0_y = nodes_y(shortest_path);
        x0_local = [[1:n]' round(interp1(x0_x, x0_y, 1:n))'];
        
        %plot(x0_x, x0_y, 'g');

        sol(1).shortest_path = shortest_path;
        sol(1).x0_x = x0_x;
        sol(1).x0_y = x0_y;
        sol(1).x0_local = x0_local;
        n_sols = 1;
        jj = 1;
        total_sols = 6;
        while n_sols < total_sols
            
            try
                dists = diff(x0_x).^2 + diff(x0_y).^2;
                [mm idx] = max(dists);
                adj(shortest_path(idx), shortest_path(idx+1)) = 0;        
                
                [dp pred] = dijkstra_sp(adj,1);
                clear shortest_path;
                shortest_path(1) = pred(N+2);

                ii = 1;
                while (shortest_path(ii) > 1)
                    shortest_path(ii+1) = pred(shortest_path(ii));
                    ii = ii +1;
                end
                shortest_path = fliplr(shortest_path);
                shortest_path(1) = [];

                %Obtain x0 and width values of shortest path nodes :- the local
                %solutions to the ellipse model
                x0_x = nodes_x(shortest_path);
                x0_y = nodes_y(shortest_path);
                x0_local = [[1:n]' round(interp1(x0_x, x0_y, 1:n))'];
                
                for ii = 1:n_sols
                    mean_diff(ii) = ...
                        mean(sqrt(sum((sol(ii).x0_local' - x0_local').^2)));
                end
                
                if min(mean_diff) > 1
                    %figure, gplot(adj, [nodes_x nodes_y]); hold on;
                    %plot(sol(n_sols).x0_x, sol(n_sols).x0_y, 'r');
                    %plot([sol(n_sols).x0_x(idx) sol(n_sols).x0_x(idx+1)],...
                    %    [sol(n_sols).x0_y(idx) sol(n_sols).x0_y(idx+1)], 'y');
                    %plot(x0_x, x0_y, 'g');
                    n_sols = n_sols + 1;
                    sol(n_sols).shortest_path = shortest_path;
                    sol(n_sols).x0_x = x0_x;
                    sol(n_sols).x0_y = x0_y;
                    sol(n_sols).x0_local = x0_local;
                end
            catch
                total_sols = n_sols;
                display('error occurred finding new solution');
                display(lasterr);
            end
        end
        
        for ii = 1:n_sols
            mean_x(ii) = mean(abs(sol(ii).x0_y));
        end
        [mm idx] = min(mean_x);
        shortest_path = sol(idx).shortest_path;
        x0_x = sol(idx).x0_x;
        x0_y = sol(idx).x0_y;
        x0_local = sol(idx).x0_local;
        clear sol;
        
        %Get widths of shortes paths nodes
        widths(1,:) = solutions(1).profile(1:n,2)';
        widths(2,:) = solutions(2).profile(1:n,2)';
        widths(3,:) = solutions(3).profile(1:n,2)';
        widths(4,:) = solutions(4).profile(1:n,2)';  
        widths = widths(:);
        w_local = round(interp1(x0_x, widths(shortest_path-1), 1:c))';
        
        brights(1,:) = solutions(1).profile(1:n,3)';
        brights(2,:) = solutions(2).profile(1:n,3)';
        brights(3,:) = solutions(3).profile(1:n,3)';
        brights(4,:) = solutions(4).profile(1:n,3)';  
        brights = brights(:);
        b_local = round(interp1(x0_x, brights(shortest_path-1), 1:c))';
        
    end
 
end

    %function profile = calculate_error
    %    
    %    figure, imagesc(e_profile); axis image; colormap gray; hold on;
    %    [r c] = size(e_profile);
    %    
    %    for ii = 1:c;
    %        for w = 1:20;
    %            for x0 = 20:60;
    %                
    %                x = [-w:0.2:w]';
    %                N = 10*w + 1;
    %                p = improfile(e_profile, [ii, ii], [x0-w, x0+w], N);
    %                b = calculate_b(p, x, w); %calculate_b
    %                
    %                E(x0-19, w, ii) = sqrt(sum((p - g(x, 0, b, w)).^2)) / (N*b);
    %                
    %            end
    %        end                       
    %    end
    %    [min1 ind1] = min(reshape(E, [41*20, c]));
    %    [i_x0 i_w] = ind2sub([41 20], ind1);
    %    profile = [i_x0' + 19, i_w'];
    %    plot(1:c, i_x0 - i_w, 'r');
    %    plot(1:c, i_x0 + i_w, 'g');
    %    plot(1:c, i_x0, 'y');
    %end