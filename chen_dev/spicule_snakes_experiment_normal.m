%%
%--------------------------------------------------------------------------
%-- Experimental code using snakes to update annotated spicule position
%--------------------------------------------------------------------------

%Load in probability images and spicules
mass_idx = 56;
prob_image = 1-u_load(['M:\chen\data\predict_masses512x512_chen\probability_image', zerostr(mass_idx,3), '.mat']);
%prob_image = u_load(['M:\chen\data\masses512x512_nms\mass_nms', zerostr(mass_idx,3), '.mat']);
spicules = u_load(['M:\chen\data\masses512x512_spicules\mass_spicules', zerostr(mass_idx,3), '.mat']);
% mass_image = u_load(['M:\chen\data\masses512x512\mass', zerostr(mass_idx,3), '.mat']);
% save_dir = 'M:\chen\data\masses512x512_spicules\mass_newspicules_chen';

% set(0,'DefaultFigureWindowStyle','docked')


if ~isempty(spicules)
    %Set snake parameters
    alpha = 0; %Weight to minimise length (set to zero as don't want to shrink spicule)
    beta = .01; %Weight to minimise curvature
    search_width = 4;
    norm_width = 3*search_width;
    resolution = 1;

    feat_img = prob_image; %Feature image = line probability
    %
    spic_axes = cell(length(spicules), 1);
    for ii = 1:length(spicules)

        %Select i-th spicule
        spicule = spicules{ii};

        %Initialise energies
        e_old = Inf;
        e_new = 1e6;
        thresh = 1e-6;
        iter = 1;

        
        %Iterate until energy not sufficient reduced from previous iteration
        while e_new+thresh < e_old;

            %Display energy image (i.e line prob) and original spicule in green
            figure('windowstyle', 'docked'); imagesc(feat_img); axis image; colormap(gray(256)); hold on;
            plot(spicule(:,1), spicule(:,2), 'g', 'LineWidth',2);
            
            %>>>>>code
            len = size(spicule,1);
            normal_p = zeros(len, norm_width);
            normal_x = zeros(len, norm_width);
            normal_y = zeros(len, norm_width);

            %Compute the normal vectors at each point on the inner border
            [fx, fy] = gradient(spicule);

            %normalise fy
            fy = fy ./ [sqrt(sum(fy.^2, 2)), sqrt(sum(fy.^2, 2))];

            %Compute normal profiles of the image at every point
            for jj = 1:len %= number of rows in skin_air

                n1_x = spicule(jj, 1) - norm_width*fy(jj, 2);
                n1_y = spicule(jj, 2) + norm_width*fy(jj, 1);
                n2_x = spicule(jj, 1) + norm_width*fy(jj, 2);
                n2_y = spicule(jj, 2) - norm_width*fy(jj, 1);

                [cx, cy, cp] = improfile(feat_img, [n1_x, n2_x], [n1_y, n2_y], norm_width);
                normal_p(jj, :) = cp';
                normal_x(jj, :) = cx';
                normal_y(jj, :) = cy';

            end
            normal_p(isnan(normal_p)) = 0;
            
                     %Store old energy
            e_old = e_new;
            initial_pts = [round(norm_width/2)*ones(len,1) (1:len)'];
            [initial_pts,e_new] = mb_snake_normal(initial_pts, alpha, beta, search_width, resolution, normal_p, normal_x, normal_y);

            for kk = 1:len
                spicule(kk,1) = normal_x(kk, initial_pts(kk,1));
                spicule(kk,2) = normal_y(kk, initial_pts(kk,1));
            end
%             %%%%%%%%%%%%%%%% 
%             %Store old energy
%             e_old = e_new;
% 
%             %Apply snake to spicule and compute new energy
%             [spicule, e_new] = mb_snake(spicule, alpha, beta, max_delta_x, resol_x, max_delta_y, resol_y, feat_img);

            %Display new spciule in blue/red
            plot(spicule(:,1), spicule(:,2), 'bx', 'LineWidth',2);
            plot(spicule(:,1), spicule(:,2), 'r', 'LineWidth',2);
            % Plot the normal points
%             plot(normal_x, normal_y, 'b.');
            %Print out the new energy and increment the iteration count
            display(['Iteration ', num2str(iter), ': energy = ', num2str(e_new)]);
            spic_axes{ii}(iter) = gca;
            iter = iter + 1;
            
        end
        linkaxes(spic_axes{ii});
        
        dxy = abs(spicule(1,:)-spicule(end,:));
        dx = dxy(1);
        dy = dxy(2);
        if dx>=dy
            sort_spicule = sortrows(spicule, 1);
            xx=sort_spicule(:, 1);
            yy=sort_spicule(:, 2);
            xi = round(sort_spicule(1,1):sort_spicule(end,1));
            yi = round(interp1(xx,yy,xi, 'spline'));
        else
           sort_spicule = sortrows(spicule, 2);
           xx=sort_spicule(:, 1);
            yy=sort_spicule(:, 2);
            yi = round(sort_spicule(1,2):sort_spicule(end,2));
            xi = round(interp1(yy,xx,yi, 'spline'));
        end  
        spicule=[xi', yi'];
        plot(spicule(:,1), spicule(:,2), 'c', 'LineWidth',2);
      
        new_spicules{ii}=spicule;
        
    end
end
% spicules = new_spicules;
% save([save_dir, 'spicules', zerostr(patch_num, 3)], 'spicules');
return

