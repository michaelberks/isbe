%%
%Names 
m_names = get_mammo_info(dir('I:\thesis_dev\dev\masses\*.mat'));
[dummy missing_idx] = match_mammo_names('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta', m_names);
mkdir C:\isbe\asymmetry_project\data\spicules\2004_screening_processed\abnormals\
for ii = 1:179
    if ~ismember(ii, missing_idx);
         
        spicules = u_load(['M:\chen\data\masses512x512_spicules\mass_spicules', zerostr(ii,3), '.mat']);
    
        if isempty(spicules); continue; end
        
        %Load mask of source to get size
        mask = u_load([mask_dir m_names{ii} '_mask.mat']);
        [rows cols] = size(mask); 
        clear mask;

        %load mass outline
        mass_xy = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\'...
            m_names{ii} '_meta.mat']);
        mass_x = mass_xy(:,1) * cols;
        mass_y = mass_xy(:,2) * rows;
        mass_xc = round(mean(mass_x));
        mass_yc = round(mean(mass_y));

        for jj = 1:length(spicules) 
            %Translate the original spicule relative to whole mammo axis
            spicules{jj}(:,1) = spicules{jj}(:,1) + mass_xc - 256;
            spicules{jj}(:,2) = spicules{jj}(:,2) + mass_yc - 256;
        end

        save(['C:\isbe\asymmetry_project\data\spicules\2004_screening_processed\abnormals\'...
            m_names{ii} '_spic.mat'], 'spicules');
    end
end
    
%%
seg_dir = 'C:\isbe\asymmetry_project\data\segmentations\2004_screening\abnormals\';
mask_dir = 'C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\';

mammo_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals.mat');

if_plot = 0;
    
for ii = 1:length(mammo_names)
    
    ab_spics = u_load(['C:\isbe\asymmetry_project\data\spicules\2004_screening_processed\abnormals\'...
        mammo_names{ii} '_spic.mat']);
    
    if isempty(ab_spics); continue; end
    
    spicules = ab_spics;

    source_name = mammo_names{ii};
    target_name = mammo_names{ii};
    if mammo_names{ii}(4) == 'R'
        target_name(4) = 'L';
    else
        target_name(4) = 'R';
    end 

    %Load mask of source to get size
    mask = u_load([mask_dir source_name '_mask.mat']);
    [rows cols] = size(mask); 
    clear mask;

    %Load in segmentations
    seg_source = u_load([seg_dir source_name '_segmentation.mat']);
    seg_target = u_load([seg_dir target_name '_segmentation.mat']);

    source_breast = seg_source.breast_border(seg_source.breast_air,:);
    target_breast = seg_target.breast_border(seg_target.breast_air,:); 

    %compute rescale factor
    scale_factor = seg_source.size(1) / rows;
    
    if if_plot
        figure;
        a1 = subplot(1,2,1);
        plot(source_breast(:,1), source_breast(:,2)); hold on; axis equal ij;
        a2 = subplot(1,2,2);
        plot(target_breast(:,1), target_breast(:,2)); hold on; axis equal ij;
    end
    for jj = 1:length(ab_spics)
        
        %Scale the source points
        xy = ab_spics{jj}*scale_factor;

        %Compute corresponding target points
        [xyt] = select_corresponding_position(source_breast, target_breast, xy, 1);

        if if_plot
            plot(a1, xy(:,1), xy(:,2), 'g');
            plot(a2, xyt(:,1), xyt(:,2), 'r');
        end

        %Rescale and save the target points
        spicules{jj} = round(xyt / scale_factor);
    end
    
    save(['C:\isbe\asymmetry_project\data\spicules\2004_screening_processed\abnormals\'...
        target_name '_spic.mat'], 'spicules');

end
%%
%Set snake parameters
alpha = 0; %Weight to minimise length (set to zero as don't want to shrink spicule)
beta = .01; %Weight to minimise curvature
thresh = 1e-6;

search_width = 4;
norm_width = 3*search_width;
resolution = 1;
do_plot = 1;        
m_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals.mat');
%mkdir C:\isbe\asymmetry_project\data\spicule_snakes\2004_screening_processed\abnormals\;
%%
for ii = 2:10
    
    spicules = u_load(['C:\isbe\asymmetry_project\data\spicules\2004_screening_processed\abnormals\'...
        m_names{ii} '_spic.mat']);
    
    if isempty(spicules); continue; end
    
%     line_map = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\rf_prob\2004_screening_processed\abnormals\'...
%         m_names{ii} '_class.mat']);
    line_map = abs(load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_thin\2004_screening_processed\abnormals\'...
        m_names{ii} '_class.mat']));
%     line_map = -imresize(load_uint8(['C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\abnormals\line_strengths\'...
%         m_names{ii} '_data.mat']), 2, 'bilinear');
%     line_map = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\abnormals\'...
%         m_names{ii} '_lines.mat']);
    
    if do_plot
        figure; imagesc(line_map); axis image; colormap(gray(256)); hold on;
    end
    
    new_spicules = spicules;
    
    for ss = 1:length(spicules)

        %Select i-th spicule
        spicule = spicules{ss};
        
        if do_plot
            plot(spicule(:,1), spicule(:,2), 'r');
        end

        %Initialise energies
        e_old = Inf;
        e_new = 1e6;
        iter = 1;

        %Iterate until energy not sufficient reduced from previous iteration
        while e_new+thresh < e_old;

            %Display energy image (i.e line prob) and original spicule in
            %green

            %
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

                [cx, cy, cp] = improfile(line_map, [n1_x, n2_x], [n1_y, n2_y], norm_width);
                normal_p(jj, :) = cp';
                normal_x(jj, :) = cx';
                normal_y(jj, :) = cy';

            end
            normal_p(isnan(normal_p)) = 0;

            %Store old energy
            e_old = e_new;
            
             %Apply snake to spicule and compute new energy
            initial_pts = [round(norm_width/2)*ones(len,1) (1:len)'];
            [initial_pts,e_new] = mb_snake_normal(initial_pts, alpha, beta, search_width, resolution, normal_p, normal_x, normal_y);

            for kk = 1:len
                spicule(kk,1) = normal_x(kk, initial_pts(kk,1));
                spicule(kk,2) = normal_y(kk, initial_pts(kk,1));
            end

            display(['Iteration ', num2str(iter), ': energy = ', num2str(e_new)]);
            iter = iter + 1;

        end

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
        
        spicule = [xi', yi'];
        spicules{ii} = spicule;
        
        if do_plot
            plot(spicule(:,1), spicule(:,2), 'g');
        end

    end
    save(['C:\isbe\asymmetry_project\data\spicule_snakes\2004_screening_processed\abnormals\'...
        m_names{ii} '_spic.mat'], 'spicules');
end