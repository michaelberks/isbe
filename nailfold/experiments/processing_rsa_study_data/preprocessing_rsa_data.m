%Script for processing the data acquired from the pilot images sent out in
%GD's imaging study - 25 images, marked up by a set of experts using PT's
%software
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%%
%Get list of images
rsa_dir = 'C:\isbe\nailfold\data\rsa_study\';
rsa_images = dir([rsa_dir 'images\*.png']);
num_images = length(rsa_images);

%Get list of markers
markers = dir([rsa_dir '\markup\']);
markers = markers(3:end);
keep = false(length(markers),1);
for i_ma = 1:length(markers)
    keep(i_ma) = markers(i_ma).isdir;
end
markers = markers(keep);


num_markers = length(markers);

apex_dir = [rsa_dir 'apexes\'];
%--------------------------------------------------------------------------
%% Make mask enhanced copies of the images
mkdir([rsa_dir 'images_enhanced']);

for i_im = 1:num_images
    
    %load in image
    im_num = rsa_images(i_im).name(1:end-4);   
    nailfold = imread([rsa_dir 'images\' im_num '.png']);

    if size(nailfold,3)==3
        nailfold = nailfold(:,:,1);
    end
    nailfold = double(nailfold);
    [mask] = make_nailfold_mosaic_mask(nailfold, 240, 20);
    clims = [min(nailfold(mask)) max(nailfold(mask))];
    write_im_from_colormap(nailfold, [rsa_dir 'images_enhanced\' im_num '.png'], gray(256), clims);     
end
%--------------------------------------------------------------------------
%% Initial inspection of images and apexes
%This code displays the images and the apexes that each marker has marked
for i_im = 1:2%num_images
    
    %load in image
    im_num = rsa_images(i_im).name(1:end-4);   
    im = imread([rsa_dir 'images\' im_num '.png']);
    
    %get mask of image to set contrast
    mask = im < 250;
    mask = imerode(mask, strel('disk', 10));
    g_max = max(im(mask));
    g_min = min(im(mask));
    
    %Display image
    figure; imgray(im); caxis([g_min g_max]);
    
    %Loop through markers
    for i_ma = 1:num_markers
        
        %Check if marker has marked this image
        markup_list = dir([rsa_dir 'markup\' markers(i_ma).name '\*#' im_num '*']);         
        if isempty(markup_list); continue; end
        
        %Read in vessel markup data - where multiple markup files exist,
        %use the last one
        markup = read_markup_from([rsa_dir 'markup\' markers(i_ma).name '\' markup_list(end).name]);
        
        num_vessels = length(markup.vessels);
        
        for i_v = 1:num_vessels
            
            %Check this is a valid vessel
            anchor_xy = markup.vessels(i_v).anchor;
            
            if isempty(anchor_xy); continue; end                           
                
            %Check if distal
            is_distal = markup.vessels(i_v).ncm_vessel_properties.is_distal;
            
            if is_distal
                
                %Find nearest marked apex to the anchor point
                num_apices = length(markup.vessels(i_v).apices);
                if num_apices == 1
                    
                    if isempty(markup.vessels(i_v).apices.inner_point)
                        %plot the anchor and continue
                        text(anchor_xy(1), anchor_xy(2), num2str(i_ma), 'color', 'b');
                        continue;
                    else
                        apex_xy = (markup.vessels(i_v).apices.inner_point +...
                            markup.vessels(i_v).apices.outer_point) / 2;
                    end
                else
                    min_dist = inf;
                    for i_a = 1:num_apices
                        apex_xyi = (markup.vessels(i_v).apices(i_a).inner_point +...
                            markup.vessels(i_v).apices(i_a).outer_point) / 2;
                        dist = sum((anchor_xy - apex_xyi).^2);
                        if dist < min_dist
                            apex_xy = apex_xyi;
                            min_dist = dist;
                        end
                    end
                end
                
                %plot the apex
                text(apex_xy(1), apex_xy(2), num2str(i_ma), 'color', 'g');
            else
                %plot the anchor
                text(anchor_xy(1), anchor_xy(2), num2str(i_ma), 'color', 'r');
            end
        end
        
    end
end
%
clc
for i_m = 1:num_markers
    display(['Marker ' num2str(i_m) ': ' markers(i_m).name]);
end
%%
%--------------------------------------------------------------------------
%Cluster the apex annotations from ecah marker so we have a concept of a
%single vessel apex (with multiple markers associated with it)
mkdir([rsa_dir 'apex_clusters\'])
for i_im = 1:num_images%1:
  
    %load in image
    im_num = rsa_images(i_im).name(1:end-4);   
    
    [vessels] = cluster_vessel_apices(im_num, [rsa_dir 'markup\'], markers, 20, 1);
    num_vessels = size(vessels.cluster_centres,1);
    vessels.majority_shapes = cell(num_vessels,1);
    vessels.majority_sizes = cell(num_vessels,1);
    vessels.num_markers = zeros(num_vessels,1);
    
    for i_v = 1:num_vessels;   
        
        [size_idx sizes] = grp2idx(vessels.cluster_sizes{i_v});
        [mode_idx] = mode(size_idx);
        vessels.majority_sizes{i_v} = sizes{mode_idx};
        
        [shape_idx shapes] = grp2idx(vessels.cluster_shapes{i_v});
        [mode_idx] = mode(shape_idx);
        vessels.majority_shapes{i_v} = shapes{mode_idx};
        
        vessels.num_markers(i_v) = length(vessels.cluster_sizes{i_v});
        
    end
    save([rsa_dir 'apex_clusters\' im_num '_apex_clusters.mat'], 'vessels');
end

%--------------------------------------------------------------------------
%% Now take a copy of the patch around each apex
normal_counter = 0;
giant_counter = 0;
enlarged_counter = 0;
apex_counter = 0;
vessel_counter = 0;
use_markers = [14];

create_folder([apex_dir 'normal']);
create_folder([apex_dir 'enlarged']);
create_folder([apex_dir 'giant']);

save_patch = 1;
save_apex_pts_txt = 1;
save_apex_pts_mat = 1;
save_vessel_pts_txt = 1;
save_vessel_pts_mat = 1;
save_vessel_patch = 1;
do_plot = 0;
num_v_pts = 15;
min_v_pts_sz = 30;

plot_rows = 3;
plot_cols = 5;
%
for i_im = 1:num_images
    
    %load in image
    im_num = rsa_images(i_im).name(1:end-4);   
    nailfold_name = [rsa_dir 'images\' im_num '.png'];   
    nailfold = imread(nailfold_name);
    [rows cols] = size(nailfold);
    
    %Loop through markers
    for i_ma = use_markers
        
        marker_name = markers(i_ma).name;
        %Check if marker has marked this image
        markup_list = dir([rsa_dir 'markup\' marker_name '\*#' im_num '*']);         
        if isempty(markup_list); continue; end
        
        %Read in vessel markup data - where multiple markup files exist,
        %use the last one
        markup = read_markup_from([rsa_dir 'markup\' marker_name '\' markup_list(end).name]);
        
        num_vessels = length(markup.vessels);
        
        for i_v = 1:num_vessels
            
            %Check this is a valid vessel
            anchor_xy = markup.vessels(i_v).anchor;
            
            if isempty(anchor_xy); continue; end                           
                
            %Check if distal
            is_distal = markup.vessels(i_v).ncm_vessel_properties.is_distal;
            
            if ~is_distal; continue; end
                
            %Find nearest marked apex to the anchor point
            num_apices = length(markup.vessels(i_v).apices);
            
            if ~num_apices; continue; end

            for i_a = 1:num_apices
                
                if isempty(markup.vessels(i_v).apices(i_a).inner_point)
                    continue;
                end
                vessel_shape = markup.vessels(i_v).ncm_vessel_properties.shape;
                vessel_size = markup.vessels(i_v).ncm_vessel_properties.size;

                switch vessel_size

                    case 'Normal'
                        %Increment the apex count
                        normal_counter = normal_counter + 1;
                        save_path = [apex_dir 'normal\apex' zerostr(normal_counter, 4)];
                        box_sz = 50;
                        v_step_sz = 2;

                        color = 'g';

                    case 'Enlarged'
                        %Increment the apex count
                        enlarged_counter = enlarged_counter + 1;
                        save_path = [apex_dir 'enlarged\apex' zerostr(enlarged_counter, 4)];
                        box_sz = 100;
                        v_step_sz = 8;
                        
                        color = 'm';

                    case 'Giant'
                        %Increment the apex count
                        giant_counter = giant_counter + 1;
                        save_path = [apex_dir 'giant\apex' zerostr(giant_counter, 4)];
                        box_sz = 200;
                        v_step_sz = 16;
                        
                        color = 'r';
                    otherwise
                        display(['Not recognised: ' markup.vessels(i_v).ncm_vessel_properties.size]);
                        continue;
                end
                apex_counter = apex_counter + 1;
                
                apex_properties.nailfold_name = nailfold_name;
                apex_properties.vessel_shape = vessel_shape;
                apex_properties.vessel_size = vessel_size;
                
                %Get start, end and centre of apex
                apex_xy = ...
                    [markup.vessels(i_v).apices.inner_point;...
                    markup.vessels(i_v).apices.outer_point];
                apex_cxy = mean(apex_xy);

                %Sample a patch from the nailfold about this apex
                sr = max(1, round(apex_cxy(:,2)) - box_sz);
                er = min(rows, round(apex_cxy(:,2)) + box_sz);
                sc = max(1, round(apex_cxy(:,1)) - box_sz);
                ec = min(cols, round(apex_cxy(:,1)) + box_sz);
                apex_patch = nailfold(sr:er, sc:ec);

                apex_xy(:,1) = apex_xy(:,1) - sc + 1;
                apex_xy(:,2) = apex_xy(:,2) - sr + 1;
                apex_cxy(:,1) = apex_cxy(:,1) - sc + 1;
                apex_cxy(:,2) = apex_cxy(:,2) - sr + 1;

                apex_properties.sr = sr;
                apex_properties.sc = sc;
                
                if save_patch
                    %Sample patch from 1st deriv.
                    imwrite(apex_patch, [save_path '.png']);
                end
                
                if save_apex_pts_txt
                    %Write out a pts file we can read in to VXL
                    fid1 = fopen([save_path '.pts'], 'wt');
                    fprintf(fid1, '%s \n', '{'); 
                    fprintf(fid1,'%.2f %.2f \n', apex_cxy(1,1), apex_cxy(1,2));
                    fprintf(fid1,'%.2f %.2f \n', apex_xy(1,1), apex_xy(1,2));
                    fprintf(fid1,'%.2f %.2f \n', apex_xy(2,1), apex_xy(2,2));
                    fprintf(fid1, '%s \n', '}');
                    fprintf(fid1, 'nailfold: %s \n', nailfold_name);
                    fprintf(fid1, 'marker: %s \n', marker_name);
                    fprintf(fid1, 'shape: %s \n', vessel_shape);
                    fprintf(fid1, 'size: %s \n', vessel_size);
                    fprintf(fid1, '%s %d \n', 'start_row:', sr);
                    fprintf(fid1, '%s %d \n', 'start_col: ', sc);
                    fclose(fid1);
                end
                
                if save_apex_pts_mat
                    %Wtie out points as a matlab file
                    
                    save([save_path '.mat'], 'apex_xy', 'apex_properties');
                end
                
                if isfield(markup.vessels(i_v), 'points') ...
                        && ~isempty(markup.vessels(i_v).points)
                    
                    %Get the marked vessel points and adjust to apex frame
                    v_pts = markup.vessels(i_v).points(:,1:2);                 
                    v_pts(:,1) = v_pts(:,1) - sc + 1;  
                    v_pts(:,2) = v_pts(:,2) - sr + 1;
                    
                    %We need a minimum number of vessel points or i's not
                    %worth bothering
                    if size(v_pts, 1) < min_v_pts_sz; continue; end
                        
                    %Check if the vessel path cross the apex, if not, skip
                    %this apex, otherwise use the crossing point as the
                    %apex of the vessel points
                    [do_cross v_apex v_top_idx] = line_cross(v_pts, apex_xy);
                    
                    if ~do_cross; continue; end
                    
                    %Increment the vessel counter
                    vessel_counter = vessel_counter + 1;
                    
                    %insert the crossing point into the vessel points (if
                    %it isn't already a member)
                    if ~any(ismember(v_pts, v_apex'))
                        v_pts = [v_pts(1:v_top_idx,:); v_apex'; v_pts(v_top_idx+1:end,:)];
                        v_top_idx = v_top_idx + 1;
                    end
                    
                    %Get points to the left of the apex
                    v_idx_l = num_v_pts*v_step_sz:-v_step_sz:v_step_sz;
                    vessel_l = v_pts(v_top_idx:-1:1,:);
                    dists = cumsum([0; sqrt(sum(diff(vessel_l).^2, 2))]);
                    if dists(end) < num_v_pts*v_step_sz
                        v_idx_l = linspace(0, dists(end), num_v_pts+1);
                        v_idx_l = v_idx_l(end:-1:2);
                    end
                    vessel_l = interp1(dists, vessel_l, v_idx_l, 'spline');
                    
                    %Get points to the right of the apex
                    v_idx_r = 0:v_step_sz:num_v_pts*v_step_sz;
                    vessel_r = v_pts(v_top_idx:end,:);
                    dists = cumsum([0; sqrt(sum(diff(vessel_r).^2, 2))]);
                    if dists(end) < num_v_pts*v_step_sz
                        v_idx_r = linspace(0, dists(end), num_v_pts+1);
                    end
                    vessel_r = interp1(dists, vessel_r, v_idx_r, 'spline');

                    %Complete the vessel with both parts
                    vessel_top = [vessel_l; vessel_r];

                    if save_vessel_pts_txt
                        %Write out a pts file we can read in to VXL
                        fid1 = fopen([save_path '_v_pts.pts'], 'wt');
                        fprintf(fid1, 'num_pts: %d \n', 2*num_v_pts+1);
                        fprintf(fid1, '%s \n', '{');
                        for i_vp = 1:2*num_v_pts+1
                            fprintf(fid1,'%.2f %.2f \n', vessel_top(i_vp,1), vessel_top(i_vp,2));
                        end
                        fprintf(fid1, '%s \n', '}');
                        fprintf(fid1, 'nailfold: %s \n', nailfold_name);
                        fprintf(fid1, 'marker: %s \n', marker_name);
                        fprintf(fid1, 'shape: %s \n', vessel_shape);
                        fprintf(fid1, 'size: %s \n', vessel_size);
                        fprintf(fid1, '%s %d \n', 'start_row:', sr);
                        fprintf(fid1, '%s %d \n', 'start_col: ', sc);
                        fclose(fid1);
                        
                    end
                    if save_vessel_pts_mat
                        save([save_path '_v_pts.mat'], 'vessel_top', 'apex_properties');
                    end
                    
                    if save_vessel_patch
                        %Get the marked vessel points
                        v_pts = markup.vessels(i_v).points(:,1:2);
                        
                        sr = max(1, round(min(v_pts(:,2))) - box_sz);
                        er = min(rows, round(max(v_pts(:,2))) + box_sz);
                        sc = max(1, round(min(v_pts(:,1))) - box_sz);
                        ec = min(cols, round(max(v_pts(:,1))) + box_sz);
                        vessel_patch = nailfold(sr:er, sc:ec);
                        
                        v_pts(:,1) = v_pts(:,1) - sc + 1;  
                        v_pts(:,2) = v_pts(:,2) - sr + 1;
                        
                        vessel_properties = apex_properties;
                        vessel_properties.sr = sr;
                        vessel_properties.sc = sc;
                        save([save_path '_vessel.mat'], 'vessel_patch', 'v_pts', 'vessel_properties');
                    end                       
                    
                    if do_plot

                        plot_num = mod(vessel_counter-1, plot_rows*plot_cols)+1;
                        if plot_num == 1;
                            figure;
                        end

                        subplot(plot_rows,plot_cols,plot_num); imgray(apex_patch);
                        plot(apex_xy(:,1), apex_xy(:,2), [color '-x']);

                        plot(v_pts(:,1),v_pts(:,2), 'g-');
                        plot(vessel_top(:,1), vessel_top(:,2), 'c');
                        plot(vessel_top(1,1), vessel_top(1,2), 'ro');
                        plot(vessel_top(num_v_pts+1,1), vessel_top(num_v_pts+1,2), 'bo');
                        plot(vessel_top(2*num_v_pts+1,1), vessel_top(2*num_v_pts+1,2), 'yo');
                        
                        title(vessel_shape);
                    end
                end
            end
        end
        
    end
end
%--------------------------------------------------------------------------
%%
for i_ap = 361:540    
    
    save_path = [apex_dir 'normal\apex' zerostr(i_ap, 4)];
    apex_patch = imread([save_path '.png']);
    load([save_path '.mat'], 'apex_xy');
    
    plot_num = mod(i_ap-1, 18)+1;
    if plot_num == 1;
        figure;
    end
    subplot(3,6,plot_num); imgray(apex_patch);
    plot(apex_xy(:,1), apex_xy(:,2), 'g-x');
end
%%
template_sz = 49;
scale_factors = [1 2 4];
vessel_sizes = {'normal','enlarged', 'giant'};%

mean_apex_i = zeros(template_sz, template_sz, 3, 3);
mean_apex_g1 = zeros(template_sz, template_sz, 3, 3);
mean_apex_g2 = zeros(template_sz, template_sz, 3, 3);
shape_counts = zeros(3,3);

for i_sz = 1:3
    vessel_size = vessel_sizes{i_sz};
    scale_factor = scale_factors(i_sz);
    
    apex_list = dir([apex_dir vessel_size '\apex*.png']);
    num_apexes = length(apex_list);    



    for i_ap = 1:num_apexes
        save_path = [apex_dir vessel_size '\apex' zerostr(i_ap, 4)];
        apex_patch = double(imread([save_path '.png']));

        [mag_1d] = gaussian_1st_derivative_gradient(apex_patch, 2*scale_factor);
        [mag_2d] = gaussian_2nd_derivative_line(apex_patch, 2*scale_factor);

        load([save_path '.mat'], 'apex_xy', 'apex_properties');
        apex_cxy = mean(apex_xy);

        switch apex_properties.vessel_shape
            case 'Normal'
                i_sh = 1;
            case 'Meandering'
                i_sh = 2;
            case 'Angiogenetic'
                i_sh = 3;
            otherwise
                display(['Shape: ' apex_properties.vessel_shape ' not recognised']);
        end
        shape_counts(i_sh,i_sz) = shape_counts(i_sh,i_sz) + 1;
        
        intensity_patch = sample_window(apex_patch,...
            (template_sz-1)*scale_factor+1, round(apex_cxy(:,2)), round(apex_cxy(:,1)));
        intensity_patch(isnan(intensity_patch)) = mean(intensity_patch(~isnan(intensity_patch)));
        
        g1_patch = sample_window(mag_1d, ...
            (template_sz-1)*scale_factor+1, round(apex_cxy(:,2)), round(apex_cxy(:,1)), 0);
        g2_patch = sample_window(mag_2d, ...
            (template_sz-1)*scale_factor+1, round(apex_cxy(:,2)), round(apex_cxy(:,1)), 0);

        mean_apex_i(:,:,i_sh,i_sz) = mean_apex_i(:,:,i_sh,i_sz) +...
            imresize(intensity_patch, [template_sz, template_sz]);
        mean_apex_g1(:,:,i_sh,i_sz) = mean_apex_g1(:,:,i_sh,i_sz) + ...
            imresize(g1_patch, [template_sz, template_sz]);
        mean_apex_g2(:,:,i_sh,i_sz) = mean_apex_g2(:,:,i_sh,i_sz) + ...
            imresize(g2_patch, [template_sz, template_sz]);
    end
    %
    figure;
    for i_sh = 1:3
        subplot(3,4,i_sh + 0); imgray(mean_apex_i(:,:,i_sh,i_sz));
        subplot(3,4,i_sh + 4); imgray(mean_apex_g1(:,:,i_sh,i_sz));
        subplot(3,4,i_sh + 8); imgray(mean_apex_g2(:,:,i_sh,i_sz));

    end
    subplot(3,4,4); imgray(sum(mean_apex_i(:,:,:,i_sz),3));
    subplot(3,4,8); imgray(sum(mean_apex_g1(:,:,:,i_sz),3));
    subplot(3,4,12); imgray(sum(mean_apex_g2(:,:,:,i_sz),3));

    subplot(3,4,1); ylabel('Image intensities');
    subplot(3,4,5); ylabel('1st Derivs');
    subplot(3,4,9); ylabel('2nd Derivs');

    subplot(3,4,1); title('Normal');
    subplot(3,4,2); title('Meandering');
    subplot(3,4,3); title('Angiogenic');
    subplot(3,4,4); title('All');

    subplot(3,4,9); xlabel([num2str(shape_counts(1,i_sz)) ' apexes']);
    subplot(3,4,10); xlabel([num2str(shape_counts(2,i_sz)) ' apexes']);
    subplot(3,4,11); xlabel([num2str(shape_counts(3,i_sz)) ' apexes']);
    subplot(3,4,12); xlabel([num2str(sum(shape_counts(:,i_sz))) ' apexes']);
end
%%
template_sz = 49;
scale_factors = [1 2 4];
vessel_sizes = {'normal','enlarged', 'giant'};%

apex_shapes = cell(3,3);
apex_areas = cell(3,3);
for i_sz = 1:3
    for i_sh = 1:3
        apex_shapes{i_sh,i_sz} = zeros(shape_counts(i_sh,i_sz), 62);
        apex_areas{i_sh,i_sz} = zeros(shape_counts(i_sh,i_sz),1);
    end
end
shape_counts_v = zeros(3,3);

for i_sz = 1:3
    vessel_size = vessel_sizes{i_sz};
    scale_factor = scale_factors(i_sz);
    
    apex_list = dir([apex_dir vessel_size '\apex*v_pts.mat']);
    num_apexes = length(apex_list);

    for i_ap = 1:num_apexes
        
        load([apex_dir vessel_size '\' apex_list(i_ap).name], 'vessel_top', 'apex_properties');

        switch apex_properties.vessel_shape
            case 'Normal'
                i_sh = 1;
            case 'Meandering'
                i_sh = 2;
            case 'Angiogenetic'
                i_sh = 3;
            otherwise
                display(['Shape: ' apex_properties.vessel_shape ' not recognised']);
        end
        shape_counts_v(i_sh,i_sz) = shape_counts_v(i_sh,i_sz) + 1;
        
        apex_shapes{i_sh,i_sz}(shape_counts_v(i_sh,i_sz),:) = vessel_top(:);
        apex_areas{i_sh,i_sz}(shape_counts_v(i_sh,i_sz)) = polyarea(vessel_top(:,1), vessel_top(:,2));
    end
    
    figure;
    for i_sh = 1:3
        subplot(1,3,i_sh); hold all;
        for i_ap = 1:shape_counts_v(i_sh,i_sz)
            plot(apex_shapes{i_sh,i_sz}(i_ap,1:31), apex_shapes{i_sh,i_sz}(i_ap,32:62));
        end
    end

    subplot(1,3,1); title('Normal');
    subplot(1,3,2); title('Meandering');
    subplot(1,3,3); title('Angiogenic');

end
%%
apex_rots = cell(3,3);
apex_trans = cell(3,3);
apex_scales = cell(3,3);
apex_theta = cell(3,3);
%%
for i_sz = 2:3
    for i_sh = 1:3
        
        [a_shapes, a_scales, mean_target, a_rots, a_trans a_origins] = align_shapes(apex_shapes{i_sh,i_sz}, 'area', mean(apex_areas{i_sh,i_sz}));
        mean_shape = [mean(a_shapes(:,1:31))', mean(a_shapes(:,32:62))'];

        theta = -mean(acos(squeeze(a_rots(1,1,:))));
        theta_rot = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        apex_rots{i_sh,i_sz} = a_rots;
        apex_trans{i_sh,i_sz} = a_trans;
        apex_scales{i_sh,i_sz} = a_scales;
        apex_theta{i_sh,i_sz} = theta_rot;

        figure; hold on;
        plot(a_shapes(:,1:31)', a_shapes(:,32:62)'); axis ij equal;
        plot(mean_shape(:,1), mean_shape(:,2), 'k', 'linewidth', 2)

        mean_shape = mean_target * theta_rot';
        mean_shape = [mean_shape(:,1) - mean_shape(16,1) mean_shape(:,2)-mean_shape(16,2)];
        figure; plot(mean_shape(:,1), mean_shape(:,2)); axis ij equal;
    end
end
%%
template_sz = 49;
scale_factors = [1 2 4];
vessel_sizes = {'normal','enlarged', 'giant'};%

mean_apex_rot_i = zeros(template_sz, template_sz, 3, 3);
mean_apex_rot_g1 = zeros(template_sz, template_sz, 3, 3);
mean_apex_rot_g2 = zeros(template_sz, template_sz, 3, 3);
shape_counts_v = zeros(3,3);    

for i_sz = 1:3
    vessel_size = vessel_sizes{i_sz};
    scale_factor = scale_factors(i_sz);
    
    lim = (template_sz-1)*scale_factor/2;
    x = repmat(-lim:scale_factor:lim, 49, 1);
    y = x';
    xy = [x(:) y(:)];
    
    apex_list = dir([apex_dir vessel_size '\apex*v_pts.mat']);
    num_apexes = length(apex_list);

    for i_ap = 1:num_apexes
        
        load([apex_dir vessel_size '\' apex_list(i_ap).name], 'vessel_top', 'apex_properties');
        apex_num = apex_list(i_ap).name(5:8);       
        apex_patch = double(imread([apex_dir vessel_size '\apex' apex_num '.png']));

        [mag_1d] = gaussian_1st_derivative_gradient(apex_patch, 2*scale_factor);
        [mag_2d] = gaussian_2nd_derivative_line(apex_patch, 2*scale_factor);

        switch apex_properties.vessel_shape
            case 'Normal'
                i_sh = 1;
            case 'Meandering'
                i_sh = 2;
            case 'Angiogenetic'
                i_sh = 3;
            otherwise
                display(['Shape: ' apex_properties.vessel_shape ' not recognised']);
        end
        
        shape_counts_v(i_sh,i_sz) = shape_counts_v(i_sh,i_sz) + 1;
        i_ve = shape_counts_v(i_sh,i_sz);
        
        xya = (xy * apex_rots{i_sh,i_sz}(:,:,i_ve)' / apex_scales{i_sh,i_sz}(i_ve))* apex_theta{i_sh,i_sz};
        xa = reshape(xya(:,1) + vessel_top(16,1), template_sz, template_sz);
        ya = reshape(xya(:,2) + vessel_top(16,2), template_sz, template_sz);

        intensity_patch = interp2(apex_patch, xa, ya, 'bilinear');
        intensity_patch(isnan(intensity_patch)) = mean(intensity_patch(~isnan(intensity_patch)));
        
        mean_apex_rot_i(:,:,i_sh,i_sz) = mean_apex_rot_i(:,:,i_sh,i_sz) +...
            intensity_patch;
        mean_apex_rot_g1(:,:,i_sh,i_sz) = mean_apex_rot_g1(:,:,i_sh,i_sz) +...
            interp2(mag_1d, xa, ya, 'bilinear', 0);
        mean_apex_rot_g2(:,:,i_sh,i_sz) = mean_apex_rot_g2(:,:,i_sh,i_sz) +...
            interp2(mag_2d, xa, ya, 'bilinear', 0);

    end
    %
    figure;
    for i_sh = 1:3
        subplot(3,4,i_sh + 0); imgray(mean_apex_rot_i(:,:,i_sh,i_sz));
        subplot(3,4,i_sh + 4); imgray(mean_apex_rot_g1(:,:,i_sh,i_sz));
        subplot(3,4,i_sh + 8); imgray(mean_apex_rot_g2(:,:,i_sh,i_sz));

    end
    subplot(3,4,4); imgray(sum(mean_apex_rot_i(:,:,:,i_sz),3));
    subplot(3,4,8); imgray(sum(mean_apex_rot_g1(:,:,:,i_sz),3));
    subplot(3,4,12); imgray(sum(mean_apex_rot_g2(:,:,:,i_sz),3));

    subplot(3,4,1); ylabel('Image intensities');
    subplot(3,4,5); ylabel('1st Derivs');
    subplot(3,4,9); ylabel('2nd Derivs');

    subplot(3,4,1); title('Normal');
    subplot(3,4,2); title('Meandering');
    subplot(3,4,3); title('Angiogenic');
    subplot(3,4,4); title('All');

    subplot(3,4,9); xlabel([num2str(shape_counts(1,i_sz)) ' apexes']);
    subplot(3,4,10); xlabel([num2str(shape_counts(2,i_sz)) ' apexes']);
    subplot(3,4,11); xlabel([num2str(shape_counts(3,i_sz)) ' apexes']);
    subplot(3,4,12); xlabel([num2str(sum(shape_counts(:,i_sz))) ' apexes']);
end
%%
mkdir([rsa_dir '\models']);
save([rsa_dir '\models\apex_templates.mat'], 'mean_apex_*');
%%
load([rsa_dir '\models\apex_templates.mat'], 'mean_apex_*');

mkdir([rsa_dir '\results\template_matching\original\intensity']);
mkdir([rsa_dir '\results\template_matching\aligned\intensity']); 
mkdir([rsa_dir '\results\template_matching\original\g12d']);
mkdir([rsa_dir '\results\template_matching\aligned\g12d']);

for i_im = 2:num_images
    
    %load in image
    im_num = rsa_images(i_im).name(1:end-4);   
    im = imread([rsa_dir 'images\' im_num '.png']);
    
    for i_sz = 1:3
        
        for i_sh = 1:3
            
            %Templates for Gaussian derivatives - original and aligned
            templates = cell(2,2);
            templates{1,1} = 'g1d';
            templates{2,1} = 'g2d';
            
            %Original
            templates{1,2} = mean_apex_g1(:,:,i_sh,i_sz);
            templates{2,2} = mean_apex_g2(:,:,i_sh,i_sz);            
            [maxima_pos, maxima_vals] = template_match_apexes(im, templates, 'plot', 0);

            save([rsa_dir '\results\template_matching\original\g12d\' im_num ...
                '_apex_candidates' num2str(i_sh) '_' num2str(i_sz) '.mat'],...
                'maxima_pos', 'maxima_vals');
            
            %Aligned
            templates{1,2} = mean_apex_rot_g1(:,:,i_sh,i_sz);
            templates{2,2} = mean_apex_rot_g2(:,:,i_sh,i_sz);           
            [maxima_pos, maxima_vals] = template_match_apexes(im, templates, 'plot', 0);
            
            save([rsa_dir '\results\template_matching\aligned\g12d\' im_num ...
                '_apex_candidates' num2str(i_sh) '_' num2str(i_sz) '.mat'],...
                'maxima_pos', 'maxima_vals');
            
            %Templates for intensity
            templates = cell(1,1);
            templates{1,1} = 'intensity';
            
            %Original
            templates{1,2} = mean_apex_i(:,:,i_sh,i_sz);            
            [maxima_pos, maxima_vals] = template_match_apexes(im, templates, 'plot', 0);

            save([rsa_dir '\results\template_matching\original\intensity\' im_num ...
                '_apex_candidates' num2str(i_sh) '_' num2str(i_sz) '.mat'],...
                'maxima_pos', 'maxima_vals');

            %aligned
            templates{1,2} = mean_apex_rot_i(:,:,i_sh,i_sz);
            [maxima_pos, maxima_vals] = template_match_apexes(im, templates, 'plot', 0);

            save([rsa_dir '\results\template_matching\aligned\intensity\' im_num ...
                '_apex_candidates' num2str(i_sh) '_' num2str(i_sz) '.mat'],...
                'maxima_pos', 'maxima_vals');
                   
        end
    end
end
%%
for i_im = 1:2%num_images
    
    %load in image
    im_num = rsa_images(i_im).name(1:end-4);   
    im = imread([rsa_dir 'images\' im_num '.png']);
    
    figure; imgray(im);
    for i_sz = 1:3
        for i_sh = 1:3
            
            load([rsa_dir '\results\template_matching\original\' im_num ...
                '_apex_candidates' num2str(i_sh) '_' num2str(i_sz) '.mat'],...
                'maxima_pos', 'maxima_vals');
            
            plot(maxima_pos(1:20, 1), maxima_pos(1:20, 2),...
                'x', 'markersize', 3*i_sz);
            
            load([rsa_dir '\results\template_matching\aligned\' im_num ...
                '_apex_candidates' num2str(i_sh) '_' num2str(i_sz) '.mat'],...
                'maxima_pos', 'maxima_vals');
            plot(maxima_pos(1:20, 1), maxima_pos(1:20, 2),...
                'o', 'markersize', 3*i_sz);
        end
    end
end

    
    
    

    
    

    


            
            
        



