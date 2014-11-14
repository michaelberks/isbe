%Get list of images
pilot_dir = 'C:\isbe\nailfold\data\pilot_study\';
pilot_images = dir([pilot_dir 'images\*.png']);
num_images = length(pilot_images);

%Get list of markers
markers = dir([pilot_dir '\markup\']);
markers = markers(3:end);

num_markers = length(markers);
%%
%Loop through images
for i_im = 1:num_images
    
    %load in image
    im_num = pilot_images(i_im).name(1:end-4);   
    im = imread([pilot_dir 'images\' im_num '.png']);
    
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
        markup_list = dir([pilot_dir 'markup\' markers(i_ma).name '\*#' im_num '*']);         
        if isempty(markup_list); continue; end
        
        %Read in vessel markup data - where multiple markup files exist,
        %use the last one
        markup = read_markup_from([pilot_dir 'markup\' markers(i_ma).name '\' markup_list(end).name]);
        
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
%%
clc
for i_m = 1:num_markers
    display(['Marker ' num2str(i_m) ': ' markers(i_m).name]);
end
%%
all_vessel_marker_counts = zeros(0,1);
all_vessel_sizes = cell(0,1);
all_vessel_shapes = cell(0,1);

all_vessel_marker_counts_d = zeros(0,1);
all_vessel_sizes_d = cell(0,1);
all_vessel_shapes_d = cell(0,1);
    
all_vessel_marker_counts_e = zeros(0,1);
all_vessel_sizes_e = cell(0,1);
all_vessel_shapes_e = cell(0,1);
%
%Loop through images
for i_im = 1:num_images%1:
    
    %load in image
    im_num = pilot_images(i_im).name(1:end-4);   
    im = imread([pilot_dir 'images\' im_num '.png']);
    
    %get mask of image to set contrast
    mask = im < 250;
    mask = imerode(mask, strel('disk', 10));
    g_max = max(im(mask));
    g_min = min(im(mask));  
    
    %Display image
    figure; 
    subplot(2,1,1); imgray(im); caxis([g_min g_max]); a1 = gca;
    subplot(2,1,2); imgray(im); caxis([g_min g_max]); a2 = gca;
    
    %------------------------------------------------------------------
    [vessels] = cluster_vessel_apices(im_num, [pilot_dir 'markup\'], markers, 20);
    for i_v = 1:size(vessels.cluster_centres,1);
        
        
        all_vessel_marker_counts(end+1,:) = length(vessels.cluster_shapes{i_v}); %#ok
        
        [shape_idx shapes] = grp2idx(vessels.cluster_shapes{i_v});
        
        if length(shapes) >= 3
            display(['Im ' num2str(i_im) ', vessel ' num2str(i_v) ', ' num2str(length(shapes)) ' shapes.']);
            plot(vessels.cluster_xy{i_v}(:,1), vessels.cluster_xy{i_v}(:,2), 'x');
        end
        
        [mode_idx, ~, all_modes] = mode(shape_idx);
        if length(all_modes) > 1;
            display('we have a tie for shapes!');
            display(vessels.cluster_shaps{i_v});
        end
        
        switch shapes{mode_idx}
            case 'Normal'    
                m_color = 'r';
            case 'Meandering'
                m_color = 'g';
            case 'Angiogenetic'
                m_color = 'b';
            case 'Undefined'
                m_color = 'k';
            otherwise
                display(['Shape: ' shapes{shape_idx} ' not recognised']);
        end
        all_vessel_shapes{end+1,:} = shapes{mode_idx}; %#ok
        
        [size_idx sizes] = grp2idx(vessels.cluster_sizes{i_v});
        
        if length(sizes) >= 3
            display(['Im ' num2str(i_im) ', vessel ' num2str(i_v) ', ' num2str(length(sizes)) ' sizes.']);
            plot(vessels.cluster_xy{i_v}(:,1), vessels.cluster_xy{i_v}(:,2), 'x');
        end
        
        [mode_idx, ~, all_modes] = mode(size_idx);
        if length(all_modes) > 1;
            display('we have a tie for sizes!');
            display(vessels.cluster_sizes{i_v});
        end
        
        switch sizes{mode_idx}

            case 'Normal'
                m_type = '+';
            case 'Enlarged'    
                m_type = 'o';
            case 'Giant'    
                m_type = '^';
            case 'Undefined'
                m_type = '*';
            otherwise
                display(['Size: ' sizes{shape_idx} ' not recognised']);
        end
        all_vessel_sizes{end+1,:} = sizes{mode_idx}; %#ok
        plot(a1, vessels.cluster_centres(i_v,1), vessels.cluster_centres(i_v,2), [m_color m_type], 'markersize', 10);
    end

    %------------------------------------------------------------------
    [vessels_d] = cluster_vessel_apices(im_num, [pilot_dir 'markup\'], markers, 20, 1);    
    for i_v = 1:size(vessels_d.cluster_centres,1);   
        all_vessel_marker_counts_d(end+1,:) = length(vessels_d.cluster_shapes{i_v}); %#ok
        
        [shape_idx shapes] = grp2idx(vessels_d.cluster_shapes{i_v});
        [mode_idx, ~, all_modes] = mode(shape_idx);
        if length(all_modes) > 1;
            display('we have a tie for shapes!');
            display(vessels_d.cluster_shaps{i_v});
        end
        
        switch shapes{mode_idx}
            case 'Normal'    
                m_color = 'r';
            case 'Meandering'
                m_color = 'g';
            case 'Angiogenetic'
                m_color = 'b';
            case 'Distal'
                m_color = 'c';
            case 'Undefined'
                m_color = 'k';
            otherwise
                display(['Shape: ' shapes{shape_idx} ' not recognised']);
        end
        all_vessel_shapes_d{end+1,:} = shapes{mode_idx}; %#ok
        
        [size_idx sizes] = grp2idx(vessels_d.cluster_sizes{i_v});
        [mode_idx, ~, all_modes] = mode(size_idx);
        if length(all_modes) > 1;
            display('we have a tie for sizes!');
            display(vessels_d.cluster_sizes{i_v});
        end
        
        switch sizes{mode_idx}

            case 'Normal'
                m_type = '+';
            case 'Enlarged'    
                m_type = 'o';
            case 'Giant'    
                m_type = '^';
            case 'Distal'
                m_type = 'v';
            case 'Undefined'
                m_type = '*';
            otherwise
                display(['Size: ' sizes{shape_idx} ' not recognised']);
        end
        all_vessel_sizes_d{end+1,:} = sizes{mode_idx}; %#ok
        
        %plot(vessels.cluster_xy{i_v}(:,1), vessels.cluster_xy{i_v}(:,2), 'x');
        plot(a2, vessels_d.cluster_centres(i_v,1), vessels_d.cluster_centres(i_v,2), [m_color m_type], 'markersize', 10);
        
        %------------------------------------------------------------------
        distal_idx = strcmp(vessels_d.cluster_shapes{i_v}, 'Distal');
        vessels_d.cluster_shapes{i_v}(distal_idx) = [];
        vessels_d.cluster_sizes{i_v}(distal_idx) = [];
        
        if isempty(vessels_d.cluster_shapes{i_v}); continue; end
        
        all_vessel_marker_counts_e(end+1,:) = length(distal_idx); %#ok
        
        [shape_idx shapes] = grp2idx(vessels_d.cluster_shapes{i_v});
        [mode_idx, ~, all_modes] = mode(shape_idx);
        if length(all_modes) > 1;
            display('we have a tie for shapes!');
            display(vessels_d.cluster_shaps{i_v});
        end
        
        switch shapes{mode_idx}
            case 'Normal'    
                m_color = 'r';
            case 'Meandering'
                m_color = 'g';
            case 'Angiogenetic'
                m_color = 'b';
            case 'Distal'
                m_color = 'c';
            case 'Undefined'
                m_color = 'k';
            otherwise
                display(['Shape: ' shapes{shape_idx} ' not recognised']);
        end
        all_vessel_shapes_e{end+1,:} = shapes{mode_idx}; %#ok
        
        [size_idx sizes] = grp2idx(vessels_d.cluster_sizes{i_v});
        [mode_idx, ~, all_modes] = mode(size_idx);
        if length(all_modes) > 1;
            display('we have a tie for sizes!');
            display(vessels_d.cluster_sizes{i_v});
        end
        
        switch sizes{mode_idx}

            case 'Normal'
                m_type = '+';
            case 'Enlarged'    
                m_type = 'o';
            case 'Giant'    
                m_type = '^';
            case 'Distal'
                m_type = 'v';
            case 'Undefined'
                m_type = '*';
            otherwise
                display(['Size: ' sizes{shape_idx} ' not recognised']);
        end
        all_vessel_sizes_e{end+1,:} = sizes{mode_idx}; %#ok
        
    end
end
%%
figure; 
hist(all_vessel_marker_counts, 1:10);
title('All vessels')
xlabel('Number of markers that selected the vessel');
%%
figure;
[shape_idx shapes] = grp2idx(all_vessel_shapes);
[size_idx sizes] = grp2idx(all_vessel_sizes);
include_shapes = [1 2 4];
include_sizes = [3 2 1];

for i_sh = 1:3%length(shapes)
    for i_sz = 1:3%length(sizes)
        subplot(3,3,(i_sh-1)*3 + i_sz); 
        hist(all_vessel_marker_counts((shape_idx == include_shapes(i_sh)) & (size_idx == include_sizes(i_sz))), 1:10);
        title({[shapes{include_shapes(i_sh)} ' shape, '  sizes{include_sizes(i_sz)} ' size vessels'];});
        ylabel(['Total: ' num2str(sum((shape_idx == include_shapes(i_sh)) & (size_idx == include_sizes(i_sz)))) ' apices']);
        set(gca, 'xlim', [.5 10.5]);
    end
end
%%
figure;
[shape_idx shapes] = grp2idx(all_vessel_shapes_d);
[size_idx sizes] = grp2idx(all_vessel_sizes_d);
include_shapes = [2 3 5];
include_sizes = [4 3 2];

for i_sh = 1:3%length(shapes)
    for i_sz = 1:3%length(sizes)
        subplot(3,3,(i_sh-1)*3 + i_sz); 
        hist(all_vessel_marker_counts_d((shape_idx == include_shapes(i_sh)) & (size_idx == include_sizes(i_sz))), 1:10);
        title({[shapes{include_shapes(i_sh)} ' shape, '  sizes{include_sizes(i_sz)} ' size vessels'];});
        ylabel(['Total: ' num2str(sum((shape_idx == include_shapes(i_sh)) & (size_idx == include_sizes(i_sz)))) ' apices']);
        set(gca, 'xlim', [.5 10.5]);
    end
end
%%
figure;
[shape_idx shapes] = grp2idx(all_vessel_shapes_e);
[size_idx sizes] = grp2idx(all_vessel_sizes_e);
include_shapes = [1 2 4];
include_sizes = [3 2 1];

for i_sh = 1:3%length(shapes)
    for i_sz = 1:3%length(sizes)
        subplot(3,3,(i_sh-1)*3 + i_sz); 
        hist(all_vessel_marker_counts_e((shape_idx == include_shapes(i_sh)) & (size_idx == include_sizes(i_sz))), 1:10);
        title({[shapes{include_shapes(i_sh)} ' shape, '  sizes{include_sizes(i_sz)} ' size vessels'];});
        ylabel(['Total: ' num2str(sum((shape_idx == include_shapes(i_sh)) & (size_idx == include_sizes(i_sz)))) ' apices']);
        set(gca, 'xlim', [.5 10.5]);
    end
end
%%
shapes = {'Distal','Normal','Meandering','Undefined','Angiogenetic'};
sizes = {'Distal','Giant','Enlarged','Normal','Undefined'};

vessel_shape_mixes = zeros(0,5);
vessel_size_mixes = zeros(0,5);
vessel_marker_counter = zeros(0,1);
for i_im = 1:num_images%1:
  
    %load in image
    im_num = pilot_images(i_im).name(1:end-4);   
    
    vessels_d = u_load([pilot_dir 'apex_clusters\' im_num '_apex_clusters.mat']);
    for i_v = 1:size(vessels_d.cluster_centres,1);   
        
        vessel_size_mixes(end+1,:) = zeros(1,5); %#ok
        vessel_shape_mixes(end+1,:) = zeros(1,5); %#ok
        vessel_marker_counter(end+1,:) = length(vessels_d.cluster_sizes{i_v}); %#ok
        for i_sz = 1:5
            vessel_size_mixes(end,i_sz) = ...
                any(strcmp(vessels_d.cluster_sizes{i_v}, sizes{i_sz}));
            vessel_shape_mixes(end,i_sz) = ...
                any(strcmp(vessels_d.cluster_shapes{i_v}, shapes{i_sz}));
        end
        
    end
end

%%
for i_sz1 = 1:5
    
    for i_sz2 = i_sz1:5
        if i_sz2==i_sz1
            num_vessels = sum(vessel_size_mixes(:,i_sz1) == sum(vessel_size_mixes,2));
        else
            num_vessels = sum(vessel_size_mixes(:,i_sz1) & vessel_size_mixes(:,i_sz2));
        end
        display([sizes{i_sz1} ' with ' sizes{i_sz2} ': ' num2str(num_vessels) '  vessels']);
    end
end
%%
for i_sh1 = 1:5
    
    for i_sh2 = i_sh1:5
        if i_sh2==i_sh1
            num_vessels = sum(vessel_shape_mixes(:,i_sh1) == sum(vessel_shape_mixes,2));
        else
            num_vessels = sum(vessel_shape_mixes(:,i_sh1) & vessel_shape_mixes(:,i_sh2));
        end
        display([shapes{i_sh1} ' with ' shapes{i_sh2} ': ' num2str(num_vessels) '  vessels']);
    end
end
%%
vessel_size_mixes2 = vessel_size_mixes;
vessel_shape_mixes2 = vessel_shape_mixes;
vessel_marker_counter2 = vessel_marker_counter;

distal_idx = sum(vessel_shape_mixes2(:,1) == sum(vessel_shape_mixes2,2));
vessel_size_mixes2(distal_idx,:) = [];
vessel_shape_mixes2(distal_idx,:) = [];
vessel_marker_counter2(distal_idx,:) = [];

for i_m = 1:10
    vessel_size_mixes3 = vessel_size_mixes2;
    vessel_shape_mixes3 = vessel_shape_mixes2;
    vessel_marker_counter3 = vessel_marker_counter2;

    at_least_idx = vessel_marker_counter2 ~= i_m;
    vessel_size_mixes3(at_least_idx,:) = [];
    vessel_shape_mixes3(at_least_idx,:) = [];
    vessel_marker_counter3(at_least_idx,:) = [];

    vessel_size_mixes3(:,1) = [];
    vessel_shape_mixes3(:,1) = [];

    shapes_purity = sum(vessel_shape_mixes3,2);
    sizes_purity = sum(vessel_size_mixes3,2);

    display(['Number of markers ' num2str(i_m)]);
    for ii = 1:4
        display(['Shapes with ' num2str(ii) ' label(s): ' num2str(sum(shapes_purity==ii))]);
        display(['Sizes with ' num2str(ii) ' label(s): ' num2str(sum(sizes_purity==ii))]);
    end
end
%%
%Loop through images
colors = lines(10);
for i_im = 1:num_images%1:
    
    %load in image
    im_num = pilot_images(i_im).name(1:end-4);   
    im = imread([pilot_dir 'images\' im_num '.png']);
    [rows cols] = size(im); 
    
    [vessels] = cluster_vessel_apices(im_num, [pilot_dir 'markup\'], markers, 20, 1);    
    for i_v = 1:size(vessels.cluster_centres,1);
        
        distal_idx = strcmp(vessels.cluster_shapes{i_v}, 'Distal');
      
        if all(distal_idx); continue; end
        
        [size_idx sizes] = grp2idx(vessels.cluster_sizes{i_v});
        [mode_idx] = mode(size_idx);
        
        switch sizes{mode_idx}

            case 'Normal'
                box_sz = 50;
            case 'Enlarged'    
                box_sz = 100;
            case 'Giant'    
                box_sz = 200;
            case 'Distal'
                box_sz = 100;
            case 'Undefined'
                box_sz = 100;
            otherwise
                display(['Size: ' sizes{shape_idx} ' not recognised']);
        end
        
        sr = max(1, round(vessels.cluster_centres(i_v,2)) - box_sz);
        er = min(rows, round(vessels.cluster_centres(i_v,2)) + box_sz);
        sc = max(1, round(vessels.cluster_centres(i_v,1)) - box_sz);
        ec = min(cols, round(vessels.cluster_centres(i_v,1)) + box_sz);
        apex_patch = im(sr:er, sc:ec);
        

        f1 = figure(...
            'windowstyle', 'normal',...
            'Units', 'pixels',...
            'position', [0 0 1300 900],...
            'PaperPositionMode','auto',...
            'Visible', 'off');
        a1 = axes(...
            'Units', 'pixels',...
            'position', [0 0 900 900]); 
        imgray(apex_patch);       
        
        a2 = axes(...
            'Units', 'pixels',...
            'position', [900 0 400 900],...
            'Xlim', [0 400],...
            'Ylim', [0 900]); 
        axis off; hold on;
        
        for i_ma = 1:10
            idx = find(vessels.cluster_members{i_v} == i_ma);
            if ~isempty(idx)
                
                axes(a2);
                text(10, 950 - i_ma*90, ['Marker ' num2str(i_ma) ': '...
                    vessels.cluster_shapes{i_v}{idx} ', ' vessels.cluster_sizes{i_v}{idx}],...
                    'FontSize', 16, 'Color', colors(i_ma,:));
                
                plot(a1, vessels.cluster_xy{i_v}(idx,1) - sc, vessels.cluster_xy{i_v}(idx,2) - sr,...
                    'x', 'MarkerSize', 20, 'MarkerEdgeColor', colors(i_ma,:));
        
            else
                axes(a2);
                text(10, 950 - i_ma*90, ['Marker ' num2str(i_ma) ': '],...
                    'FontSize', 16, 'Color', [0.5 0.5 0.5]);
            end
        end
        fname = ['C:\isbe\nailfold\data\temp\vessel_display_images\vessel' zerostr(i_im,2) '_' zerostr(i_v,3)];
        
        print('-depsc2', [fname '.eps'], f1);
        eval(['!convert -density 150 "', fname,'.eps" "',fname,'.png"']);
        delete([fname '.eps']);
        close(f1);
        
        %------------------------------------------------------------------
        
    end
end
%%
mkdir([pilot_dir 'apex_clusters\display_images\']);
for i_im = 1:num_images%1:
  
    %load in image
    im_num = pilot_images(i_im).name(1:end-4);   
    im = imread([pilot_dir 'images\' im_num '.png']);
    [rows cols] = size(im);
    
    %get mask of image to set contrast
    mask = im < 250;
    mask = imerode(mask, strel('disk', 10));
    g_max = max(im(mask));
    g_min = min(im(mask));  
    
    %Display image
    f1 = figure(...
        'windowstyle', 'normal',...
        'Units', 'pixels',...
        'position', [0 0 1024 1024*cols/rows],...
        'PaperPositionMode','auto',...
        'Visible', 'off'); 
    
    axes( 'Units', 'normalized', 'position', [0 0 1 1]);
        
    imgray(im); caxis([g_min g_max]); axis off;
    
    vessels_d = u_load([pilot_dir 'apex_clusters\' im_num '_apex_clusters.mat']);
    for i_v = 1:size(vessels_d.cluster_centres,1);   
        
        switch vessels_d.majority_sizes{i_v}
            case 'Normal'
                marker_size = 8;
            case 'Enlarged'
                marker_size = 16;
            case 'Giant'
                marker_size = 24;
            case 'Distal'
                marker_size = 4;
            case 'Undefined'
                marker_size = 16;
        end
        
        switch vessels_d.majority_shapes{i_v}
            
            case 'Normal'
                plot(...
                    vessels_d.cluster_centres(i_v,1), vessels_d.cluster_centres(i_v,2),...
                    'ro', 'MarkerSize', marker_size);
                
            case 'Meandering'
%                 text(...
%                     vessels_d.cluster_centres(i_v,1) - marker_size, vessels_d.cluster_centres(i_v,2),...
%                     '~', 'Color', 'r', 'fontsize', 2*marker_size);
                plot(...
                    vessels_d.cluster_centres(i_v,1), vessels_d.cluster_centres(i_v,2),...
                    'r^', 'MarkerSize', marker_size);
                
            case 'Angiogenetic'
%                 text(...
%                     vessels_d.cluster_centres(i_v,1) - marker_size/2, vessels_d.cluster_centres(i_v,2),...
%                     'm', 'Color', 'r', 'fontsize', marker_size);
                plot(...
                    vessels_d.cluster_centres(i_v,1), vessels_d.cluster_centres(i_v,2),...
                    'rs', 'MarkerSize', marker_size);
                
            case 'Distal'
                plot(...
                    vessels_d.cluster_centres(i_v,1), vessels_d.cluster_centres(i_v,2),...
                    'cx', 'MarkerSize', marker_size);
                
            case 'Undefined'
                text(...
                    vessels_d.cluster_centres(i_v,1) - marker_size/2, vessels_d.cluster_centres(i_v,2),...
                    '?', 'Color', 'g', 'fontsize', marker_size);
        end  
        
    end
    
    fname = [pilot_dir 'apex_clusters\display_images\' im_num '_apex_clusters'];

    print('-depsc2', [fname '.eps'], f1);
    eval(['!convert -density 150 "', fname,'.eps" "',fname,'.png"']);
    delete([fname '.eps']);
    close(f1);
end
%%
%Display image
f1 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [0 0 150 200],...
    'PaperPositionMode','auto',...
    'Visible', 'on'); 

axes( 'Units', 'normalized', 'position', [0 0 1 1]);
hold on;
axis off ij equal;

plot(10, 28, 'ro', 'markersize', 8);
plot(60, 24, 'ro', 'markersize', 16);
plot(110, 20, 'ro', 'markersize', 24);

plot(10, 68, 'r^', 'markersize', 8);
plot(60, 64, 'r^', 'markersize', 16);
plot(110, 60, 'r^', 'markersize', 24);

plot(10, 108, 'rs', 'markersize', 8);
plot(60, 104, 'rs', 'markersize', 16);
plot(110, 100, 'rs', 'markersize', 24);

plot(10, 140, 'cx', 'markersize', 8);
text(5, 180, '?', 'Color', 'g', 'fontsize', 16);
axis([0 150 0 200]);

fname = [pilot_dir 'apex_clusters\display_images\key'];

print('-depsc2', [fname '.eps'], f1);
eval(['!convert -density 150 "', fname,'.eps" "',fname,'.png"']);
delete([fname '.eps']);
close(f1);
