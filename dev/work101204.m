data_type = {'abnormals', 'normals'};

for ii = 1:2

    mkdir(['C:\isbe\asymmetry_project\data\k_line_strengths\2004_screening_processed\' data_type{ii} '\']);
    mkdir(['C:\isbe\asymmetry_project\data\k_line_maps\2004_screening_processed\' data_type{ii} '\']);
    mkdir(['C:\isbe\asymmetry_project\data\k_ori_maps\2004_screening_processed\' data_type{ii} '\']);
    mkdir(['C:\isbe\asymmetry_project\data\k_scale_maps\2004_screening_processed\' data_type{ii} '\']);
    
    mammo_dir = ['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\' data_type{ii} '\'];
    mask_dir = ['C:\isbe\asymmetry_project\data\masks\2004_screening\' data_type{ii} '\'];
    
    mammo_list = dir([mammo_dir '*.mat']);
    mammo_names = get_mammo_info(mammo_list);
    
    for jj = 1:length(mammo_list)
        display(['processing image ' num2str(jj) ' of ' num2str(length(mammo_list))]);
        mammo = imresize(u_load([mammo_dir mammo_list(jj).name]), 0.5, 'bilinear');
        mask = imresize(u_load([mask_dir mammo_names{jj} '_mask.mat']), 0.5);
        
        [line_strength, orientation_map, scale_map] = gaussian_2nd_derivative_line(mammo, [1 2 4 8]);
        [grad_strength, grad_orientation] = gaussian_1st_derivative_gradient(mammo, 5);
  
        line_map = ...
            ((abs(mb_mod(orientation_map - grad_orientation + pi/2, pi)) < pi/6) &...
            (grad_strength > 25)) | ...
            (line_strength > -1) | ...
            ~mask;
        
        line_map = ~line_map;
        orientation_map = mod(180*orientation_map/pi, 180);
    
        save(['C:\isbe\asymmetry_project\data\k_line_strengths\2004_screening_processed\' data_type{ii} '\'...
            mammo_names{jj} '_data.mat'], 'line_strength');
        save(['C:\isbe\asymmetry_project\data\k_ori_maps\2004_screening_processed\' data_type{ii} '\'...
            mammo_names{jj} '_data.mat'], 'orientation_map');
        save(['C:\isbe\asymmetry_project\data\k_scale_maps\2004_screening_processed\' data_type{ii} '\'...
            mammo_names{jj} '_data.mat'], 'scale_map');
        save(['C:\isbe\asymmetry_project\data\k_line_maps\2004_screening_processed\' data_type{ii} '\'...
            mammo_names{jj} '_data.mat'], 'line_map');
        
        %figure; imagesc(line_map); axis image;
        clear mammo line_strength line_map orientation_map grad_map grad_orientation mask;
    end
end
%%
[line_strength, orientation_map, scale_map] = gaussian_2nd_derivative_line(bg, [1 2 4 8]);
[grad_strength, grad_orientation] = gaussian_1st_derivative_gradient(bg, 5);
  
line_map = ...
    ((abs(mb_mod(orientation_map - grad_orientation + pi/2, pi)) < pi/6) &...
    (grad_strength > 25)) | ...
    (line_strength > 0);

line_map = ~line_map;
orientation_map = mod(180*orientation_map/pi, 180);

figure; colormap(hsv(180));
subplot(1,2,1); imagesc(line_map); axis image; caxis([0 2]);
subplot(1,2,2); imagesc(orientation_map); axis image; caxis([0 180]); colorbar;
%%
data_type = {'abnormals', 'normals'};
ii = 2;
mammo_dir = ['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\' data_type{ii} '\'];
mask_dir = ['C:\isbe\asymmetry_project\data\masks\2004_screening\' data_type{ii} '\'];

mammo_list = dir([mammo_dir '*ML.mat']);
mammo_names = get_mammo_info(mammo_list);
%
for jj = 0*10 + (1:10);

    load(['C:\isbe\asymmetry_project\data\k_ori_maps\2004_screening_processed\' data_type{ii} '\'...
        mammo_names{jj} '_data.mat'], 'orientation_map');

    load(['C:\isbe\asymmetry_project\data\k_line_maps\2004_screening_processed\' data_type{ii} '\'...
        mammo_names{jj} '_data.mat'], 'line_map');

    figure('Name', mammo_names{jj}); colormap(hsv(180));
    subplot(1,2,1); imagesc(line_map); axis image; caxis([0 2]);
    subplot(1,2,2); imagesc(orientation_map); axis image; caxis([0 180]);
    clear line_map orientation_map;
end
%%
for jj = 18*10 + (1:10);

    template_map = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening_processed\' data_type{ii} '\'...
        mammo_names{jj} '_template.mat']);

    scales_map = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening_processed\' data_type{ii} '\scales\'...
        mammo_names{jj} '_scales.mat']);

    figure('Name', mammo_names{jj});
    subplot(1,2,1); imagesc(template_map); axis image;
    subplot(1,2,2); imagesc(scales_map); axis image;
    clear scales_map template_map;
end
%%
for jj = 0*10 + (1:10);

    line_map = load_uint8(['Z:\asymmetry_project\data\mammograms\2004_screening_processed\' data_type{ii} '\results\191905\'...
        mammo_names{jj} '_class.mat']);

    orientation_map = load_uint8(['Z:\asymmetry_project\data\mammograms\2004_screening_processed\' data_type{ii} '\results\191934\'...
        mammo_names{jj} '_class.mat']);

    figure('Name', mammo_names{jj});
    subplot(1,2,1); imagesc(line_map); axis image;
    subplot(1,2,2); imagesc(orientation_map); axis image;
    clear line_map orientation_map;
end
%%
for jj = 9*10 + (1:10);

    mammo = imresize(u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\' data_type{ii} '\'...
        mammo_names{jj}]), 0.5, 'bilinear');

    segmentation = u_load(['C:\isbe\asymmetry_project\data\segmentations\2004_screening\' data_type{ii} '\'...
        mammo_names{jj} '_segmentation.mat']);
    
    [breast_border breast_air pectoral_edge] = segment_breast_resize(size(mammo), segmentation);

    figure('Name', mammo_names{jj});
    imagesc(mammo); axis image; colormap(gray(256)); hold all;
    plot(breast_border(:,1), breast_border(:,2), 'g');
    
    if ~isempty(pectoral_edge)
        plot(pectoral_edge(:,1), pectoral_edge(:,2), 'r');
    end
    clear mammo;
end
%%
data_type = {'abnormals', 'normals'};
    
for ii = 1:2
    mammo_dir = ['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\' data_type{ii} '\'];

    mammo_list = dir([mammo_dir '*ML.mat']);
    mammo_names = get_mammo_info(mammo_list);
    
    for jj = 1:length(mammo_names);
        
        %load segmentation
        segmentation = u_load(['C:\isbe\asymmetry_project\data\segmentations\2004_screening\' data_type{ii} '\'...
            mammo_names{jj} '_segmentation.mat']);
        
        %load K-line map
        line_map = u_load(['C:\isbe\asymmetry_project\data\k_line_maps\2004_screening_processed\' data_type{ii} '\'...
            mammo_names{jj} '_data.mat']);
        [r c] = size(line_map);
        
        %resize the segmentation to the map
        [breast_border breast_air pectoral_edge] = segment_breast_resize([r c], segmentation);

        %Make a mask of the pectoral triangle
        if strcmpi(mammo_names{jj}(4), 'R')
            sx = c;
            sy = 1;
        else
            sx = 1;
            sy = 1;
        end
        mask = poly2mask([sx; pectoral_edge(:,1)], [sy; pectoral_edge(:,2)], r, c);

        %Set the line map to zero inside the mask
        line_map(mask) = 0;
        
%         figure('Name', mammo_names{jj});
%         imagesc(line_map); axis image; colormap(gray(256));
        save(['C:\isbe\asymmetry_project\data\k_line_maps\2004_screening_processed\' data_type{ii} '\'...
            mammo_names{jj} '_data.mat'], 'line_map');
%         clear line_map;
%         
%         %Repeat for the RF maps
%         %load K-line map
%         line_map = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\2004_screening_processed\' data_type{ii} '\'...
%             mammo_names{jj} '_class.mat']);
%         [r c] = size(line_map);
%         
%         %resize the segmentation to the map
%         [breast_border breast_air pectoral_edge] = segment_breast_resize([r c], segmentation);
% 
%         %Make a mask of the pectoral triangle
%         if strcmpi(mammo_names{jj}(4), 'R')
%             sx = c;
%             sy = 1;
%         else
%             sx = 1;
%             sy = 1;
%         end
%         mask = poly2mask([sx; pectoral_edge(:,1)], [sy; pectoral_edge(:,2)], r, c);
% 
%         %Set the line map to zero inside the mask
%         line_map(mask) = 0;
%         
% %         figure('Name', mammo_names{jj});
% %         imagesc(line_map); axis image; colormap(gray(256));
%         save_uint8(['C:\isbe\asymmetry_project\data\line_maps\2004_screening_processed\' data_type{ii} '\'...
%             mammo_names{jj} '_class.mat'], line_map);
        clear line_map;
    end
end
%%
%Get path of all sub-directories
dir_path = genpath('Z:\asymmetry_project');

%Workout where each sub-directory starts and ends
dir_ends = strfind(dir_path, pathsep);
num_dirs = length(dir_ends);

%Loop through the sub-directories looking for dicom files
id2 = -1;
for ii = 1:num_dirs
    id1 = id2 + 2;
    id2 = dir_ends(ii) - 1;
    dir_name = [dir_path(id1:id2) filesep];
    o04_list = dir([dir_name 'o04_*.mat']);
    
    %For each dicom file found, try converting it's name and moving to the
    %output directory
    for jj = 1:length(o04_list)
        movefile([dir_name o04_list(jj).name], [dir_name o04_list(jj).name(5:end)]);
        display(['Moved ' dir_name o04_list(jj).name ' to ' o04_list(jj).name(5:end)]);
    end
end
%%
abnormal_names = get_mammo_info(dir('C:\isbe\asymmetry_project\data\mammograms2\2004_screening\abnormals\meta\*.mat'));
full_names = get_mammo_info(dir('C:\isbe\asymmetry_project\data\mammograms2\2004_screening\abnormals\*.mat'));
contra_names = setdiff(full_names, abnormal_names);
%%
counts = cell(6,1);
bins = cell(6,1);
for ii = 1:6
    [counts{ii} bins{ii}] = mass_maps_hist('2004_screening/abnormals', 'weighted_radial_maps', 'view', zerostr(2^(ii+2),3));
end
%%
figure; hold all;
for ii = 1:6
    plot(bins{ii}, counts{ii});
end
legend(num2str(2.^((1:6)'+3)));

%%
for ii = 1:6
    figure; hold all;
    plot(bins{ii}, 100*counts{ii} / sum(counts{ii}));
    [dummy idx] = max(counts{ii});
    
    mu = bins{ii}(idx);
    mus(ii) = mu;
    sigma2 = sum(counts{ii}(idx+1:end).*(bins{ii}(idx+1:end) - mu).^2) / sum(counts{ii}(idx+1:end));
    
    gx = exp(-((bins{ii}-mu).^2)/sigma2) / sqrt(2*pi*sigma2);
    plot(bins{ii}, gx, 'r');
    
    mu = sum(bins{ii} .* counts{ii}) / sum(counts{ii});
    sigma2 = sum(counts{ii}.*(bins{ii} - mu).^2) / sum(counts{ii});
    
    gx = exp(-((bins{ii}-mu).^2)/sigma2) / sqrt(2*pi*sigma2);
    plot(bins{ii}, gx, 'g');
end
%%
mam_names = get_mammo_info(dir('C:\isbe\asymmetry_project\data\mammograms2\2004_screening\abnormals\*.mat'));
mkdir C:\isbe\asymmetry_project\data\combined_wr_maps\2004_screening\abnormals\
mkdir C:\isbe\asymmetry_project\data\combined_wr_maps\2004_screening\abnormals\scales
for ii = 1:length(mam_names)
    map = load_uint8(['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\' mam_names{ii} '_rad_map_008.mat']);
    map = imfilter(map, fspecial('gaussian', 49, 8), 'replicate');
    scales_map = 8*ones(size(map));
    for jj = 1:5
        map_j = load_uint8(['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\' mam_names{ii} '_rad_map_' zerostr(2^(jj+3),3) '.mat']);
        map_j = imfilter(map_j, fspecial('gaussian', 49, 8), 'replicate');
        swap_idx = map_j > map;
        map(swap_idx) = map_j(swap_idx);
        scales_map(swap_idx) = 2^(jj+3);
    end
%     figure;
%     subplot(1,2,1); imagesc(map); axis image;
%     subplot(1,2,2); imagesc(log(scales_map) / log(2)); axis image;
    save_uint8(['C:\isbe\asymmetry_project\data\combined_wr_maps\2004_screening\abnormals\' mam_names{ii} '_rad_map'], map);
    save_uint8(['C:\isbe\asymmetry_project\data\combined_wr_maps\2004_screening\abnormals\scales\' mam_names{ii} '_scales'], scales_map);
end
%%
mkdir C:\isbe\asymmetry_project\data\scaled_wr_maps\2004_screening\abnormals\
scales = [30 45 60 75 90];
%
load C:\isbe\asymmetry_project\experiments\radial_maps\wr_dists.mat counts bins
cdf = cell(5,1);
for ii = 1:5
    cdf{ii} = cumsum(counts{ii+1}) / sum(counts{ii+1});
end
bins(1) = [];
%
mam_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals_contra');
for ii = 1:length(mam_names)

    map = load_uint8(['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\' mam_names{ii} '_rad_map_064.mat']);
    map = interp1(bins{3}, cdf{3}, map, 'linear');
    %map = imfilter(map, fspecial('gaussian', 49, 8), 'replicate');
    
    scale_map = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening_processed\abnormals\scales\' mam_names{ii} '_scales.mat']);
    template_map = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening_processed\abnormals\' mam_names{ii} '_template.mat']);
    template_map = template_map > 0.4;
    
    for jj = [1 2 4 5]
        map_j = load_uint8(['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\' mam_names{ii} '_rad_map_' zerostr(2^(jj+3),3) '.mat']);
%         figure;
%         subplot(1,2,1); imagesc(map_j); axis image; colorbar;
        
        map_j = interp1(bins{jj}, cdf{jj}, map_j, 'linear');
%         subplot(1,2,2); imagesc(map_j); axis image; colorbar;
        
        %map_j = imfilter(map_j, fspecial('gaussian', 49, 8), 'replicate');
        swap_idx = (scale_map == scales(jj)) & template_map;
        map(swap_idx) = map_j(swap_idx);
        
    end
%     figure; imagesc(map); axis image; colorbar;
    
    save_uint8(['C:\isbe\asymmetry_project\data\scaled_wr_maps\2004_screening\abnormals\' mam_names{ii} '_rad_map'], map);
end
%%
mam_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals_contra');
scales = [30 45 60 75 90];
for ii = 1:3%length(mam_names)
    if strcmp(mam_names{ii}(4), 'R')
        r_name = mam_names{ii};
        l_name = [mam_names{ii}(1:3) 'L' mam_names{ii}(5:6)];
    else
        l_name = mam_names{ii};
        r_name = [mam_names{ii}(1:3) 'R' mam_names{ii}(5:6)];
    end

    map_left = load_uint8(['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\' l_name '_rad_map_064.mat']);
    map_left = imfilter(map_left, fspecial('gaussian', 49, 8), 'replicate');
    map_right = load_uint8(['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\' r_name '_rad_map_064.mat']);
    map_right = imfilter(map_right, fspecial('gaussian', 49, 8), 'replicate');
    
    scale_map_left = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening_processed\abnormals\scales\' l_name '_scales.mat']);
    template_map_left = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening_processed\abnormals\' l_name '_template.mat']);
    template_map_left = template_map_left > 0.4;
    
    scale_map_right = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening_processed\abnormals\scales\' r_name '_scales.mat']);
    template_map_right = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening_processed\abnormals\' r_name '_template.mat']);
    template_map_right = template_map_right > 0.4;
    
    max_map = max(max(map_left(:)), max(map_right(:)));
    map_left = map_left / max_map;
    map_right = map_right / max_map;
    for jj = [1 2 4 5]
        map_left_j = load_uint8(['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\' l_name '_rad_map_' zerostr(2^(jj+3),3) '.mat']);
        map_left_j = imfilter(map_left_j, fspecial('gaussian', 49, 8), 'replicate');
        
        map_right_j = load_uint8(['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\' r_name '_rad_map_' zerostr(2^(jj+3),3) '.mat']);
        map_right_j = imfilter(map_right_j, fspecial('gaussian', 49, 8), 'replicate');

        max_map = max(max(map_left_j(:)), max(map_right_j(:)));
        map_left_j = map_left_j / max_map;
        map_right_j = map_right_j / max_map;
        
        
        swap_idx = (scale_map_left == scales(jj)) & template_map_left;
        map_left(swap_idx) = map_left_j(swap_idx);
        
        swap_idx = (scale_map_right == scales(jj)) & template_map_right;
        map_right(swap_idx) = map_right_j(swap_idx);
        
    end
    figure; 
    subplot(1,2,1); imagesc(map_right); axis image;
    subplot(1,2,2); imagesc(map_left); axis image;
    
    %save_uint8(['C:\isbe\asymmetry_project\data\scaled_wr_maps\2004_screening\abnormals\' mam_names{ii} '_rad_map'], map);
end
%%
map_types = {'k_stellate_maps_rf', 'k_stellate_maps_wrf', 'k_stellate_maps_g2'};
mam_types = {'\2004_screening', '\2004_screening', '\2004_screening_processed'};
% mkdir C:\isbe\asymmetry_project\data\k_stellate_maps_rf_scaled\2004_screening\normals\
% mkdir C:\isbe\asymmetry_project\data\k_stellate_maps_wrf_scaled\2004_screening\normals\
% mkdir C:\isbe\asymmetry_project\data\k_stellate_maps_g2_scaled\2004_screening_processed\normals\

%
mam_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals');
%
for kk = 1:3
    for ii = 1:length(mam_names)

        f1_map = load_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} mam_types{kk}...
            '\abnormals\' mam_names{ii} '_f1_120.mat']);
        f2_map = load_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} mam_types{kk}...
            '\abnormals\' mam_names{ii} '_f2_120.mat']);

        scale_map = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening_processed\abnormals\scales\'...
            mam_names{ii} '_scales.mat']);
        template_map = load_uint8(['C:\isbe\asymmetry_project\data\template_maps\2004_screening_processed\abnormals\'...
            mam_names{ii} '_template.mat']);
        template_map = template_map > 0.3;

        for jj = [30 45 75 90]
            f1_j = load_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} mam_types{kk}...
                '\abnormals\' mam_names{ii} '_f1_' zerostr(2*jj,3) '.mat']);
            f2_j = load_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} mam_types{kk}...
                '\abnormals\' mam_names{ii} '_f2_' zerostr(2*jj,3) '.mat']);
            
            swap_idx = (scale_map == jj) & template_map;
            
            f1_map(swap_idx) = f1_j(swap_idx);
            f2_map(swap_idx) = f2_j(swap_idx);

        end
    %     figure; imagesc(map); axis image; colorbar;

        save_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} '_scaled' mam_types{kk}...
            '\abnormals\' mam_names{ii} '_f1.mat'], f1_map);
        save_uint8(['C:\isbe\asymmetry_project\data\' map_types{kk} '_scaled' mam_types{kk}...
            '\abnormals\' mam_names{ii} '_f2.mat'], f2_map);
    end
end