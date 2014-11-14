rad_list = dir('\\isbe-san1\mberks\dev\ad\*.mat');

max_rad = 0;
for ii = 1:length(rad_list)
    load(['\\isbe-san1\mberks\dev\ad\' rad_list(ii).name], 'angle_bands');
    max_rad = max(max_rad, max(angle_bands(:)));
end
save \\isbe-san1\mberks\dev\ad\max_rad.mat max_rad
save C:\isbe\dev\ad\max_rad.mat max_rad
%%
rad_list = dir('\\isbe-san1\mberks\dev\ad\*.mat');
for ii = 1:20
    load(['C:\isbe\dev\ad\' rad_list(ii).name], 'line_prob');
    load(['C:\isbe\mammograms\new_CAD\BMP_2004_half\o04_' rad_list(ii).name(1:6) '.mat']);
    figure; 
    subplot(1,2,1); imagesc(mammogram); axis image;
    subplot(1,2,2); imagesc(line_prob); axis image; colormap(gray(256));
    
    title(['Line probabilities for ' rad_list(ii).name(1:6)]);
    clear mammogram angle_bands;
    pack
end
%%
load C:\isbe\dev\ad\max_rad.mat max_rad
rad_list = dir('\\isbe-san1\mberks\dev\ad\*.mat');
for ii = 15:20
    load(['C:\isbe\dev\ad\' rad_list(ii).name], 'angle_bands');
    load(['C:\isbe\mammograms\new_CAD\BMP_2004_half\o04_' rad_list(ii).name(1:6) '.mat']);
    figure; 
    subplot(1,2,1); imagesc(mammogram); axis image;
    subplot(1,2,2); imagesc(angle_bands); axis image; caxis([0 max_rad]);
    
    title(['Radial projections for ' rad_list(ii).name(1:10)]);
    clear mammogram angle_bands;
    pack
end
%%
probability_dir = 'M:\chen\data\line_detection_mammo\results\DTCWT_W3L5\';
angle_dir = 'M:\chen\data\line_detection_mammo\results\regtree_W3L5_chen\';
max_band = 0;
max_sum = 0;
max_dist = 0;
for im = 1:5
    roi_prob = u_load([probability_dir, 'probability_image', zerostr(im,3), '.mat']);
    roi_ori = u_load([angle_dir, 'probability_image', zerostr(im,3), '.mat']);
    
    [angle_bands dist_sum] = radial_line_projection(roi_prob, roi_ori, [40 8]);
    all_bands = sum(angle_bands,3);
    
    max_band = max(max_band, max(angle_bands(:)));    
    max_sum = max(max_sum, max(all_bands(:)));
    max_dist = max(max_dist, max(dist_sum(:)));
    
    save(['C:\isbe\dev\ad\regions\region' zerostr(im, 2) '_rad_map.mat'], 'angle_bands', 'dist_sum');
%     figure; 
%     subplot(1,2,1); imagesc(all_bands); axis image;
%     subplot(1,2,2); imagesc(dist_sum); axis image;
end
%%
for im = 1:5
    
    load(['C:\isbe\dev\ad\regions\region' zerostr(im, 2) '_rad_map.mat'], 'angle_bands', 'dist_sum');
    all_bands = sum(angle_bands,3);
    
    for ang = 1:8
        write_im_from_colormap(angle_bands(:,:,ang), ...
            ['C:\isbe\dev\ad\figures\region' zerostr(im, 2) '_rad_map_A' zerostr(ang,2) '.bmp'], ...
            jet(256), [0 max_band]);
    end
    
    write_im_from_colormap(all_bands, ...
        ['C:\isbe\dev\ad\figures\region' zerostr(im, 2) '_rad_map.bmp'], ...
        jet(256), [0 max_sum]);
    write_im_from_colormap(dist_sum, ...
        ['C:\isbe\dev\ad\figures\region' zerostr(im, 2) '_dist_map.bmp'], ...
        jet(256), [0 max_dist]);
    
    for level = 3:6
        blur = 2^level;
        write_im_from_colormap(imfilter(all_bands, fspecial('disk', blur), 'symmetric'), ...
        ['C:\isbe\dev\ad\figures\region' zerostr(im, 2) '_rad_map_L' zerostr(blur, 2) '.bmp'], ...
        jet(256), [0 max_sum]);
    end
end
%%
probability_dir = 'M:\chen\data\line_detection_mammo\results\DTCWT_W3L5\';
angle_dir = 'M:\chen\data\line_detection_mammo\results\regtree_W3L5_chen\';

max_dist = 0;
for im = [1 4]
    roi_prob = rot90(u_load([probability_dir, 'probability_image', zerostr(im,3), '.mat']));
    roi_ori = rot90(u_load([angle_dir, 'probability_image', zerostr(im,3), '.mat']))+90;
    roi_ori(roi_ori > 180) = roi_ori(roi_ori > 180) - 180;
    
    [angle_bands dist_sum] = radial_line_projection(roi_prob, roi_ori, [40 8]);
    max_band = max(angle_bands(:));
    figure;
    for ang = 1:8
        subplot(2,4,ang); imagesc(angle_bands(:,:,ang)); axis image; caxis([0 max_band]);
    end 

end
%%
for im = 1:5
    
    load(['C:\isbe\dev\ad\regions\region' zerostr(im, 2) '_dist_map.mat'], 'dist_sum');
    write_im_from_colormap(dist_sum, ...
        ['C:\isbe\dev\ad\figures\region' zerostr(im, 2) '_dist_map_2.bmp'], ...
        jet(256), [0 max_dist]);
%     figure; 
%     subplot(1,2,1); imagesc(all_bands); axis image;
%     subplot(1,2,2); imagesc(dist_sum); axis image;
end
%%
probability_dir = 'M:\chen\data\line_detection_mammo\results\DTCWT_W3L5\';
angle_dir = 'M:\chen\data\line_detection_mammo\results\regtree_W3L5_chen\';

max_dist = 0;
for im = 1:5
    roi_prob = u_load([probability_dir, 'probability_image', zerostr(im,3), '.mat']);
    roi_ori = u_load([angle_dir, 'probability_image', zerostr(im,3), '.mat']);
    
    [angle_bands dist_sum] = radial_line_projection(roi_prob, roi_ori, [36 1]);
    max_dist = max(max_dist, max(dist_sum(:)));
    
    save(['C:\isbe\dev\ad\regions\region' zerostr(im, 2) '_dist_map.mat'], 'dist_sum');
%     figure; 
%     subplot(1,2,1); imagesc(all_bands); axis image;
%     subplot(1,2,2); imagesc(dist_sum); axis image;
end
%%
probability_dir = 'M:\chen\data\line_detection_mammo\results\DTCWT_W3L5\';
angle_dir = 'M:\chen\data\line_detection_mammo\results\regtree_W3L5_chen\';

max_dist = 0;
for im = [1 4]
    roi_prob = u_load([probability_dir, 'probability_image', zerostr(im,3), '.mat']);
    roi_ori = u_load([angle_dir, 'probability_image', zerostr(im,3), '.mat']);

    for g_width = [16 32 64 128]

        [angle_bands dist_sum] = radial_line_projection(roi_prob, roi_ori, [36 1], fspecial('gaussian', [1 5*g_width], g_width));
        max_int = max(dist_sum(:));

        figure; ii = 1;
        for f_width = [1 8 16 32]


            dist_sum_f = imfilter(dist_sum, fspecial('gaussian', 5*f_width, f_width));

            write_im_from_colormap(dist_sum_f, ...
                ['C:\isbe\dev\ad\figures\region' zerostr(im, 2) '_dist_map_' zerostr(g_width,3) '_' zerostr(f_width,3) '.bmp'], ...
                jet(256), [0 max_int]);

            subplot(2,2,ii); imagesc(dist_sum_f); axis image;
            ii = ii+1;
        end
    end
end
%%
load('M:\chen\data\regtree_W3L5_chen\tree_combine\random_forest.mat');
load('M:\chen\data\line_detection_mammo\bg006.mat');
%%
line_ori = classify_image(...
        'image_in', bg,...
        'forest', random_forest,...
        'forest_type', 'regression',...
        'num_levels', 5);
%%
line_ori_90 = classify_image(...
        'image_in', rot90(bg),...
        'forest', random_forest,...
        'forest_type', 'regression',...
        'num_levels', 5);
%%
load('C:\isbe\dev\ad\mammograms\024LML_data.mat', 'line_ori', 'line_prob');
line_ori = imresize(line_ori, 0.5, 'bilinear');
line_prob = imresize(line_prob, 0.5, 'bilinear');
pack;

for g_width = [8 16 32 64]

    [angle_bands dist_sum] = radial_line_projection(line_prob, line_ori, [36 1], fspecial('gaussian', [1 5*g_width], g_width));

    figure; ii = 1;
    for f_width = [1 4 8 16]


        dist_sum_f = imfilter(dist_sum, fspecial('gaussian', 5*f_width, f_width));

        save(...
            ['C:\isbe\dev\ad\024LML_dist_map_' zerostr(g_width,3) '_' zerostr(f_width,3) '.mat'], ...
            'dist_sum_f');

        subplot(2,2,ii); imagesc(dist_sum_f); axis image;
        ii = ii+1;
    end
end
clear all;
%
load('C:\isbe\dev\ad\mammograms\024RML_data.mat', 'line_ori', 'line_prob');
line_ori = imresize(line_ori, 0.5, 'bilinear');
line_prob = imresize(line_prob, 0.5, 'bilinear');
pack;

for g_width = [8 16 32 64]

    [angle_bands dist_sum] = radial_line_projection(line_prob, line_ori, [36 1], fspecial('gaussian', [1 5*g_width], g_width));

    figure; ii = 1;
    for f_width = [1 4 8 16]


        dist_sum_f = imfilter(dist_sum, fspecial('gaussian', 5*f_width, f_width));

        save(...
            ['C:\isbe\dev\ad\024RML_dist_map_' zerostr(g_width,3) '_' zerostr(f_width,3) '.mat'], ...
            'dist_sum_f');

        subplot(2,2,ii); imagesc(dist_sum_f); axis image;
        ii = ii+1;
    end
end