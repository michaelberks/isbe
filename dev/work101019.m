ab_sum = 0;
ab_N = 0;
norm_sum = 0;
norm_N = 0;
ab_means = [];
norm_means = [];
    
ab_list = dir('Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\abnormals\results\lines\*.mat');
norm_list = dir('Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\normals\results\lines\*.mat');

for ii = 1:length(norm_list)

    ab_map = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\abnormals\results\lines\' ab_list(ii).name]);
    norm_map = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\normals\results\lines\' norm_list(ii).name]);

    if strcmpi(ab_list(ii).name(4), 'R')
        ab_map = fliplr(ab_map);
    else
        norm_map = fliplr(norm_map);
    end
    ab_mask = u_load(['Z:\asymmetry_project\data\masks\2004_screening\contralateral_roi\abnormals\' ab_list(ii).name(1:6) '_mask.mat']);
    norm_mask = u_load(['Z:\asymmetry_project\data\masks\2004_screening\contralateral_roi\normals\' norm_list(ii).name(1:6) '_mask.mat']);       

    ab_im_sum = sum(1 - ab_map(ab_mask));
    ab_im_N = sum(ab_mask(:));
    norm_im_sum = sum(1 - norm_map(norm_mask));
    norm_im_N = sum(norm_mask(:));

    %display(['Image ' num2str(ii) ': ' num2str([ab_im_sum/ab_im_N norm_im_sum/norm_im_N], 3) ]); 

    ab_sum = ab_sum + ab_im_sum;
    ab_N = ab_N + ab_im_N;
    norm_sum = norm_sum + norm_im_sum;
    norm_N = norm_N + norm_im_N;

    ab_means(end+1) = ab_im_sum/ab_im_N; %#ok
    norm_means(end+1) = norm_im_sum/norm_im_N; %#ok

end

display(ab_sum / ab_N);
display(norm_sum / norm_N);
[h p] = ttest(ab_means, norm_means) %#ok
sum(ab_means > norm_means) %#ok
%%
roi_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\contralateral_abnormal_roi\*.mat');
ab_list = dir('Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\abnormals\results\lines\*.mat');
norm_list = dir('Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\normals\results\lines\*.mat');

ab_names = get_mammo_info(roi_list);
for ii = 11:30;
    load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\contralateral_abnormal_roi\' roi_list(ii).name]);
    
    ab_line = u_load(['C:\isbe\asymmetry_project\data\line_maps\2004_screening\abnormals\' ab_names{ii} '_data.mat']);
    
    contralateral_pair.abnormal_pos = round(contralateral_pair.abnormal_pos/2);
    
    ra1 = contralateral_pair.abnormal_pos(1,2);
    ra2 = contralateral_pair.abnormal_pos(2,2);
    ca1 = contralateral_pair.abnormal_pos(1,1);
    ca2 = contralateral_pair.abnormal_pos(2,1);

    ab_line = ab_line(ra1:ra2, ca1:ca2);
    
    norm_name = ab_names{ii};
    if contralateral_pair.right
        norm_name(4) = 'L';
    else
        norm_name(4) = 'R';
    end
    
    norm_line = u_load(['C:\isbe\asymmetry_project\data\line_maps\2004_screening\abnormals\' norm_name '_data.mat']);
    
    contralateral_pair.normal_pos = round(contralateral_pair.normal_pos/2);
    
    rn1 = contralateral_pair.normal_pos(1,2);
    rn2 = contralateral_pair.normal_pos(2,2);
    cn1 = contralateral_pair.normal_pos(1,1);
    cn2 = contralateral_pair.normal_pos(2,1);

    norm_line = norm_line(rn1:rn2, cn1:cn2);

    if contralateral_pair.right
        ab_line = fliplr(ab_line);
    else
        norm_line = fliplr(norm_line);
    end
    
    ab_map = 1-u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\abnormals\results\lines\' ab_list(ii).name]);
    norm_map = 1-u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\normals\results\lines\' norm_list(ii).name]);
    
    [m n] = size(ab_map);
    ab_hsv = ones(m, n, 3);
    ab_hsv(:,:,2) = ab_map.^2;
    ab_hsv(:,:,3) = ab_line.^2;
    ab_rgb = hsv2rgb(ab_hsv);
    
    [m n] = size(norm_map);
    norm_hsv = ones(m, n, 3);
    norm_hsv(:,:,2) = norm_map.^2;
    norm_hsv(:,:,3) = norm_line.^2;
    norm_rgb = hsv2rgb(norm_hsv);
    
    figure;
    subplot(1,2,1); image(ab_rgb); axis image;
    subplot(1,2,2); image(norm_rgb); axis image;
 
%     figure; colormap(jet(256));
%     subplot(2,2,1); imagesc(ab_map); axis image; caxis([0 1]);
%     subplot(2,2,2); imagesc(norm_map); axis image; caxis([0 1]);
%     subplot(2,2,3); imagesc(ab_line); axis image;
%     subplot(2,2,4); imagesc(norm_line); axis image;
    
    
end
%%

mam_list = dir('Z:\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\*mat');
mam_names = get_mammo_info(mam_list);

for ii = 1:20
    ori_map = u_load(['Z:\asymmetry_project\data\orientation_maps\2004_screening\abnormals\' mam_names{ii} '_data.mat']);
    line_map = u_load(['Z:\asymmetry_project\data\line_maps\2004_screening\abnormals\' mam_names{ii} '_data.mat']);
    weight_map = u_load(['Z:\asymmetry_project\data\relevance_maps\2004_screening\abnormals\' mam_names{ii} '_class.mat']);

    meta = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\' mam_names{ii} '_meta.mat']);

    combined_map = weight_map .* line_map .* exp(i*pi*ori_map/180);

    
%     %Display quiver plots
%     mammo = u_load('Z:\asymmetry_project\data\mammograms\2004_screening\abnormals\o04_024RML.mat');
%     mask = u_load('Z:\asymmetry_project\data\masks\2004_screening\abnormals\o04_024RML_mask.mat');
%     display_orientation(mammo, ori_map, 32, mask);
% 
%     fun = @(x) max(x(:));
%     spacing = 16;
%     max_map = blkproc(combined_map, [spacing spacing], fun);
% 
%     figure; imagesc(mammo); axis image; hold on;
%     [y x] = size(mammo);
%     quiver(1:spacing:x, 1:spacing:y, real(max_map), -imag(max_map), 2);

    figure; image(complex2rgb(combined_map(1:4:end, 1:4:end).^2)); axis image; hold on;
    plot(meta(:,1)/4, meta(:,2)/4, 'r', 'linewidth', 1);
    
    clear *map
end
%%
forest = u_load('Z:\asymmetry_project\data\line_width_rfs\208194\random_forest.mat');
sampling_args = u_load('Z:\asymmetry_project\data\line_width_rfs\208194\sampling_args.mat');
sampling_args_c.num_levels = sampling_args.num_levels;
sampling_args_c.feature_shape = sampling_args.feature_shape;
sampling_args_c.feature_type = sampling_args.feature_type;
sampling_args_c.do_max = sampling_args.do_max;
sampling_args_c.rotate = sampling_args.rotate;
sampling_args_c.win_size = sampling_args.win_size;
sampling_args_c.use_nag = sampling_args.use_nag;
forest.tree_root = 'Z:\asymmetry_project\data\line_width_rfs\';
    
bg = u_load('C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\bg00013.mat');
%%        
for width = 8
    ori = ceil(rand*180);
    [bar label label_centre] = create_sin_bar(width/2, 8, ori, 512, 512, 0.5, 256, 256);
    bar_bg = bar+bg;
    
    [label_width] = classify_image(...
            'image_in', bar_bg, ...
            'forest', forest,...
            'sampling_args', sampling_args_c,...
            'forest_type', 'regression',...
            'decomp_type', 'dt',...
            'mask', label_centre,...
            'num_trees', [], ...
            'max_size', 512,...
            'use_probs', 0);
        
    err = mean((width-label_width(label_width>0)).^2);
    display(['Error for width ' num2str(width) ' = ' num2str(err)]);
end
%%
centre_pts = [];
edge_pts1 = [];
edge_pts2 = [];
centre_pts_idx = find(label_centre); 
[centre_pts(:,2) centre_pts(:,1)] = ind2sub(size(label_centre), centre_pts_idx);

edge_pts1(:,1) = centre_pts(:,1) - label_width(centre_pts_idx)*sin(pi*ori/180)/2;
edge_pts1(:,2) = centre_pts(:,2) - label_width(centre_pts_idx)*cos(pi*ori/180)/2;

edge_pts2(:,1) = centre_pts(:,1) + label_width(centre_pts_idx)*sin(pi*ori/180)/2;
edge_pts2(:,2) = centre_pts(:,2) + label_width(centre_pts_idx)*cos(pi*ori/180)/2;

figure; imagesc(bar_bg); axis image; hold on;
plot(centre_pts(:,1), centre_pts(:,2), 'g.');
plot(edge_pts1(:,1), edge_pts1(:,2), 'r.');
plot(edge_pts2(:,1), edge_pts2(:,2), 'r.');
%%
forest = u_load('Z:\asymmetry_project\data\line_width_rfs\208194\random_forest.mat');
forest.tree_root = 'Z:\asymmetry_project\data\line_width_rfs\';

sampling_args = u_load('Z:\asymmetry_project\data\line_width_rfs\208194\sampling_args.mat');
sampling_args_c.num_levels = sampling_args.num_levels;
sampling_args_c.feature_shape = sampling_args.feature_shape;
sampling_args_c.feature_type = sampling_args.feature_type;
sampling_args_c.do_max = sampling_args.do_max;
sampling_args_c.rotate = sampling_args.rotate;
sampling_args_c.win_size = sampling_args.win_size;
sampling_args_c.use_nag = sampling_args.use_nag;

%%    
for ii = 1:10
    load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\image' zerostr(ii,3) '.mat']);
    load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\labels\label' zerostr(ii,3) '.mat']);
    ori_map = u_load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191934\probability_image' zerostr(ii,3) '.mat']);
    
    [label_width] = classify_image(...
            'image_in', test_image, ...
            'forest', forest,...
            'sampling_args', sampling_args_c,...
            'forest_type', 'regression',...
            'decomp_type', 'dt',...
            'mask', label_centre,...
            'num_trees', [], ...
            'max_size', 512,...
            'use_probs', 0);
        
    centre_pts = [];
    edge_pts1 = [];
    edge_pts2 = [];
    centre_pts_idx = find(label_centre); 
    [centre_pts(:,2) centre_pts(:,1)] = ind2sub(size(label_centre), centre_pts_idx);

    edge_pts1(:,1) = centre_pts(:,1) - label_width(centre_pts_idx).*sin(pi*ori_map(centre_pts_idx)/180)/2;
    edge_pts1(:,2) = centre_pts(:,2) - label_width(centre_pts_idx).*cos(pi*ori_map(centre_pts_idx)/180)/2;

    edge_pts2(:,1) = centre_pts(:,1) + label_width(centre_pts_idx).*sin(pi*ori_map(centre_pts_idx)/180)/2;
    edge_pts2(:,2) = centre_pts(:,2) + label_width(centre_pts_idx).*cos(pi*ori_map(centre_pts_idx)/180)/2;
    
    figure; imagesc(test_image); axis image; hold on;
    plot(centre_pts(:,1), centre_pts(:,2), 'g.', 'markersize', 2);
    plot(edge_pts1(:,1), edge_pts1(:,2), 'r.', 'markersize', 2);
    plot(edge_pts2(:,1), edge_pts2(:,2), 'r.', 'markersize', 2); 
end
%%
for ii = 1:10
    load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\image' zerostr(ii,3) '.mat']);
    load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\labels\label' zerostr(ii,3) '.mat']);
    ori_map = u_load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191934\probability_image' zerostr(ii,3) '.mat']);
    line_map = u_load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191905\probability_image' zerostr(ii,3) '.mat']);
    
    nms_map = mb_non_maximal_supp(line_map, ori_map, 1);
    final_map = hysterisis(line_map, [], 0.95);
    
    [label_width] = classify_image(...
            'image_in', test_image, ...
            'forest', forest,...
            'sampling_args', sampling_args_c,...
            'forest_type', 'regression',...
            'decomp_type', 'dt',...
            'mask', final_map,...
            'num_trees', [], ...
            'max_size', 512,...
            'use_probs', 0);
        
    centre_pts = [];
    edge_pts1 = [];
    edge_pts2 = [];
    centre_pts_idx = find(final_map); 
    [centre_pts(:,2) centre_pts(:,1)] = ind2sub(size(line_map), centre_pts_idx);

    edge_pts1(:,1) = centre_pts(:,1) - label_width(centre_pts_idx).*sin(pi*ori_map(centre_pts_idx)/180)/2;
    edge_pts1(:,2) = centre_pts(:,2) - label_width(centre_pts_idx).*cos(pi*ori_map(centre_pts_idx)/180)/2;

    edge_pts2(:,1) = centre_pts(:,1) + label_width(centre_pts_idx).*sin(pi*ori_map(centre_pts_idx)/180)/2;
    edge_pts2(:,2) = centre_pts(:,2) + label_width(centre_pts_idx).*cos(pi*ori_map(centre_pts_idx)/180)/2;
    
    figure;
    subplot(2,2,1); imagesc(line_map); axis image;
    subplot(2,2,2); imagesc(nms_map); axis image;
    subplot(2,2,3); imagesc(final_map); axis image;
    subplot(2,2,4); imagesc(label); axis image; hold on;
    plot(centre_pts(:,1), centre_pts(:,2), 'g.', 'markersize', 2);
    plot(edge_pts1(:,1), edge_pts1(:,2), 'r.', 'markersize', 2);
    plot(edge_pts2(:,1), edge_pts2(:,2), 'r.', 'markersize', 2); 
end
%%
for ii = 1:10
    %load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\image' zerostr(ii,3) '.mat']);
    %load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\labels\label' zerostr(ii,3) '.mat']);
    ori_map = u_load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191934\probability_image' zerostr(ii,3) '.mat']);
    line_map = u_load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191905\probability_image' zerostr(ii,3) '.mat']);
    
    nms_map = mb_non_maximal_supp(line_map, ori_map, 1);
    final_map = hysterisis(line_map, 0.9, 0, 0.6);
    
    
    figure;
    subplot(1,2,1); imagesc(nms_map); axis image;
    subplot(1,2,2); imagesc(final_map); axis image;
 
end
%%
mkdir C:\isbe\mammograms\new_CAD\BMP_2004\match
prior_list = dir('C:\isbe\mammograms\new_CAD\BMP_2001\*.bmp');
[prior_num d d] = get_mammo_info(prior_list);
prior_num = unique(prior_num);
for ii = 1:length(prior_num)
    curr_list = dir(['C:\isbe\mammograms\new_CAD\BMP_2004\*' prior_num{ii} '*.bmp']);
    
    for jj = 1:length(curr_list)
        movefile(...
            ['C:\isbe\mammograms\new_CAD\BMP_2004\' curr_list(jj).name],...
            ['C:\isbe\mammograms\new_CAD\BMP_2004\match\' curr_list(jj).name]);
    end
end
%%
mkdir C:\isbe\mammograms\new_CAD\BMP_2004\lcc
mkdir C:\isbe\mammograms\new_CAD\BMP_2004\lml
mkdir C:\isbe\mammograms\new_CAD\BMP_2004\rcc
mkdir C:\isbe\mammograms\new_CAD\BMP_2004\rml

curr_list = dir('C:\isbe\mammograms\new_CAD\BMP_2004\*LCC*.bmp');
for jj = 1:length(curr_list)
    movefile(...
        ['C:\isbe\mammograms\new_CAD\BMP_2004\' curr_list(jj).name],...
        ['C:\isbe\mammograms\new_CAD\BMP_2004\lcc\' curr_list(jj).name]);
end     
curr_list = dir('C:\isbe\mammograms\new_CAD\BMP_2004\*LML*.bmp');
for jj = 1:length(curr_list)
    movefile(...
        ['C:\isbe\mammograms\new_CAD\BMP_2004\' curr_list(jj).name],...
        ['C:\isbe\mammograms\new_CAD\BMP_2004\lml\' curr_list(jj).name]);
end        
curr_list = dir('C:\isbe\mammograms\new_CAD\BMP_2004\*RCC*.bmp');
for jj = 1:length(curr_list)
    movefile(...
        ['C:\isbe\mammograms\new_CAD\BMP_2004\' curr_list(jj).name],...
        ['C:\isbe\mammograms\new_CAD\BMP_2004\rcc\' curr_list(jj).name]);
end   
curr_list = dir('C:\isbe\mammograms\new_CAD\BMP_2004\*RML*.bmp');
for jj = 1:length(curr_list)
    movefile(...
        ['C:\isbe\mammograms\new_CAD\BMP_2004\' curr_list(jj).name],...
        ['C:\isbe\mammograms\new_CAD\BMP_2004\rml\' curr_list(jj).name]);
end    