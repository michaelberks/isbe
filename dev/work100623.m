grp1 = mvnrnd([ 1  1],   eye(2), 100);
grp2 = mvnrnd([-1 -1], 2*eye(2), 100);

figure; hold on;
plot(grp1(:,1), grp1(:,2), 'rx');
plot(grp2(:,1), grp2(:,2), 'b+');

% K = 10;
% 
% num_grp1 = size(grp1,1);
% d_sorted_grp1 = zeros(num_grp1,K);
% d_index_grp1 = zeros(num_grp1,K);
% 
% d_sorted_grp2 = zeros(num_grp1,K);
% d_index_grp2 = zeros(num_grp1,K);
% 
% for i = 1:num_grp1
%     Dk = sum(bsxfun(@minus, grp1([1:i-1 i+1:end],:), grp1(i,:)).^2, 2);
%     if K>1
%         [sorted,index] = sort(Dk);
%         d_sorted_grp1(i,:) = sorted(1:K);
%         d_index_grp1(i,:) = index(1:K);
%     else
%         [d_sorted_grp1(i,:),d_index_grp1(i,:)] = min(Dk);
%     end
%     Dk = sum(bsxfun(@minus, grp2, grp1(i,:)).^2, 2);
%     if K>1
%         [sorted,index] = sort(Dk);
%         d_sorted_grp2(i,:) = sorted(1:K);
%         d_index_grp2(i,:) = index(1:K);
%     else
%         [d_sorted_grp2(i,:),d_index_grp2(i,:)] = min(Dk);
%     end
% end
% 
% grp1_to_grp1 = mean(d_sorted_grp1,2);
% grp1_to_grp2 = mean(d_sorted_grp2,2);
% 
% dist_ratio = grp1_to_grp2 ./ grp1_to_grp1;

[dist_ratio, grp1_to_grp1, grp1_to_grp2] = knn_2_sample_distance(grp1, grp2, 10);

idx = find(grp1_to_grp1 > grp1_to_grp2);
plot(grp1(idx,1), grp1(idx,2), 'go');
%%
sample_args.pair_name = 'M:\asymmetry_project\data\contralateral\002LCC';
sample_args.win_size = 3;
sample_args.num_levels = 6;
sample_args.do_max = 0;
sample_args.num_train = 2e3;
sample_args.num_test = 1e3;
sample_args.feature_type = 'all';

[training_data training_labels test_data test_labels] = sample_contralateral_data(sample_args);
%%
sample_args.pair_name = 'M:\asymmetry_project\data\contralateral\002LCC';
sample_args.win_size = 3;
sample_args.num_levels = 6;
sample_args.do_max = 0;
sample_args.num_train = 2e3;
sample_args.num_test = 1e3;
sample_args.feature_type = 'all';

forest_args.sampling_method = 'sample_contralateral_data';
forest_args.sampling_method_args = sample_args;
forest_args.d = 20;
forest_args.split_min = 20;
forest_args.tree_dir = 'M:\asymmetry_project\results\contralateral_rfs\002LCC_rf\';
forest_args.save_path = 'M:\asymmetry_project\results\contralateral_rfs\002LCC_rf\random_forest.mat';
forest_args.do_err = 1;
forest_args.n_trees = 5;

[random_forest] = mb_random_forest_class_train(forest_args);
%%
contra_pair = u_load('M:\asymmetry_project\data\contralateral\002LCC');
map_args.region1 = imresize(contra_pair.abnormal_roi, 0.5, 'bilinear');
map_args.region2 = imresize(contra_pair.normal_roi, 0.5, 'bilinear');
map_args.num_samples2 = 1e3;
map_args.k = 10;
clear contra_pair; pack;
dist_map = knn_2_region_map(map_args);
save('M:\asymmetry_project\results\contralateral_data\knn_maps\002LCC_knn_map.mat', 'dist_map');
%%
contra_pair = u_load('M:\asymmetry_project\data\contralateral\002LCC');
dist_map = u_load('M:\asymmetry_project\results\contralateral_data\knn_maps\002LCC.mat');
base_im = double(imresize(contra_pair.abnormal_roi, 0.5, 'bilinear'));
base_im = repmat((base_im-min(base_im(:))) / (max(base_im(:)) - min(base_im(:))), [1 1 3]);
hot_colors = hot(256);
hot_map = (dist_map - min(dist_map(:))) / (max(dist_map(:)) - min(dist_map(:)));
hot_map = round(255*hot_map)+1;
hot_im = reshape(hot_colors(hot_map, :), [size(hot_map,1) size(hot_map,2) 3]);

figure; image(0.3*hot_im + 0.7*base_im); axis image;
%%
% Load in part maps computed on hydra
for ii = 1:25
    dist_map_part = u_load(['\\isbe-san1\mberks\asymmetry_project\results\knn_maps\024RML_knn_map' zerostr(ii, 2) '.mat']);
    if ii == 1
        dist_map_all = zeros(size(dist_map_part));
        
    end
    dist_map_all = dist_map_all + dist_map_part;
end
figure; imagesc(dist_map_all); axis image;

contra_pair = u_load('M:\asymmetry_project\data\contralateral\024RML');
base_im = double(imresize(contra_pair.abnormal_roi, 0.5, 'bilinear'));
base_im = repmat((base_im-min(base_im(:))) / (max(base_im(:)) - min(base_im(:))), [1 1 3]);
hot_colors = hot(256);
hot_map = (dist_map_all - min(dist_map_all(:))) / (max(dist_map_all(:)) - min(dist_map_all(:)));
hot_map = round(255*hot_map)+1;
hot_im = reshape(hot_colors(hot_map, :), [size(hot_map,1) size(hot_map,2) 3]);

figure; image(0.3*hot_im + 0.7*base_im); axis image;

imwrite(hot_im, 'M:\asymmetry_project\results\contralateral_data\knn_maps\figures\024RML_knn_map.bmp');
imwrite(0.4*hot_im + 0.6*base_im, 'M:\asymmetry_project\results\contralateral_data\knn_maps\figures\024RML_knn_overlay.bmp');
%%
% Load in part maps computed on hydra
for ii = 1:25
    dist_map_part = u_load(['\\isbe-san1\mberks\asymmetry_project\results\knn_maps\024RMLn_knn_map_W3L5N1e4K10_' zerostr(ii, 2) '.mat']);
    if ii == 1
        dist_map_all = zeros(size(dist_map_part));
        
    end
    dist_map_all = dist_map_all + dist_map_part;
end
figure; imagesc(dist_map_all); axis image; colormap(hot(256));

contra_pair = u_load('M:\asymmetry_project\data\contralateral\024RML');
base_im = double(imresize(contra_pair.normal_roi, 0.5, 'bilinear'));
base_im = repmat((base_im-min(base_im(:))) / (max(base_im(:)) - min(base_im(:))), [1 1 3]);
hot_colors = hot(256);
hot_map = (dist_map_all - min(dist_map_all(:))) / (max(dist_map_all(:)) - min(dist_map_all(:)));
hot_map = round(255*hot_map)+1;
hot_im = reshape(hot_colors(hot_map, :), [size(hot_map,1) size(hot_map,2) 3]);

figure; image(0.3*hot_im + 0.7*base_im); axis image;
%%
rf_map = u_load('M:\asymmetry_project\results\predict_results\image_pair027_50kMSL200\probability_image027.mat');
rf_map = imresize(rf_map, 0.5, 'bilinear');
hot_map = (rf_map - min(rf_map(:))) / (max(rf_map(:)) - min(rf_map(:)));
hot_map = round(255*hot_map)+1;
hot_im = reshape(hot_colors(hot_map, :), [size(hot_map,1) size(hot_map,2) 3]);
figure; image(0.3*hot_im + 0.7*base_im); axis image;

imwrite(hot_im, 'M:\asymmetry_project\results\contralateral_data\knn_maps\figures\024RML_rf_map.bmp');
imwrite(0.3*hot_im + 0.7*base_im, 'M:\asymmetry_project\results\contralateral_data\knn_maps\figures\024RML_rf_overlay.bmp');
%%
for ii = [1 27 54 81]
    contra_pair = u_load(['M:\asymmetry_project\data\contralateral\' c_list(ii).name]);
    write_im_from_colormap(contra_pair.abnormal_roi,...
        ['M:\asymmetry_project\data\contralateral\figures\' c_list(ii).name(1:end-4) '_abnormal.jpg'], gray(256));
    write_im_from_colormap(contra_pair.normal_roi,...
        ['M:\asymmetry_project\data\contralateral\figures\' c_list(ii).name(1:end-4) '_normal.jpg'], gray(256));
end
