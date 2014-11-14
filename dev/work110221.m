%%
load C:\isbe\asymmetry_project\data\misc\reg_data1e5.mat y
uv = complex(cosd(2*y), sind(2*y));
N = length(uv);
[uv_sorted uv_sort_idx] = sort(abs(angle(conj(uv(1)) * uv)));%

colors = hsv(180);
colors = [colors; colors(1,:)];
figure;
%%
f1 = figure; 
subplot(1,2,1); 
hist(abs(tree_old.class(leaves_old)), linspace(0,1,50));
title('Sum of L and R squares');
xlabel(['Mean = ' num2str(mean(abs(tree_old.class(leaves_old)))) ...
    ' Median = ' num2str(median(abs(tree_old.class(leaves_old))))]);


subplot(1,2,2); hist(abs(tree_new.class(leaves_new)), linspace(0,1,50));
title('Difference of L and R magnitudes');
xlabel(['Mean = ' num2str(mean(abs(tree_new.class(leaves_new)))) ...
    ' Median = ' num2str(median(abs(tree_new.class(leaves_new))))]);
%%
for ii = 1:200
    
    %Make dir to save sampled to
    save_dir = ['Z:/data/synthetic_data/real512_dt/' zerostr(ii,3) '/'];

    y = u_load([save_dir 'y.mat']);
    for jj = 1:20
        keep_rows = (1:1e4) + (jj-1)*1e4;
        y_10k = y(keep_rows,:); %#ok
        save([save_dir 'y_' zerostr(jj,2) '.mat'], 'y_10k');
    end
    delete([save_dir 'y.mat']);
end
%%
for ii = 1:10
    for jj = 1:20
        tree = u_load(['Z:\data\line_orientation_rfs\290903\' zerostr(ii,2) '_trees\rf_tree' zerostr(jj,4) '.mat']);
        tree = rmfield(tree, {'goodness', 'outputs', 'nodeopts'});
        save(['Z:\data\line_orientation_rfs\290903\rf_tree' zerostr((ii-1)*20+jj,4) '.mat']);
        delete(['Z:\data\line_orientation_rfs\290903\' zerostr(ii,2) '_trees\rf_tree' zerostr(jj,4) '.mat']);
    end
end
%%
for kk = 11:200
    load(['Z:\data\line_orientation_rfs\290903\rf_tree' zerostr(kk,4) '.mat']);
    save(['Z:\data\line_orientation_rfs\290903\rf_tree' zerostr(kk,4) '.mat'], 'tree');
end
%%
for kk = 1:10
    load(['Z:\data\line_orientation_rfs\290899\rf_tree' zerostr(kk,4) '.mat']);
    tree = rmfield(tree, {'goodness', 'outputs', 'nodeopts'});
    save(['Z:\data\line_orientation_rfs\290903\rf_tree' zerostr(kk,4) '.mat'], 'tree');
end
%%
rf_mag = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\243629\random_forest.mat');
rf_phase = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\243630\random_forest.mat');
rf_all = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\243311\random_forest.mat');
rf_or = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\286712\random_forest.mat');

[ang_mag d_mag dang_mag] = forest_leaf_output(rf_mag, 0);
[ang_phase d_phase dang_phase] = forest_leaf_output(rf_phase, 0);
[ang_all d_all dang_all] = forest_leaf_output(rf_all, 0);
[ang_or d_or dang_or] = forest_leaf_output(rf_or, 0);
%
figure; hold on;
plot(1:180, ang_mag ./ sum(ang_mag), 'r');
plot(1:180, ang_phase ./ sum(ang_phase), 'g');
plot(1:180, ang_all ./ sum(ang_all), 'b');
plot(1:180, ang_or ./ sum(ang_or), 'm');
y_lim = get(gca, 'ylim');
for bb = 15:30:165
    plot([bb bb], y_lim, 'k:');
end

figure; hold on;
plot(1:180, dang_mag ./ sum(dang_mag), 'r');
plot(1:180, dang_phase ./ sum(dang_phase), 'g');
plot(1:180, dang_all ./ sum(dang_all), 'b');
plot(1:180, dang_or ./ sum(dang_or), 'm');
y_lim = get(gca, 'ylim');
for bb = 15:30:165
    plot([bb bb], y_lim, 'k:');
end

figure; hold on;
plot(1:180, dang_mag ./ ang_mag, 'r');
plot(1:180, dang_phase ./ ang_phase, 'g');
plot(1:180, dang_all ./ ang_all, 'b');
plot(1:180, dang_or ./ ang_or, 'm');
y_lim = get(gca, 'ylim');
for bb = 15:30:165
    plot([bb bb], y_lim, 'k:');
end

figure; bar(linspace(0,1,50), d_mag);
figure; bar(linspace(0,1,50), d_phase);
figure; bar(linspace(0,1,50), d_all);
figure; bar(linspace(0,1,50), d_or);
%%
random_forest = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\243311\random_forest.mat');
bg = u_load('C:\isbe\asymmetry_project\data\synthetic_backgrounds\real512\train\bg00001.mat');
mask = false(512);
mask(256,:) = 1;
sampling_args_c.num_levels = 4;
sampling_args_c.feature_shape = 'rect';
sampling_args_c.feature_type = 'conj';
sampling_args_c.do_max = 0;
sampling_args_c.rotate = 0;
sampling_args_c.win_size = 3;
sampling_args_c.use_nag = 0;
%%
tic;
[probability_image2] = classify_image(...
    'image_in', bg,...
    'sampling_args', sampling_args_c,...
    'forest', random_forest,...
    'decomp_type', 'dt',...
    'forest_type', 'orientation',...
    'use_probs', 0,...
    'mask', mask,...
    'num_trees', 1,...
    'max_size', 512);
toc;