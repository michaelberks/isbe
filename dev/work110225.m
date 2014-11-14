rf_dabs = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\299400\random_forest.mat');
rf_midp = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\299402\random_forest.mat');
rf_ssq = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\299407\random_forest.mat');
rf_ssq2 = u_load('C:\isbe\asymmetry_project\data\line_orientation_rfs\299430\random_forest.mat');

rf_dabs.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
rf_midp.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
rf_ssq.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
rf_ssq2.tree_root = 'C:\isbe\asymmetry_project\data\line_orientation_rfs\';
%
[leaf_dabs] = forest_leaf_output(rf_dabs);
[leaf_midp] = forest_leaf_output(rf_midp);
[leaf_ssq] = forest_leaf_output(rf_ssq);
[leaf_ssq2] = forest_leaf_output(rf_ssq2);

%%
figure; hold on;
plot(1:180, leaf_dabs.ang_hist ./ sum(leaf_dabs.ang_hist), 'r');
plot(1:180, leaf_midp.ang_hist ./ sum(leaf_midp.ang_hist), 'g');
plot(1:180, leaf_ssq.ang_hist ./ sum(leaf_ssq.ang_hist), 'b');
%plot(1:180, leaf_ssq2.ang_hist ./ sum(leaf_ssq2.ang_hist), 'm');
y_lim = get(gca, 'ylim');
for bb = 15:30:165
    plot([bb bb], y_lim, 'k:');
end
%%
figure; hold on;
plot(1:180, leaf_dabs.dang_hist ./ leaf_dabs.ang_hist, 'r');
plot(1:180, leaf_midp.dang_hist ./ leaf_midp.ang_hist, 'g');
plot(1:180, leaf_ssq.dang_hist ./ leaf_ssq.ang_hist, 'b');
plot(1:180, leaf_ssq2.dang_hist ./ leaf_ssq2.ang_hist, 'm');
y_lim = get(gca, 'ylim');
for bb = 15:30:165
    plot([bb bb], y_lim, 'k:');
end
%
figure; 
subplot(2,2,1); bar(linspace(0,1,50), leaf_dabs.d_hist);
subplot(2,2,2); bar(linspace(0,1,50), leaf_midp.d_hist);
subplot(2,2,3); bar(linspace(0,1,50), leaf_ssq.d_hist);
subplot(2,2,4); bar(linspace(0,1,50), leaf_ssq2.d_hist);

figure; 
subplot(2,2,1); bar(linspace(2,200,100), leaf_dabs.size_hist);
subplot(2,2,2); bar(linspace(2,200,100), leaf_midp.size_hist);
subplot(2,2,3); bar(linspace(2,200,100), leaf_ssq.size_hist);
subplot(2,2,4); bar(linspace(2,200,100), leaf_ssq2.size_hist);