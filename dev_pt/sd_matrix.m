clc; clear all;
figpath = [asymmetryroot,'figs/'];

clear tree;

% load most recently created tree if none specified
forest_root		= [asymmetryroot,'data/line_orientation_rfs/'];
forest_dir		= dir([forest_root,'pc*']);
forest_job		= [forest_dir(end).name,'/']
forest_dir		= dir([forest_root,forest_job,'/random_forest*.mat']);
forest_fname	= [forest_root,forest_job,forest_dir(end).name];
forest			= u_load(forest_fname);
tree			= u_load([forest.tree_root,forest.tree_dir,forest.trees{1}]);

load([forest.tree_root,forest.tree_dir,'traindata.mat']);

[N,d] = size(X);
angy = angle(y);

bins	= linspace(-pi,pi,d+1);
sd0		= std(X)'*ones(1,length(bins)-1);
sd_mat	= zeros(d,length(bins)-1);
for b = 1:length(bins)-1
	inds = (angy>bins(b)) & (angy<bins(b+1));
	Xi = X(inds,:);
	sd_mat(:,b) = std(Xi)';
end
sd_mat = sd_mat ./ sd0;
max_sd = max(sd_mat(:));
i_zero = round(255/max_sd);
cmap = zeros(255,3);
cmap(1:i_zero,1) = linspace(1,0,i_zero);
cmap(i_zero+1:end,2) = linspace(0,1,255-i_zero);

image(255*sd_mat/max_sd); colormap(cmap); axis image;
xlabel('orientation'); ylabel('feature');

