clc; %clear all;
figroot = [asymmetryroot,'figs/'];

clear tree;

% profile clear; profile on;

% load most recently created forest if none specified
forest_root		= [asymmetryroot,'data/line_orientation_rfs/'];
forest_dir		= dir([forest_root,'pc*']);
forest_job		= [forest_dir(end).name,'/'];
% forest_job		= ['pc20110307T105747/']; % 20x4000
forest_job		= ['pc20110307T145627/']; % 40x16000
% forest_job		= ['pc20110307T163301/']; % 200x8000, dabs
% forest_job		= ['pc20110308T123013/']; % 200x8000, mabs
display(forest_job);
forest_dir		= dir([forest_root,forest_job,'/random_forest*.mat']);
forest_fname	= [forest_root,forest_job,forest_dir(end).name];
forest			= u_load(forest_fname);

n_trees = length(forest.trees);
n_trees = 1;
display(n_trees);

figpath = [forest_root,forest_job,'figs/'];
% make figure dir if not already there
if ~exist(figpath,'dir')
	mkdir(figpath);
end

% get sampling arguments
sampling_args	= u_load([forest_root,forest_job,'/sampling_args.mat']);

%% generate test data
if (1)
	if ~exist('X','var')
		% sample training data myself if it doesn't exist already
		fprintf('Generating training data...');
		sampling_args = rmfield(sampling_args,'detection_type');
		sampling_args.num_samples = 8000;
		sampling_args.task_id = 1;
		[X,y] = sample_saved_dt_line_data(sampling_args);
		fprintf('done\n');
	end

	if ~exist('X_test','var')
		% sample test data myself if it doesn't exist already
		fprintf('Generating test data...');
% 		sampling_args = rmfield(sampling_args,'detection_type');
		sampling_args.num_samples = 4000;
		sampling_args.task_id = 1;
		[X_test,y_test] = sample_saved_dt_line_data(sampling_args);
		fprintf('done\n');
	end
else
	% use Mike's sampled data
	X_test = [];
end

%% build boosted regressor
fprintf('Building...');
for t = 1:n_trees
% 	load([forest.tree_root,forest.tree_dir,sprintf('traindata%04d.mat',t)]);

	% build regressor
% 	boostreg = pt_boosted_reg_train(X,y,...
% 		'shrinkage',0.1,'n_levels',100,'n_bins',4);
	linreg = pt_lin_reg_train(X,y);
end
fprintf('done\n');

if exist('boostreg','var')
	inps = [boostreg.levels(:).input];
	figure(3); hist(inps,1:size(X,2));
end

%% predict

% boosted predictions for training data
% y_train = boosted_regressor_predict(boostreg,X);
y_train = linear_regressor_predict(linreg,X);
figure(1); clf; hold on;
	plot(y_train,'b.');
	plot(0.8*y,'r.');
	n = 10;
	plot(conj([y_train(1:n)'; 0.8*y(1:n)']),'c-');
	axis('equal',1.05*[-1,1,-1,1]);
% graph(1); exportfig([figpath,'boostreg_train',err_suffix]);

return

% tree predictions for training data
y_train = mb_random_forest_reg_predict(forest,X);
figure(1); clf; hold on;
	plot(y_train,'b.');
	plot(y,'r.');
	axis('equal');
graph(1); exportfig([figpath,'rf_train',err_suffix]);

fprintf('Predicting...');
y_tree = []; y_boosted = [];
if isempty(X_test)
	% use Mike's sampled data
	
	y_test = []; % all test targets in one place
	datapath = 'A:\data\synthetic_data\real512_dt\test\';
	datadir = dir([datapath,'X_*.mat']);
	tb = timebar('title','Processing datasets','limit',length(datadir));
	for i = 1:length(datadir)
		% load test data (X and y)
		X_test = u_load([datapath,datadir(i).name]);
		X_test = convert_dt_representation(X_test(:,5,:,:),...
					'feature_type',sampling_args.feature_type,...
					'win_size',sampling_args.win_size);
				
		% this data has outputs in range [0..pi]
		% double to [-pi..pi]
		y_tmp = u_load([datapath,strrep(datadir(i).name,'X_','y_')]);
		y_test = [y_test; exp(sqrt(-1)*angle(y_tmp)*2)];

% 		% predict with tree
% 		y_tree = [y_tree; mb_random_forest_reg_predict(forest,X_test)];

		% predict with boosted regressor
		y_boosted = [y_boosted; boosted_regressor_predict(boostreg,X_test)];
		
		timebar(tb,'advance');
	end
	timebar(tb,'close'); clear tb;
	
	% set error file suffix
	err_suffix = '_mb';
else
	% use my sampled data
	
	% predict with tree
	y_tree = [y_tree; mb_random_forest_reg_predict(forest,X_test)];

	% predict with boosted regressor
	y_boosted = [y_boosted; boosted_regressor_predict(boostreg,X_test)];
	
	% set error file suffix
	err_suffix = '';
end
fprintf('done\n');

% save the errors themselves
save([figpath,'outputs',err_suffix,'.mat'],...
	 'y_test','y_tree','y_boosted');

[e_tree,err_stats_tree] = ori_error(y_test,y_tree);
[e_boosted,err_stats_boosted] = ori_error(y_test,y_boosted);

%% output

% display on screen
display(err_stats_tree);
display(err_stats_boosted);

% dump errors to log file
filename = [figpath,'errors',err_suffix,'.txt'];
fid = fopen(filename,'w');
	fprintf(fid,'User = %s\n',get_username());
	fprintf(fid,'%s\n',evalc('display(err_stats_tree);'));
	fprintf(fid,'%s\n',evalc('display(err_stats_boosted);'));
fclose(fid);

figure(1); clf;
	plot(err_stats_tree.abs_percentiles,(1:100)/100,'b-',...
		 err_stats_boosted.abs_percentiles,(1:100)/100,'r-');
legend({'Tree','boosted'},'location','southeast');
graph(1); exportfig([figpath,'cum_freq',err_suffix]);

return

% plot betas (good and bad)
figure(2); clf;
	plot(abs(beta));
graph(2); exportfig([figpath,'betas_good']);
	plot(abs(beta_bad));
graph(2); exportfig([figpath,'betas_bad']);

% get histogram over split variables and plot
n_trees = 1;
hists = zeros(size(X_test,2),n_trees);
for t = 1:n_trees
    tree = u_load([forest.tree_root,forest.tree_dir,forest.trees{t}]);
	hists(:,t) = hist(tree.var(tree.var>0),1:60);
	hists(:,t) = hists(:,t)/sum(hists(:,t));
end
figure(2); clf;
	plot(hists);
graph(2); exportfig([figpath,'betas_treehists']);

profstat = profile('status');
if strcmp(profstat.ProfilerStatus,'on')
	profile off; profile report;
end
